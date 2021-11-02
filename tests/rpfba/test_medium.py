import csv
import libsbml
import rr_cache
import tempfile

import numpy as np
import pandas as pd

from copy import deepcopy
from logging import (
    Logger,
    getLogger
)
from os import (
    path as os_path,
    remove
)
from tempfile import NamedTemporaryFile
from typing import (
    Dict,
    Iterable,
)
from cobra import io as cobra_io
from cobra.io.sbml import _f_reaction
from cobra.medium.annotations import (
    excludes,
    sbo_terms
)
from rptools.rplibs.rpCompound import rpCompound
from rptools.rplibs.rpSBML import rpSBML
from rptools.rpfba import medium
from rptools.rpfba.medium import (
    build_minimal_medium,
    is_df_medium_defined,
    load_medium_file,
    read_medium_ids,
    load_compounds,
    create_rp_compound,
    crossref_medium_id,
    merge_medium,
    merge_medium_exchange,
    df_to_medium,
    add_missing_specie
)

from main_rpfba import Main_rpfba

class Test_medium(Main_rpfba):

    # TODO: import directly from module
    __MEDIUM_DEFAULT_ID = 'not_predefined_model'
    __MEDIUM_HEADER_NAME = 'medium_name'
    __MEDIUM_HEADER_COMPOUND_ID = 'compound_id'
    __MEDIUM_HEADER_BOUND = 'upper_bound'
    __MEDIUM_HEADER_OPTIONAL = ['compound_annotation', 'compound_group']
    __MEDIUM_HEADER = __MEDIUM_HEADER_OPTIONAL + [__MEDIUM_HEADER_BOUND, __MEDIUM_HEADER_COMPOUND_ID, __MEDIUM_HEADER_NAME]

    def load_medium_file(filename):
        medium = pd.read_csv(filename)
        return create_rp_compound(
            df=medium,
            logger=self.logger
        )
        
    def setUp(self):
        super().setUp()
        # self.logger.setLevel('DEBUG')
        # objects below have to be created for each test instance
        # since some tests can modified them

    def test_is_df_medium_defined(self):
        # Return type
        self.assertTrue(isinstance(is_df_medium_defined(None), bool))
        # Values
        self.assertFalse(is_df_medium_defined(None))
        self.assertFalse(is_df_medium_defined(np.nan))
        self.assertFalse(is_df_medium_defined(pd.DataFrame()))
        self.assertFalse(is_df_medium_defined(pd.DataFrame(columns=['a'])))
        self.assertFalse(is_df_medium_defined(pd.DataFrame(index=[0])))
        self.assertTrue(is_df_medium_defined(pd.DataFrame(data=[1], columns=['a'], index=[0])))

    def test_load_medium_file(self):
        df = load_medium_file(os_path.join(self.medium_path, 'medium.io.a.tsv'))
        # Return type
        self.assertTrue(isinstance(df, pd.DataFrame))
        # Basic io profile
        self.assertTrue(is_df_medium_defined(df))
        df = load_medium_file(os_path.join(self.medium_path, 'medium.io.d.csv'))
        self.assertFalse(is_df_medium_defined(df))
        df = load_medium_file(os_path.join(self.medium_path, 'medium.io.a.xlsx'))
        self.assertFalse(is_df_medium_defined(df))
        df = load_medium_file(os_path.join(self.medium_path, 'medium.io.a.csv'))
        self.assertTrue(is_df_medium_defined(df))
        # Type expected
        self.assertTrue(pd.api.types.is_float_dtype(df[self.__MEDIUM_HEADER_BOUND]))
        self.assertEqual(
            sum(
                df['rp_compound'].apply(lambda x: isinstance(x, rpCompound) or pd.isna(x))
            ),
            len(df['rp_compound'])
        )
        # Challenge on column labels
        df_columns = df.columns.tolist()
        df_columns.remove('rp_compound')
        self.assertEqual(
            sorted(df_columns),
            sorted(self.__MEDIUM_HEADER)
        )
        tmp_file = tempfile.NamedTemporaryFile(
                suffix='.csv', 
                dir=self.temp_d,
                delete=False
        )
        for ix in range(len(self.__MEDIUM_HEADER)):
            tmp_header = deepcopy(self.__MEDIUM_HEADER)
            tmp_header = tmp_header.pop(ix)
            df_tmp = df[tmp_header]
            df_tmp.to_csv(
                tmp_file.name, 
                index=False
            )
            df_tmp = load_medium_file(tmp_file.name)
            self.assertFalse(is_df_medium_defined(df_tmp))
            tmp_file.close()
            remove(tmp_file.name)

    def test_read_medium_ids(self):
        ids = read_medium_ids(os_path.join(self.medium_path, 'medium.io.b.csv'))
        # Return type.
        self.assertTrue(ids, Iterable)
        # Values.
        self.assertEqual(
            sorted([x for x in ids if not pd.isna(x)]),
            sorted(['m9', 'lb', 'lc'])
        )

    def test_load_compounds(self):
        df = pd.read_csv(os_path.join(self.medium_path, 'medium.io.c.csv'))
        # Return type.
        dfs = load_compounds(
            self.__MEDIUM_DEFAULT_ID,
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertTrue(isinstance(dfs, Iterable))
        self.assertEqual(len(dfs), 2)
        self.assertTrue(isinstance(dfs[0], pd.DataFrame))
        self.assertTrue(isinstance(dfs[1], pd.DataFrame))
        # Base.
        df_base, df_user = load_compounds(
            self.__MEDIUM_DEFAULT_ID,
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertEqual(df_base.shape, (0,0))
        self.assertEqual(df_user.shape[0], df.shape[0])
        self.assertEqual(df_user.shape[1]-1, df.shape[1])
        # Select by id
        df_base, df_user = load_compounds(
            'm9',
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertEqual(df_base.shape[0], 2)
        self.assertEqual(df_base.shape[1]-1, len(self.__MEDIUM_HEADER))
        self.assertEqual(df_user.shape[0], df.shape[0])
        self.assertEqual(df_user.shape[1]-1, df.shape[1])
 
        # Challenge
        df_base, df_user = load_compounds(
            '',
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertFalse(is_df_medium_defined(df_base))
        df_base, df_user = load_compounds(
            np.nan,
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertFalse(is_df_medium_defined(df_base))
        df_base, df_user = load_compounds(
            self.__MEDIUM_DEFAULT_ID,
            self.temp_d,
            os_path.join(self.medium_path, 'medium.io.c.csv')
        )
        self.assertFalse(is_df_medium_defined(df_base))
        self.assertTrue(is_df_medium_defined(df_user))
        df_base, df_user = load_compounds(
            self.__MEDIUM_DEFAULT_ID,
            os_path.join(self.medium_path, 'medium.io.c.csv'),
            self.temp_d
        )
        self.assertFalse(is_df_medium_defined(df_base))
        self.assertFalse(is_df_medium_defined(df_user))

    def test_create_rp_compound(self):
        df = pd.read_csv(os_path.join(self.medium_path, 'medium.annotation.a.csv'))
        df = create_rp_compound(
            df=df,
            logger=self.logger
        )
        # Return type.
        self.assertTrue(isinstance(df, pd.DataFrame))
        # Values.
        self.assertTrue(isinstance(df.loc[0, 'rp_compound'], rpCompound))
        self.assertTrue(isinstance(df.loc[1, 'rp_compound'], rpCompound))
        self.assertTrue(isinstance(df.loc[2, 'rp_compound'], rpCompound))
        self.assertEqual(
            df.loc[0, 'rp_compound'].get_id(), 
            df.loc[2, 'rp_compound'].get_id()
        )
        self.assertTrue(pd.isna(df.loc[3, 'rp_compound']))
        self.assertTrue(pd.isna(df.loc[4, 'rp_compound']))

    def test_crossref_medium_id(self):
        # Load.
        medium = load_medium_file(os_path.join(self.medium_path, 'medium.annotation.b.csv'))
        # Return type.
        df = crossref_medium_id(
            df=None,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        self.assertEqual(df, None)
        df = crossref_medium_id(
            df=medium,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        self.assertTrue(isinstance(df, pd.DataFrame))
        # Values.
        self.assertEqual(
            ['M_pi_e', 'M_fe3_e', 'M_pi_e', 'M_pi_e'],
            df['model_id'].tolist()[:4]
        )
        self.assertEqual(
            df['model_id'].isna().sum(),
            2
        )
        df = crossref_medium_id(
            df=medium,
            model=self.rpsbml,
            compartment_id='MNXC'
        )
        self.assertEqual(
            df['model_id'].isna().sum(),
            6
        )
    ###########
    ##  Fmt  ##
    ###########
    def test_merge_medium(self):
        medium_a = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.a.csv'))
        medium_b = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.b.csv'))
        dfa = crossref_medium_id(
            df=medium_a,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        dfb = crossref_medium_id(
            df=medium_b,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        df = merge_medium(first=None, second=dfa)
        self.assertIsInstance(df, pd.DataFrame)
        df = merge_medium(first=dfa, second=None)
        self.assertIsInstance(df, pd.DataFrame)
        df = merge_medium(first=dfa, second=pd.DataFrame())
        self.assertEqual(
            df.shape, 
            dfa.shape
        )
        df = merge_medium(first=pd.DataFrame(), second=dfa)
        self.assertEqual(
            df.shape, 
            dfa.shape
        )
        # Values
        df = merge_medium(dfa, dfb)
        self.assertIsInstance(df, pd.DataFrame) # Return type
        self.assertEqual(
            df.shape[0],
            3
        )
        self.assertNotIn('MNXM9', df[self.__MEDIUM_HEADER_COMPOUND_ID])
        df = merge_medium(dfb, dfa)
        self.assertEqual(
            df.shape[0],
            3
        )
        self.assertNotIn('14791', df[self.__MEDIUM_HEADER_COMPOUND_ID])

    def test_merge_medium_exchange(self):
        # Load
        medium = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.a.csv'))
        medium = crossref_medium_id(
            df=medium,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        exchange = self.rpsbml.build_exchange_reaction('c')
        df = merge_medium_exchange(
            medium = medium,
            exchange_reaction = exchange
        )
        # Return type.
        self.assertIsInstance(
            df,
            pd.DataFrame
        )
        # Values.
        self.assertEqual(
            df.shape,
            (331, len(exchange.columns)+len(df.columns)-2)
        )
        self.assertEqual(
            df[self.__MEDIUM_HEADER_COMPOUND_ID].isna().sum(),
            331-2
        )
        self.assertTrue(pd.api.types.is_float_dtype(df[self.__MEDIUM_HEADER_BOUND]))
 
    def test_df_to_medium(self):
        # Load.
        medium_a = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.a.csv'))
        dfa = crossref_medium_id(
            df=medium_a,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        exchange = self.rpsbml.build_exchange_reaction('c')
        dfa = merge_medium_exchange(
            medium = dfa,
            exchange_reaction = exchange
        )
 
        # Return Type.
        dfad = dfa.drop(columns=['reaction_name'])
        medium = df_to_medium(dfad)
        self.assertIsInstance(medium, Dict)
        dfad = dfa.drop(columns=['reaction_name'])
        medium = df_to_medium(dfad)
        self.assertIsInstance(medium, Dict)
        medium = df_to_medium(dfa)
        self.assertIsInstance(medium, Dict)
        # Values.
        self.assertEqual(
            sum([x.startswith('EX') for x in medium.keys()]),
            dfa.shape[0]
        )

    ################
    ##  Interact  ##
    ################

    def test_add_missing_specie(self):
        # Load regular data
        specie_missing_id = 'MNXM83'
        medium = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.c.csv'))
        medium = crossref_medium_id(
            df=medium,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        exchange = self.rpsbml.build_exchange_reaction('c')
        df = merge_medium_exchange(
            medium = medium,
            exchange_reaction = exchange
        )
        rpsbml = add_missing_specie(
            self.rpsbml,
            df,
            'c'
        )
        # Return type
        self.assertIsInstance(rpsbml, rpSBML)
        # Values
        species = rpsbml.getModel().getListOfSpecies()
        self.assertEqual(
            len([x for x in species if x.getId() == specie_missing_id]),
            1
        )
        specie = rpsbml.getModel().getSpecies(specie_missing_id)
        self.assertTrue(specie.getBoundaryCondition())

    def test_global_eff(self):
        # Load regular data
        expected_medium = {
            'EX_MNXM83': 2000.0, 
            'EX_pi_e': 2000.0, 
            'EX_fe3_e': 1000.0
        }
        specie_missing_id = 'MNXM83'
        medium = load_medium_file(os_path.join(self.medium_path, 'medium.fmt.c.csv'))
        medium = crossref_medium_id(
            df=medium,
            model=self.rpsbml,
            compartment_id='MNXC2'
        )
        exchange = self.rpsbml.build_exchange_reaction('c')
        df = merge_medium_exchange(
            medium = medium,
            exchange_reaction = exchange
        )
        rpsbml = add_missing_specie(
            self.rpsbml,
            df,
            'c'
        )
        cobra_model = None
        with NamedTemporaryFile(dir=self.temp_d, delete=False) as tempf:
            rpsbml.write_to_file(tempf.name)
            cobra_model=cobra_io.read_sbml_model(tempf.name, use_fbc_package=True)
        cobra_model.medium = df_to_medium(df)
        self.assertEqual(
            cobra_model.medium,
            expected_medium
        )
        tempf.close()
        remove(tempf.name)

    def test_build_minimal_medium(self):
        # Load regular data
        cobra_model = None
        with NamedTemporaryFile(dir=self.temp_d, delete=False) as tempf:
            self.rpsbml.write_to_file(tempf.name)
            cobra_model=cobra_io.read_sbml_model(tempf.name, use_fbc_package=True)
        cobra_solution = cobra_model.optimize()

        df1 = build_minimal_medium(
            model=cobra_model,
            solution=cobra_solution
        )
        df2 = build_minimal_medium(
            model=cobra_model
        )
        # Return values.
        self.assertIsInstance(
            df1,
            pd.DataFrame
        )
        # Values
        self.assertEqual(
            df2.shape[1],
            2
        )
        # Close.
        tempf.close()
        remove(tempf.name)
