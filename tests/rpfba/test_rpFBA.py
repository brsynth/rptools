from glob import glob
from unittest import TestCase
from shutil import rmtree
from tempfile import (
    TemporaryDirectory,
    mkdtemp
)
from zipfile import ZipFile

from json import load as json_load
from os import path as os_path

from cobra.core.solution import Solution  as cobra_solution

from brs_utils import (
    create_logger,
    extract_gz
)
from rptools.rpfba.rpFBA import (
    check_SBML_compartment,
    check_SBML_rxnid,
    runFBA
)
from rptools.rplibs import (
    rpPathway,
    rpSBML
)
from main_rpfba import Main_rpfba

class Test_rpFBA(Main_rpfba):


    def setUp(self):
        super().setUp()
        input_zip = ZipFile(self.cr_path)
        input_zip.extractall(
            path=os_path.join(
                self.temp_d,
                'cr_fba'
            )
        )
        output_zip = ZipFile(self.fba_path)
        output_zip.extractall(
            path=os_path.join(
                self.temp_d,
                'lycopene_fba'
            )
        )
    def tearDown(self):
        super().tearDown()

    def test_runFBA(self):
        #TODO: taking accound extra args, like medium
        def _extract_var(dirname):
            files = glob(
                os_path.join(
                    dirname, 
                    '*xml'
                )
            )
            basenames = [os_path.basename(x) for x in files]
            names = [x.split('.')[0] for x in basenames]
            sims = [x.split('.')[1] for x in basenames]
            assert len(files) == len(names) == len(sims)
            return (files, names, sims)

        def _extract_res_from_file(filename):
            rp_pathway = rpPathway.from_rpSBML(
                infile=filename
            )
            res = {}
            res['pathway'] = rp_pathway.get_fba()
            res['reactions'] = {}
            for rid in rp_pathway.get_reactions():
                res['reactions'][rid] = rp_pathway.get_reaction(rid).get_fba()
            return res
       
        def _format_dict(old, new={}):
            print('old', old, 'new', new)
            for k, v in old.items():
                if isinstance(v, dict):
                    new[k] = _format_dict(old.get(k, {}), v)
                else:
                    if isinstance(v, str):
                        new[k] = v
                    else:
                        new[k] = round(float(v), 2)
            return new

        files, names, sims = _extract_var(os_path.join(self.temp_d, 'lycopene_fba'))

        for ix in range(len(files)):
            pathway_cr = rpPathway.from_rpSBML(
                infile=os_path.join(
                    self.temp_d,
                    'cr_fba',
                    names[ix]+'.xml'
                )
            )
            res = runFBA(
                pathway=pathway_cr,
                gem_sbml_path=self.e_coli_model_path,
                compartment_id='c',
                objective_rxn_id='rxn_target',
                biomass_rxn_id='biomass',
                sim_type=sims[ix]
            )

            res_previous = _format_dict(
                _extract_res_from_file(
                    files[ix]
                )
            )
            res_run_fba = _format_dict(
                {x:y for x,y in res.items() if x in ['pathway', 'reactions']}
            )

            self.assertDictEqual(
                res_previous,
                res_run_fba
            )

    def test_check_SBML_compartment(self):
        rpsbml = rpSBML(self.e_coli_model_path)
        # Return types.
        comp_id = 'cytosol'
        res = check_SBML_compartment(
            rpsbml=rpsbml,
            compartment_id=comp_id
        )
        self.assertIsInstance(
            res,
            str
        )
        # Values
        self.assertEqual(res, comp_id)
        # Challenge - 1
        comp_id = 'periplasm'
        res = check_SBML_compartment(
            rpsbml=rpsbml,
            compartment_id=comp_id
        )
        self.assertEqual(res, comp_id)
        # Challenge - 2
        comp_id = 'x'
        res = check_SBML_compartment(
            rpsbml=rpsbml,
            compartment_id=comp_id
        )
        self.assertIs(res, None)

    def test_check_SBML_rxnid(self):
        rpsbml = rpSBML(self.e_coli_model_path)
        # Return types.
        res = check_SBML_rxnid(
            rpsbml=rpsbml,
            rxn_id='biomass'
        )
        self.assertIsInstance(
            res,
            str
        )
        # Values
        self.assertEqual(
            res,
            'biomass'
        )
        # Challenge - 1
        res = check_SBML_rxnid(
            rpsbml=rpsbml,
            rxn_id='undefined'
        )
        self.assertIs(
            res,
            None
        )
 
#def build_rpsbml(
#def build_hidden_species(
#def build_results(
#def create_target_consumption_reaction(
#def write_results_to_pathway(
#def complete_heterologous_pathway(
#def rp_fraction(
#def runCobra(
#def build_cobra_model(
#def write_results_to_rpsbml(
#def write_fluxes_to_objectives(
#def write_fluxes_to_reactions(
#def write_objective_to_pathway(
#def create_ignored_species_group(
