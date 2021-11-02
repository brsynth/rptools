import csv
import libsbml
import rr_cache

import numpy as np
import pandas as pd

from copy import (
    copy as copy_copy
)
from logging import (
    Logger,
    getLogger
)
from os import (
    path as os_path
)
from typing import (
    Dict,
    List,
    Tuple
)
from cobra import (
    Model as cobra_model,
    Solution as cobra_solution
)
from cobra.io.sbml import _f_reaction
from cobra.medium import minimal_medium as cobra_minimal_medium

from rptools.rplibs.rpCompound import rpCompound
from rptools.rplibs.rpSBML import rpSBML

__CURRENT_PATH = os_path.dirname(
    os_path.abspath(__file__)
)
__MEDIUM_PATH = os_path.join(
    __CURRENT_PATH,
    'medium',
    'medium.csv'
)
__MEDIUM_DEFAULT_ID = 'not_predefined_model'

__MEDIUM_HEADER_NAME = 'medium_name'
__MEDIUM_HEADER_COMPOUND_ID = 'compound_id'
__MEDIUM_HEADER_BOUND = 'upper_bound'
__MEDIUM_HEADER_OPTIONAL = ['compound_annotation', 'compound_group', __MEDIUM_HEADER_NAME]
__MEDIUM_HEADER = __MEDIUM_HEADER_OPTIONAL + [__MEDIUM_HEADER_BOUND, __MEDIUM_HEADER_COMPOUND_ID]

#############
##  Utils  ##
#############
def is_df_medium_defined(
    df: pd.DataFrame=None
) -> bool:
    '''Verify if a dataframe is initialized.

    :param df: The Object to test

    :type df: pd.DataFrame

    :return: True if the dataframe exists and is not empty, False otherwise
    :rtype: bool
    '''
    if isinstance(df, pd.DataFrame): 
        if not df.empty:
            return True
    return False

##########
##  IO  ##
##########
class HeaderMalformated(Exception):
    """A simple Exception when header of a medium file is not properly set."""
    pass

def load_medium_file(
    filename:str,
    logger: Logger = getLogger(__name__)
) -> pd.DataFrame:    
    '''Load a file to fit in a dataframe supporting informations linked to a medium

    :param filename: Path of the file

    :type filename: str

    :return: A dataframe formatted
    :rtype: pd.DataFrame
    '''
 
    df = pd.DataFrame()
    try:
        if not os_path.isfile(filename):
            raise FileNotFoundError('File %s not found'%(filename,))
        # Find delimiter.
        with open(filename) as fid:
            dialect = csv.Sniffer().sniff(fid.readline())
        # Load.
        df = pd.read_csv(filename, sep=dialect.delimiter)
        # Check header.
        header_required = copy_copy(__MEDIUM_HEADER)
        for header in __MEDIUM_HEADER_OPTIONAL:
            if header in header_required:
                header_required.remove(header)
        for header in header_required:
            if not header in df.columns:
                raise HeaderMalformated('Header is malformated')        
        # Fmt.
        df[__MEDIUM_HEADER_BOUND] = df[__MEDIUM_HEADER_BOUND].astype(float)
        # Translate to rp_compound
        df = create_rp_compound(
            df, 
            logger
        )
    except Exception as e:
        logger.error(str(e))
        df = pd.DataFrame()
    return df

def read_medium_ids(
    filename: str,
    logger: Logger = getLogger(__name__)
) -> List[str]: 
    '''Read medium ids from a file describing a medium

    :param filename: Path of the file

    :type filename: str

    :return: Ids of interest
    :rtype: List[str]
    '''
 
    # Load medium file.
    df_medium = load_medium_file(filename)
    # Get medium ids.
    medium_ids = []
    if __MEDIUM_HEADER_NAME in df_medium.columns:
        medium_ids = df_medium[__MEDIUM_HEADER_NAME].unique().tolist()
    return medium_ids

def load_compounds(
    medium_id: str=__MEDIUM_DEFAULT_ID,
    file_base: str=__MEDIUM_PATH,
    file_user: str=None,
    logger: Logger = getLogger(__name__)
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    '''Load from files, information describing medium composition

    :param medium_id: An id of a medium to subset from file_base
    :param file_base: A file describing medium composition providing by rpfba
    :param file_user: A file describing medium composition provided by the user
    :param logger: A logging object to output information

    :type medium_id: str
    :type file_base: str
    :type file_user: str
    :type logger: Logger

    :return: Dataframes provided by rpfba or by the user
    :rtype: Tuple[pd.DataFrame, pd.DataFrame]
    '''
 
    # Load compounds - file base.
    df_medium_base = pd.DataFrame()
    if medium_id and medium_id != __MEDIUM_DEFAULT_ID:
        df_medium_base = load_medium_file(file_base, logger)
        if not df_medium_base.empty:
            df_medium_base = df_medium_base[df_medium_base[__MEDIUM_HEADER_NAME] == medium_id]
    # Load compounds - file user.
    df_medium_custom = pd.DataFrame()
    if file_user:
        df_medium_custom = load_medium_file(file_user, logger)

    return (df_medium_base, df_medium_custom)

##################
##  Annotation  ##
##################

def create_rp_compound(
    df: pd.DataFrame,
    logger: Logger=getLogger(__name__)
) -> pd.DataFrame:
    '''Associate an rpCompound object with each entry in a dataframe describing medium composition.
    Use rrCache to retrieve information about inchi, inchi_keys for a compound.

    :param df: a dataframe describing medium composition. Use a "compound_id" column.
    :param logger: A logging object to oyyutput information
    
    :type df: pd.DataFrame
    :type logger: Logger

    :return: a dataframe with an additional "rp_compound" column.
    :rtype: pd.DataFrame
    '''
 
    cache = rr_cache.rrCache(['cid_name'])
    def _create_rp_compound(x, cache):
        # Init.
        rp_compound = np.nan
        d_compound = {}
        # Search in rrcache.
        if not pd.isna(x[__MEDIUM_HEADER_COMPOUND_ID]) \
            and x[__MEDIUM_HEADER_COMPOUND_ID] in cache.get_list_of_compounds():
            d_compound = cache.get_compound(x[__MEDIUM_HEADER_COMPOUND_ID])
        # Fmt.
        if len(d_compound.keys()) > 0:
            if 'cid' in d_compound.keys():
                d_compound['id'] = d_compound.pop('cid')
            try:
                rp_compound = rpCompound(**d_compound)
            except BaseException as e:
                logger.error(str(e))
 
        return rp_compound
       
    df['rp_compound'] = df.apply(_create_rp_compound, axis=1, args=(cache,))
    return df

def crossref_medium_id(
    df: pd.DataFrame,
    model: rpSBML,
    compartment_id: str,
    logger: Logger=getLogger(__name__)
) -> pd.DataFrame:
    '''Perform a traduction between compounds of interest and compounds in the model.
    
    :param df: a dataframe describing medium composition. Use a "rp_compound" column which contains rpCompound objects.
    :param model: A model to retrieve compounds id.
    :param compartment_id: a compartment used to select compounds from the model to build traduction.

    :type df: pd.DataFrame
    :type model: rpSBML
    :type compartment_id: str

    :return: a dataframe with an additional "model_id" column.
    :rtype: pd.DataFrame
    '''
 
    if not is_df_medium_defined(df):
        return df

    # Get compound ids.
    for ix in df.index:
        if 'rp_compound' in df.columns and isinstance(df.loc[ix, 'rp_compound'], rpCompound):
            df.at[ix, 'tmp_compound_id'] = df.loc[ix, 'rp_compound'].get_id()
        elif __MEDIUM_HEADER_COMPOUND_ID in df.columns and not pd.isna(df.loc[ix, __MEDIUM_HEADER_COMPOUND_ID]):
            df.at[ix, 'tmp_compound_id'] = df.loc[ix, __MEDIUM_HEADER_COMPOUND_ID]
    compound_ids = [x for x in df['tmp_compound_id'] if not pd.isna(x)]

    # CrossRef componds.
    corr_species, miss_species = rpSBML.speciesMatchWith(
        compound_ids,
        model,
        compartment_id
    )
    if len(miss_species) > 0:
        logger.warning('These species provided for modify medium are not found in the model and they will not be included: %s'%(','.join(miss_species),))

    # Update id in df.
    def _crossref_medium_id(x, corr_species):
        compound_id = x['tmp_compound_id']
        if pd.isna(compound_id):
            return np.nan
        return corr_species.get(compound_id, np.nan)
    df['model_id'] = df.apply(_crossref_medium_id, axis=1, args=(corr_species,))
    
    # Clean
    df.drop('tmp_compound_id', axis=1, inplace=True, errors='ignore')

    return df


########################
##  Merging df - Fmt  ##
########################
def merge_medium(
    first: pd.DataFrame=None,
    second: pd.DataFrame=None
) -> pd.DataFrame:    
    '''Merge two dataframes describing medium composition.
    The second dataframe erase information of the first based on a column "model_id".

    :param first: a dataframe describing medium composition
    :param second: a dataframe describing medium composition

    :type first: pd.DataFrame
    :type second: pd.DataFrame

    :return: a dataframe resulting of merging of two dataframes
    :rtype: pd.DataFrame
    '''
 
    if first is None or second is None:
        return pd.DataFrame()
    df = first.append(second)
    if not df.empty and 'model_id' in df.columns:
        df.drop_duplicates('model_id', keep='last', inplace=True)
    return df

def merge_medium_exchange(
    medium: pd.DataFrame,
    exchange_reaction: pd.DataFrame
) -> pd.DataFrame:
    '''Merge two dataframes.
    The first one describes medium composition.
    The second one supports information about libsbml.Reaction and SpecieId associate for exchange reactions.
    The merge keeps all information.

    :param medium: a dataframe describing medium composition
    :param exchange_reaction: a dataframe describing exchange reaction

    :type first: pd.DataFrame
    :type second: pd.DataFrame

    :return: a dataframe resulting of merging of two dataframes
    :rtype: pd.DataFrame
    '''
 

    def _reaction_id_format(
        reaction: libsbml.Reaction,
        logger: Logger= getLogger(__name__)
    ) -> str:
        if pd.isna(reaction) or \
            not isinstance(reaction, libsbml.Reaction):
            return np.nan
        return _f_reaction(reaction.getIdAttribute())
    # Merge dfs
    df = pd.merge(medium, exchange_reaction, on='model_id', how='outer')
    # Fmt.
    df[__MEDIUM_HEADER_BOUND] = df[__MEDIUM_HEADER_BOUND].fillna(0.0)
    df['reaction_name'] = df['libsbml_reaction'].apply(_reaction_id_format)

    return df

def df_to_medium(
    df: pd.DataFrame
) -> Dict[str, float]:    
    '''Extract a dictionary from a dataframe, fitting with cobrapy needs.

    :param df: a dataframe with 'reaction_name', 'upper_bound', 'rp_compound', 'libsbml_reaction' columns

    :type df: pd.DataFrame

    :return: a dictionary with keys supporting 'reaction_name' and values supporting 'upper_bounds'
    :rtype: Dict[str, float]
    '''
    if not 'reaction_name' in df.columns \
        or not 'upper_bound' in df.columns \
        or not 'rp_compound' in df.columns \
        or not 'libsbml_reaction' in df.columns:
        return {}
    # Replace na.
    def _replace_na(x):
        if not pd.isna(x['rp_compound']) \
            and pd.isna(x['reaction_name']):
            return 'EX_%s'%(x['rp_compound'].get_id(),)
        return x['reaction_name']
    df['reaction_name'] = df.apply(_replace_na, axis=1)
    df = df[~pd.isna(df['reaction_name'])]
    # To dict.
    data = df[['reaction_name', __MEDIUM_HEADER_BOUND]].to_dict('records')
    data = {x['reaction_name']:x[__MEDIUM_HEADER_BOUND] for x in data}
    return data

###########################
##  Interact with Model  ##
###########################

def add_missing_specie(
    model: rpSBML,
    df: pd.DataFrame,
    compartment_id: str,
    logger: Logger = getLogger(__name__)
) -> rpSBML:
    '''Add species implied in an exchange reaction not present in a model

    :param model: a model into add specie
    :param df: a dataframe supporting information to add 
    :param compartment_id: id of the external compartment in the model
    :param logger: a logger object

    :type model: rpSBML
    :type df: pd.DataFrame
    :type compartment_id: str
    :type logger: Logger

    :return: the model modified
    :rtype: rpSBML
    '''

    # Add missing specie to the model.
    for ix in df.index:
        if pd.isna(df.loc[ix, 'libsbml_reaction']):
            rp_compound = df.loc[ix, 'rp_compound']
            if pd.isna(rp_compound):
                continue
            model.createSpecies(
                species_id=rp_compound.get_id(),
                species_name=rp_compound.get_name(),
                inchi=rp_compound.get_inchi(),
                inchikey=rp_compound.get_inchikey(),
                smiles=rp_compound.get_smiles(),
                compartment=compartment_id,
                is_boundary=True
            )
    return model

######################
##  Minimal Medium  ##
######################

def build_minimal_medium(
    model: cobra_model=None,
    solution: cobra_solution=None,
    logger: Logger=getLogger(__name__)
) -> pd.DataFrame:
    '''Finds the minimal growth medium for the `model` which allows for
    model as well as individual growth. Here, a minimal medium can either
    be the medium requiring the smallest total import flux and the medium
    requiring the least components (ergo ingredients).

    :param model: The model into perform search
    :param solution: The minimum growth rate (Optional: default as minimal_medium used by cobra)

    :type model: cobra.Model
    :type solution: cobra.Solution

    :return: Values for compounds which minimize the smallest total import flux and the minimial components
    :rtype: pd.DataFrame
    '''

    # Init.
    df = pd.DataFrame()
    if model is None:
        return data
    objective_value = 0.1
    if solution:
        objective_value = solution.objective_value
    # Minimize flux.
    df_flux = cobra_minimal_medium(
        model=model,
        min_objective_value=objective_value,
        exports=False,
        minimize_components=False,
        open_exchanges=False
    )
    # Minimize components.
    df_component = cobra_minimal_medium(
        model=model,
        min_objective_value=objective_value,
        exports=False,
        minimize_components=True,
        open_exchanges=False
    )
    # Concatenate results.
    if (not df_flux is None
        and isinstance(df_flux, pd.Series)):
        df = df_flux.to_frame('minimize_flux')
    if (not df_component is None
        and isinstance(df_component, pd.Series)):
        df_component = df_component.to_frame('minimize_component')
        df = pd.concat([df, df_component], axis=1)

    return df
