from logging import (
    Logger,
    getLogger
)
from typing import (
    Dict
)
from copy import deepcopy
from rptools.rplibs import rpPathway
from brs_utils import Cache

at_pattern = '__64__'
bigg_prefix = 'M_'

def to_cobra(string: str) -> str:
    if string.startswith(bigg_prefix):
        string = string[len(bigg_prefix):]
    return string.replace(at_pattern, '@')

def cobra_suffix(comp_id: str) -> str:
    return at_pattern+comp_id

def cobraize(
    string: str,
    comp_id: str
) -> str:
    suffix = cobra_suffix(comp_id)
    if string.endswith(suffix):
        return string
    else:
        return string+cobra_suffix(comp_id)

def uncobraize(
    string: str
) -> str:
    return string.split(at_pattern)[0]

# def cobraize_reactions(pathway: rpPathway) -> None:
#     '''Make the Pathway compliant with what Cobra expects
#     Remove <@compartmentID> from all compounds in reactions
#     '''
#     for rxn_id in pathway.get_reactions_ids():
#         rxn = pathway.get_reaction(rxn_id)
#         reactants = deepcopy(rxn.get_reactants())
#         rxn.set_reactants({})
#         for spe_id, spe_sto in reactants.items():
#             rxn.add_reactant(
#                 compound_id=cobraize_string(spe_id, pathway),
#                 stoichio=spe_sto
#             )
#         products = deepcopy(rxn.get_products())
#         rxn.set_products({})
#         for spe_id, spe_sto in products.items():
#             rxn.add_product(
#                 compound_id=cobraize_string(spe_id, pathway),
#                 stoichio=spe_sto
#             )
#     print(Cache.get('RetroPath_Pathway_030_0001_rxn_1').to_dict())
#     print(Cache.get('CMPD_0000000025__64__MNXC3').to_dict())
#     print(Cache.get_list_of_objects())

# def uncobraize_reactions(pathway: rpPathway) -> None:
#     '''Make the Pathway compliant with what Cobra expects
#     Remove <@compartmentID> from all compounds in reactions
#     '''
#     reactions = pathway.get_reactions()

#     target_id = None
#     for rxn in reactions:

#         # REACTANTS
#         reactants = rxn.get_reactants()
#         rxn.set_reactants({})
#         for spe_id, spe_sto in reactants:
#             compound = pathway.get_specie(spe_id)
#             compound.set_id(
#                 uncobraize_string(spe_id, pathway)
#             )
#             rxn.add_reactant(
#                 compound=compound,
#                 stoichio=spe_sto
#             )

#         # PRODUCTS
#         products = rxn.get_products()
#         rxn.set_products({})
#         for spe_id, spe_sto in products:
#             compound = pathway.get_specie(spe_id)
#             compound.set_id(
#                 uncobraize_string(spe_id, pathway)
#             )
#             if spe_id == pathway.get_target_id():
#                 target_id = compound.get_id()
#             rxn.add_product(
#                 compound=compound,
#                 stoichio=spe_sto
#             )

#         pathway.add_reaction(
#             reaction=rxn,
#             replace=True,
#             target_id=target_id
#         )
#         if target_id is not None:
#             target_id = None

def uncobraize_results(
    results: Dict,
    cobra_suffix: str,
    logger: Logger = getLogger(__name__)
) -> None:

  res = {
    'species': {},
    'reactions': {},
    'pathway': {},
    'rpfba_ignored_species': []
  }
  logger.debug(cobra_suffix)
  logger.debug(results)
  # Uncobraize species results
  for spe_id, score in results['species'].items():
    res['species'][spe_id.replace(cobra_suffix, '')] = score
  # Uncobraize rpfba_ignored_species results
  for spe_id in results['rpfba_ignored_species']:
    res['rpfba_ignored_species'] += [spe_id.replace(cobra_suffix, '')]
  # Copy other results
  for key in ['reactions', 'pathway']:
      res[key] = deepcopy(results[key])
  return res



# def uncobraize_scores(pathway: Pathway) -> None:
#     '''Make the Pathway compliant with what Cobra expects
#     Remove <@compartmentID> from all compounds in scores
#     '''
#     uncobraize_scores_fba(pathway)

# def uncobraize_scores_fba(pathway: Pathway) -> None:
#     '''Make the Pathway compliant with what Cobra expects
#     Remove <@compartmentID> from all compounds in FBA scores
#     '''
#     # SPECIES
#     species = pathway.get_fba()['species'].keys()
#     for spe_id in species:
#         print(spe_id, pathway.get_fba_species(spe_id))
#         for objective_id, value in pathway.get_fba_species(spe_id).items():
#             pathway.set_fba_species(
#                 spe_id=uncobraize_string(spe_id, pathway),
#                 objective_id=objective_id,
#                 value=value
#             )
#     exit()
