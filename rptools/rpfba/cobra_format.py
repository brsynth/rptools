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

def cobra_suffix(pathway: rpPathway) -> str:
    return f'__64__{pathway.get_info("rpSBML")["compartments"][0]["id"]}'

def cobraize_string(
    string: str,
    pathway: rpPathway,
) -> str:
    return string+cobra_suffix(pathway)

def uncobraize_string(
    string: str,
    pathway: rpPathway,
) -> str:
    return string.replace(cobra_suffix(pathway), '')

def cobraize(pathway: rpPathway) -> None:
    '''Make the Pathway compliant with what Cobra expects
    Add <@compartmentID> to all compounds in species and reactions
    '''
    # SPECIES
    for spe_id in pathway.get_species_ids():
        pathway.rename_compound(spe_id, cobraize_string(spe_id, pathway))
        # pathway.get_specie(spe_id).set_id(cobraize_string(spe_id, pathway))

    # REACTIONS
    # cobraize_reactions(pathway)

def cobraize_reactions(pathway: rpPathway) -> None:
    '''Make the Pathway compliant with what Cobra expects
    Remove <@compartmentID> from all compounds in reactions
    '''
    for rxn_id in pathway.get_reactions_ids():
        rxn = pathway.get_reaction(rxn_id)
        reactants = deepcopy(rxn.get_reactants())
        rxn.set_reactants({})
        for spe_id, spe_sto in reactants.items():
            rxn.add_reactant(
                compound_id=cobraize_string(spe_id, pathway),
                stoichio=spe_sto
            )
        products = deepcopy(rxn.get_products())
        rxn.set_products({})
        for spe_id, spe_sto in products.items():
            rxn.add_product(
                compound_id=cobraize_string(spe_id, pathway),
                stoichio=spe_sto
            )
    print(Cache.get('RetroPath_Pathway_030_0001_rxn_1').to_dict())
    print(Cache.get('CMPD_0000000025__64__MNXC3').to_dict())
    print(Cache.get_list_of_objects())

def uncobraize(pathway: rpPathway) -> None:
    '''Make the Pathway compliant with what Cobra expects
    Remove <@compartmentID> from all compounds in species, reactions and scores
    '''
    for spe_id in pathway.get_species_ids():
        pathway.rename_compound(spe_id, uncobraize_string(spe_id, pathway))

def uncobraize_reactions(pathway: rpPathway) -> None:
    '''Make the Pathway compliant with what Cobra expects
    Remove <@compartmentID> from all compounds in reactions
    '''
    reactions = pathway.get_reactions()

    target_id = None
    for rxn in reactions:

        # REACTANTS
        reactants = rxn.get_reactants()
        rxn.set_reactants({})
        for spe_id, spe_sto in reactants:
            compound = pathway.get_specie(spe_id)
            compound.set_id(
                uncobraize_string(spe_id, pathway)
            )
            rxn.add_reactant(
                compound=compound,
                stoichio=spe_sto
            )

        # PRODUCTS
        products = rxn.get_products()
        rxn.set_products({})
        for spe_id, spe_sto in products:
            compound = pathway.get_specie(spe_id)
            compound.set_id(
                uncobraize_string(spe_id, pathway)
            )
            if spe_id == pathway.get_target_id():
                target_id = compound.get_id()
            rxn.add_product(
                compound=compound,
                stoichio=spe_sto
            )

        pathway.add_reaction(
            reaction=rxn,
            replace=True,
            target_id=target_id
        )
        if target_id is not None:
            target_id = None

def uncobraize_results(
  results: Dict,
  cobra_suffix: str,
  logger: Logger = getLogger(__name__)
) -> None:
  res = {
    'species': {},
    'reactions': {},
    'pathway': {}
  }
  logger.debug(cobra_suffix)
  logger.debug(results)
  # Uncobraize species results
  for spe_id, score in results['species'].items():
    res['species'][spe_id.replace(cobra_suffix, '')] = score
  # Write reactions results
  res['reactions'] = deepcopy(results['reactions'])
  # Write pathway result
  res['pathway'] = deepcopy(results['pathway'])
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
