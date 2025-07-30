#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Utility methods to build rpviz html pages."""

__author__ = 'Thomas Duigou'
__license__ = 'MIT'

import os
import csv
import logging
from typing import Dict, Union

from rptools.rplibs import rpSBML, rpPathway
from rptools.rplibs.rpReaction import rpReaction
from rptools.rplibs.rpCompound import rpCompound
from rptools.rpfba.cobra_format import uncobraize

DEBUG = True

miriam_header = {
    'compartment': {
        'go': 'go/GO:',
        'mnx': 'metanetx.compartment/',
        'bigg': 'bigg.compartment/',
        'seed': 'seed/',
        'name': 'name/'
    },
    'reaction': {
        'metanetx': 'metanetx.reaction/',
        'rhea': 'rhea/',
        'reactome': 'reactome/',
        'bigg': 'bigg.reaction/',
        'sabiork': 'sabiork.reaction/',
        'ec-code': 'ec-code/',
        'biocyc': 'biocyc/', 'lipidmaps': 'lipidmaps/'
    },
    'species': {
        'metanetx': 'metanetx.chemical/',
        'chebi': 'chebi/CHEBI:',
        'bigg': 'bigg.metabolite/',
        'hmdb': 'hmdb/',
        'kegg_c': 'kegg.compound/',
        'kegg_d': 'kegg.drug/',
        'biocyc': 'biocyc/META:',
        'seed': 'seed.compound/',
        'metacyc': 'metacyc.compound/',
        'sabiork': 'sabiork.compound/',
        'reactome': 'reactome/R-ALL-'
    }
}


def _specie_is_target_old(specie_id: str) -> bool:
    """Detect is a specie should be considered as a target

    FIXME, this needs to be refined so that we don't rely on the specie ID.

    :param specie_id: specie ID
    :type: str
    :return: true if it is a target, otherwise false
    :rtype: bool
    """
    if specie_id.startswith('TARGET_'):
        return True
    return False


def _specie_is_intermediate_old(
    specie_id: str,
    specie_dict: dict = None,
) -> bool:
    """Detect is a specie should be considered as an intermediate compound.

    FIXME, this needs to be refined so that we don't rely on the specie ID.

    :param specie_id: specie ID
    :type: str
    :param specie_dict: dictionary about the specie
    :type specie_dict: dict
    :return: true if it is, otherwise false
    :rtype: bool
    """
    if specie_id.startswith('CMPD_'):
        return True
    return False


def _specie_is_sink_old(
    specie_id: str,
    specie_dict: dict = None,
) -> bool:
    """Detect is a specie should be considered as a sink

    FIXME, this needs to be refined so that we don't rely on the specie ID.

    :param specie_id: specie ID
    :type: str
    :param specie_dict: dictionary about the specie
    :type specie_dict: dict
    :return: true if it is, otherwise false
    :rtype: bool
    """
    if (
        not _specie_is_target_old(specie_id)
        and not _specie_is_intermediate_old(specie_id)
    ):
        return True
    return False


def _get_pathway_score(rp_pathway: rpPathway) -> dict:
    # Precompute rule score
    rscores = []
    for rxn in rp_pathway.get_reactions_ids():
        rscores.append(rp_pathway.get_reaction(rxn).get_rule_score())
    rscore = sum(rscores) / len(rscores)
    #
    scores = {
        'rule_score': rscore,
        'steps': rp_pathway.get_nb_reactions(),
        'global_score': _get_pathway_global_score(rp_pathway)
    }
    try:
        scores['thermo_dg_m_gibbs'] = rp_pathway.get_thermo_dGm_prime()['value']
    except TypeError:
        pass
    try:
        scores['fba_target_flux'] = rp_pathway.get_fba_fraction()['value']
    except TypeError:
        pass

    return scores


def _get_pathway_global_score(rp_pathway: rpPathway) -> dict:
    if rp_pathway.get_global_score() == -1:
        return None
    return rp_pathway.get_global_score()


def _get_pathway_thermo(pathway_dict: dict) -> Union[float, None]:
    try:
        return pathway_dict['brsynth']['dfg_prime_m']['value']
    except KeyError:
        return None


def _get_pathway_fba(pathway_dict: dict) -> Union[float, None]:
    try:
        return pathway_dict['brsynth']['fba_obj_fraction']['value']
    except KeyError:
        return None


def _get_reaction_node_id_old(rxn_dict: dict) -> str:
    """Return a useful ID for the reaction node.

    A reaction node could be shared between several pathways, the reaction
    SMILES is an easy way to detect identical reactions used by different
    pathways.
    """
    if _get_reaction_smiles_old(rxn_dict) is not None:
        return rxn_dict['brsynth']['smiles']
    else:
        raise NotImplementedError(
            f'Cannot assign a valid ID to reaction idx {rxn_dict["rxn_idx"]} '
        )


def _get_reaction_node_id(rxn: rpReaction) -> str:
    """Return a useful ID for the reaction node.

    A reaction node could be shared between several pathways, the reaction
    SMILES is an easy way to detect identical reactions used by different
    pathways.
    """
    if _rxn_has_smiles(rxn):
        return rxn.get_smiles()
    else:
        raise NotImplementedError(
            f'Cannot assign a valid ID to reaction idx {rxn.get_id()} '
        )


def _rxn_has_smiles(rxn: rpReaction) -> bool:
    if (rxn.get_smiles() is None) or \
        (rxn.get_smiles() == '>>') or \
        (rxn.get_smiles() == ''):
        return False
    return True


def _get_reaction_ecs(rxn_dict: dict) -> list:
    if 'ec-code' in rxn_dict['miriam'] \
            and len(rxn_dict['miriam']['ec-code']):
        return rxn_dict['miriam']['ec-code']
    else:
        return []


def _get_reaction_thermo(rxn_dict: dict) -> Union[float, None]:
    if 'dfG_prime_m' in rxn_dict['brsynth']:
        return rxn_dict['brsynth']['dfG_prime_m']
    else:
        return None


def _get_reaction_labels_old(rxn_dict: dict) -> list:
    if len(_get_reaction_ecs(rxn_dict)):
        return [*_get_reaction_ecs(rxn_dict)]
    else:
        return [rxn_dict['brsynth']['rule_ids'],]


def _get_reaction_labels(rxn: rpReaction) -> list:
    if len(rxn.get_ec_numbers()):
        return rxn.get_ec_numbers()
    elif len(rxn.get_tmpl_rxn_ids()):
        return rxn.get_tmpl_rxn_ids()
    else:
        return [rxn.get_rule_id(),]


def _get_reaction_smiles_old(rxn_dict: dict) -> Union[str, None]:
    if 'smiles' in rxn_dict['brsynth'] \
            and rxn_dict['brsynth']['smiles'] is not None \
            and rxn_dict['brsynth']['smiles'] != '':
        return rxn_dict['brsynth']['smiles']
    else:
        return None


def _get_reaction_xlinks_old(rxn_dict: dict) -> list:
    # TODO refine this method
    xlinks = []
    # Special case for EC numbers
    for ec in _get_reaction_ecs(rxn_dict):
        # Get rid of unwanted characters
        ec_tmp = []
        for digit in ec.split('.'):
            if digit in '-_' or digit == '':
                break
            ec_tmp.append(digit)
        ec_refined = '.'.join(ec_tmp)
        if ec != ec_refined:
            logging.info(f'Refined EC number from {ec} to {ec_refined}')
        # Use direct link to workaround generic ECs issue with identifiers.org
        xlinks.append({
            'db_name': 'intenz',
            'entity_id': ec_refined,
            'url': f'https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec={ec_refined}'})
        logging.debug(
            f'Replace identifiers.org to IntEnz crosslinks for EC number {ec_refined}')
    # Not EC cases
    # TODO complete me
    return xlinks


def _get_reaction_xlinks(rxn: rpReaction) -> list:
    xlinks = []
    for ec in rxn.get_ec_numbers():
        # Get rid of unwanted characters
        ec_tmp = []
        for digit in ec.split('.'):
            if digit in '-_' or digit == '':
                break
            ec_tmp.append(digit)
        ec_refined = '.'.join(ec_tmp)
        if ec != ec_refined:
            logging.info(f'Refined EC number from {ec} to {ec_refined}')
        # Use direct link to workaround generic ECs issue with identifiers.org
        xlinks.append({
            'db_name': 'intenz',
            'entity_id': ec_refined,
            'url': f'https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec={ec_refined}'})
        logging.debug(
            f'Replace identifiers.org to IntEnz crosslinks for EC number {ec_refined}')
    # Not EC cases
    # TODO complete me
    return xlinks


def _get_reaction_rule_score(rxn_dict: dict) -> Union[float, None]:
    try:
        return round(rxn_dict['brsynth']['rule_score'], 3)
    except KeyError:
        return None


def _get_specie_node_id_old(specie_dict: dict, specie_id: str = None) -> str:
    """Return a useful ID for the specie node.

    A compound/specie node could be shared between several pathways,
    the inchikey, MNXM ID and chebi ID are reliable way to detect
    identical compounds.
    """
    if _get_specie_inchikey(specie_dict) is not None:
        return _get_specie_inchikey(specie_dict)
    elif 'metanetx' in specie_dict['miriam'] \
            and len(specie_dict['miriam']['metanetx']):
        sorted_mnx_ids = sorted(
            specie_dict['miriam']['metanetx'],
            key=lambda x: int(x.replace('MNXM', ''))
        )
        return sorted_mnx_ids[0]
    elif 'chebi' in specie_dict['miriam'] \
            and len(specie_dict['miriam']['chebi']):
        sorted_chebi_ids = sorted(
            specie_dict['miriam']['chebi'],
            key=lambda x: int(x.replace('CHEBI:', ''))
        )
        return sorted_chebi_ids[0]
    elif DEBUG and specie_id is not None:
        logging.error('Workaround for node ID, this should be fixed')
        return specie_id
    else:
        raise NotImplementedError('Could not assign a valid id')


def _get_specie_inchikey(specie_dict: dict) -> Union[str, None]:
    try:
        return specie_dict['brsynth']['inchikey']
    except KeyError:
        return None


def _get_specie_smiles(specie_dict: dict) -> Union[str, None]:
    try:
        return specie_dict['brsynth']['smiles']
    except KeyError:
        return None


def _get_specie_inchi(specie_dict: dict) -> Union[str, None]:
    try:
        return specie_dict['brsynth']['inchi']
    except KeyError:
        return None


def _get_specie_xlinks_old(specie_dict: dict) -> list:
    _MIRIAM_TO_IDENTIFIERS = {
        'metanetx': 'metanetx.chemical/',
        'chebi': 'chebi/CHEBI:',
        'bigg': 'bigg.metabolite/',
        'hmdb': 'hmdb/',
        'kegg_c': 'kegg.compound/',
        'kegg_d': 'kegg.drug/',
        'biocyc': 'biocyc/META:',
        'seed': 'seed.compound/',
        'metacyc': 'metacyc.compound/',
        'sabiork': 'sabiork.compound/',
        'reactome': 'reactome/'
    }
    xlinks = []
    if 'miriam' in specie_dict:
        for db_name, db_id_list in specie_dict['miriam'].items():
            for db_id in db_id_list:
                # Direct link to metacyc because identifiers.org is bugged
                if db_name == 'metacyc':
                    url_str = f'https://metacyc.org/compound?id={db_id}'
                # KEGG cases
                elif db_name == 'kegg' and db_id[0] == 'C':
                    url_str = f'http://identifiers.org/{_MIRIAM_TO_IDENTIFIERS["kegg_c"]}{db_id}'
                elif db_name == 'kegg' and db_id[0] == 'D':
                    url_str = f'http://identifiers.org/{_MIRIAM_TO_IDENTIFIERS["kegg_d"]}{db_id}'
                else:
                    url_str = f'http://identifiers.org/{_MIRIAM_TO_IDENTIFIERS[db_name]}{db_id}'
                xlinks.append({
                    'db_name': db_name,
                    'entity_id': db_id,
                    'url': url_str
                })
    return xlinks


def _get_specie_xlinks(cmpd: rpCompound) -> dict:
    return [{
        'db_name': 'metanetx',
        'entity_id': uncobraize(cmpd.get_id()),
        'url': f'http://identifiers.org/metanetx.chemical/{uncobraize(cmpd.get_id())}'
    }]


def _nodes_seem_equal(node1: dict, node2: dict) -> bool:
    # Few basic checks
    if (
        node1['id'] == node2['id']
        and node1['type'] == node2['type']
        and node1['label'] == node2['label']
        and node1['smiles'] == node2['smiles']
        and node1['inchi'] == node2['inchi']
        and node1['inchikey'] == node2['inchikey']
        and node1['rsmiles'] == node2['rsmiles']
        and node1['rule_ids'] == node2['rule_ids']
    ):
        return True
    return False


def _edge_seem_equal(edge1: dict, edge2: dict) -> bool:
    if (
        edge1['id'] == edge2['id']
        and edge1['source'] == edge2['source']
        and edge1['target'] == edge2['target']
    ):
        return True
    return False


def _merge_nodes(node1: dict, node2: dict) -> dict:
    node3 = {}
    for key in node1.keys():
        # Only node 1 has a value
        if node1[key] is None:
            value = node2[key]
        # Only node 2 has a value
        elif node2[key] is None:
            value = node1[key]
        # Both have a value
        else:
            # list of strings
            if key in ['path_ids', 'rule_ids', 'all_labels']:
                value = list(set(node1[key] + node2[key]))
            # important str values
            elif key in ['smiles', 'inchi', 'inchikey']:                
                value = node1[key]
                if node1[key] != node2[key]:
                    logging.warning(
                        f'Not the same {key} when merging nodes: '
                        f'{node1[key]} vs {node2[key]}. '
                        f'Keeping the first one'
                    )
            # float
            elif key == 'rule_score':  # float value
                value = max(node1[key], node2[key])
            # list of dicts
            elif key == 'xlinks':
                value = []
                done = set()
                for entry in node1[key]:
                    tag = f'{entry["db_name"]}-{entry["entity_id"]}'
                    if tag not in done:
                        value.append(entry)
                        done.add(tag)
            # backup plan
            else:
                value = node1[key]
        node3[key] = value

    return node3


def _merge_edges(edge1: dict, edge2: dict) -> dict:
    edge3 = {}
    for key in edge1.keys():
        if key == 'path_ids':
            value = sorted(list(set(edge1[key] + edge2[key])))
        else:
            value = edge1[key]
        edge3[key] = value
    return edge3


def parse_one_pathway(rp_pathway: rpPathway) -> tuple:
    """Extract info from one rpSBML file

    :param sbml_path: str, path to file
    """
    nodes = {}
    edges = {}

    pathway = {
        'path_id': rp_pathway.get_id(),
        'nb_steps': rp_pathway.get_nb_reactions(),
        'node_ids': [],  # To be filled later
        'edge_ids': [],  # To be filled later
        'scores': _get_pathway_score(rp_pathway)
    }

    # Node info: reactions
    for rxn in rp_pathway.get_reactions().values():
        node = {
            'id': _get_reaction_node_id(rxn),
            'path_ids': [pathway['path_id'], ],
            'type': 'reaction',
            'label': _get_reaction_labels(rxn)[0],
            'all_labels': _get_reaction_labels(rxn),
            'svg': None,  # FIXME could add the reaction depiction here
            'xlinks': _get_reaction_xlinks(rxn),
            # Only for reaction, None for compounds
            'rsmiles': rxn.get_smiles(),
            'rule_ids': rxn.get_rule_ids(),
            'rxn_template_ids': rxn.get_tmpl_rxn_ids(),
            'ec_numbers': rxn.get_ec_numbers(),
            'thermo_dg_m_gibbs': rxn.get_thermo_dGm_prime(),
            'rule_score': rxn.get_rule_score(),
            'uniprot_ids': rxn.get_selenzy_infos(),
            # Only for compounds
            'smiles': None,
            'inchi': None,
            'inchikey': None,
            'target_chemical': None,
            'sink_chemical': None,
            'thermo_dg_m_formation': None,
            'cofactor': None,
        }
        # Refine if needed
        try:
            node['thermo_dg_m_gibbs'] = node['thermo_dg_m_gibbs']['value']
        except TypeError:
            pass
        # Collect
        if node['id'] not in nodes:
            nodes[node['id']] = node
        else:
            try:
                assert _nodes_seem_equal(node, nodes[node['id']])
            except AssertionError:
                logging.error(
                    f'Unexpected node inequality '
                    f'between 2 nodes having ID {node["id"]}.'
                )

    # Node info: compounds
    for cmpd in rp_pathway.get_compounds():
        node = {
            # 'id': _get_specie_node_id_old(specie_dict, specie_id),
            'id': cmpd.get_id(),
            'path_ids': [pathway['path_id'], ],
            'type': 'chemical',
            # 'label': _get_specie_node_id_old(specie_dict, specie_id),
            'label': cmpd.get_id(),
            # 'all_labels': [_get_specie_node_id_old(specie_dict, specie_id), ],
            'all_labels': [cmpd.get_id(), ],
            'svg': None,  # Will be filled later
            'xlinks': _get_specie_xlinks(cmpd),  # TODO: fix me
            # Only for reaction, None for compounds
            'rsmiles': None,
            'rule_ids': None,
            'rxn_template_ids': None,
            'ec_numbers': None,
            'thermo_dg_m_gibbs': None,
            'rule_score': None,
            # Only for compounds
            'smiles': cmpd.get_smiles(),
            'inchi': cmpd.get_inchi(),
            'inchikey': cmpd.get_inchikey(),
            'target_chemical': rp_pathway.get_target_id() == cmpd.get_id(),
            'sink_chemical': cmpd.get_id() in rp_pathway.get_sink(),
            'thermo_dg_m_formation': None,  # FIXME
            'cofactor': False,  # Default, refined later
        }
        # Collect
        if node['id'] not in nodes:
            nodes[node['id']] = node
        else:
            try:
                assert _nodes_seem_equal(node, nodes[node['id']])
                # TODO merge nodes even if they looks equals
            except AssertionError:
                logging.error(
                    f'Unexpected node inequality '
                    f'between 2 nodes having ID {node["id"]}.'
                )

    # Edges
    # for rxn_dict in rpsbml_dict['reactions'].values():
    for rxn in rp_pathway.get_reactions().values():
        rxn_node_id = _get_reaction_node_id(rxn)
        # Reactants
        for cmpd in rxn.get_reactants_compounds():
            cmpd_node_id = cmpd.get_id()
            edge_id = f'{cmpd_node_id}_{rxn_node_id}'
            edge = {
                'id': edge_id,
                'path_ids': [pathway['path_id'],],
                'source': cmpd_node_id,
                'target': rxn_node_id
            }
            if edge_id not in edges:
                edges[edge_id] = edge
            else:
                try:
                    assert _edge_seem_equal(edge, edges[edge_id])
                except AssertionError:
                    logging.error(
                        f'Unexpected edge inequality '
                        f'between 2 edges having ID {edge_id}.'
                    )
        # Products
        for cmpd in rxn.get_products_compounds():
            cmpd_node_id = cmpd.get_id()
            edge_id = f'{rxn_node_id}_{cmpd_node_id}'
            edge = {
                'id': edge_id,
                'path_ids': [pathway['path_id'],],
                'source': rxn_node_id,
                'target': cmpd_node_id
            }
            if edge_id not in edges:
                edges[edge_id] = edge
            else:
                try:
                    assert _edge_seem_equal(edge, edges[edge_id])
                except AssertionError:
                    logging.error(
                        f'Unexpected edge inequality '
                        f'between 2 edges having ID {edge_id}.'
                    )
    # Update pathway info
    pathway['node_ids'] = list(nodes.keys())
    pathway['edge_ids'] = list(edges.keys())

    return nodes, edges, pathway


def parse_all_pathways(input_files: list) -> tuple:
    """Parse all pathways from a list of SBML files.

    Parameters
    ----------
    input_files : list
        List of SBML file paths to parse.

    Returns
    -------
    tuple
        A tuple containing:
        - network: dict, a dictionary representing the network of elements.
        - pathways_info: dict, a dictionary containing information about each
          pathway.
    """
    network = {'elements': {'nodes': [], 'edges': []}}
    all_nodes = {}
    all_edges = {}
    pathways_info = {}

    for sbml_path in input_files:
        rpsbml = rpSBML(str(sbml_path))
        pathway = rpPathway.from_rpSBML(rpsbml=rpsbml)
        nodes, edges, pathway = parse_one_pathway(pathway)
        # Store pathway
        pathways_info[pathway['path_id']] = pathway
        # Store nodes
        for node_id, node_dict in nodes.items():
            if node_id in all_nodes:
                all_nodes[node_id] = _merge_nodes(node_dict, all_nodes[node_id])
            else:
                all_nodes[node_id] = node_dict
        # Store edges
        for edge_id, edge_dict in edges.items():
            if edge_id in all_edges:
                all_edges[edge_id] = _merge_edges(edge_dict, all_edges[edge_id])
            else:
                all_edges[edge_id] = edge_dict

    # Finally store nodes
    for node in all_nodes.values():
        network['elements']['nodes'].append({'data': node})
    for edge in all_edges.values():
        network['elements']['edges'].append({'data': edge})

    # Finally, sort node and edge IDs everywhere
    for node in network['elements']['nodes']:
        node['data']['path_ids'] = sorted(node['data']['path_ids'])
    for node in network['elements']['edges']:
        node['data']['path_ids'] = sorted(node['data']['path_ids'])
    # Finally, sort pathway_info by pathway ID
    pathways_info_ordered = {}
    path_ids_ordered = sorted(pathways_info.keys())
    for path_id in path_ids_ordered:
        pathways_info_ordered[path_id] = pathways_info[path_id]

    return network, pathways_info_ordered


def annotate_cofactors(network: Dict, cofactor_file: str) -> Dict:
    """Annotate cofactors based on structures listed in the cofactor file.

    Parameters
    ----------
    network : dict
        Network of elements as outputted by the sbml_to_json method.
    cofactor_file : str
        File path to the cofactor file.

    Returns
    -------
    dict
        Network annotated with cofactor information.
    """
    if not os.path.exists(cofactor_file):
        logging.error('Cofactor file not found: %s', cofactor_file)
        return network

    # Collect cofactor IDs and structures
    cof_inchis = set()
    cof_ids = set()
    with open(cofactor_file, 'r', encoding='utf-8') as ifh:
        reader = csv.DictReader(ifh, delimiter='\t')
        for row in reader:
            if row['ID'].startswith('#'):
                # Skip row starting with comments
                continue
            if row['INCHI'] != '':
                # InChI describing cofactors
                cof_inchis.add(row['INCHI'])
            if row['ID'] != '':
                # IDs of cofactors
                cof_ids |= set(row['ID'].split(','))

    # Match and annotate network elements
    for node in network['elements']['nodes']:
        if (
            node['data']['type'] == 'chemical'
            and node['data']['inchi'] is not None
        ):
            match = False
            for inchi in cof_inchis:
                if node['data']['inchi'].find(inchi) > -1:  # Match
                    node['data']['cofactor'] = True
                    match = True
                    continue
            if match:
                continue  # Optimisation
        if (
            node['data']['type'] == 'chemical'
            and len(node['data']['all_labels']) > 0
        ):
            match = False
            for id_ in cof_ids:
                for label in node['data']['all_labels']:
                    if label == id_:
                        node['data']['cofactor'] = True
                        match = True
                        continue
            if match:
                continue

    return network


def annotate_chemical_svg(network: Dict) -> Dict:
    """Annotate chemical nodes with SVGs depiction.

    Parameters
    ----------
    network : dict
        Network of elements as outputted by the sbml_to_json method.

    Returns
    -------
    dict
        Network annotated with SVG depictions of chemical nodes.
    """
    from rdkit.Chem import MolFromInchi
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.AllChem import Compute2DCoords
    from urllib import parse

    for node in network['elements']['nodes']:
        if node['data']['type'] == 'chemical' \
                and node['data']['inchi'] is not None \
                and node['data']['inchi'] != '':
            inchi = node['data']['inchi']
            try:
                mol = MolFromInchi(inchi)
                # if mol is None:
                #     raise BaseException('Mol is None')
                Compute2DCoords(mol)
                drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg_draft = drawer.GetDrawingText().replace("svg:", "")
                svg = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg_draft)
                node['data']['svg'] = svg
            except BaseException as e:
                msg = 'SVG depiction failed from inchi: "{}"'.format(inchi)
                logging.warning(msg)
                logging.warning("Below the RDKit backtrace...")
                logging.warning(e)
                node['data']['svg'] = None

    return network


def get_autonomous_html(ifolder):
    """Merge all needed file into a single HTML

    :param ifolder: folder containing the files to be merged
    :return html_str: string, the HTML
    """
    # find and open the index file
    html_string = open(ifolder + '/index.html', 'rb').read()
    # open and read JS files and replace them in the HTML
    js_replace = [
        'js/chroma-2.1.0.min.js',
        'js/cytoscape-3.19.0.min.js',
        'js/cytoscape-dagre-2.3.2.js',
        'js/dagre-0.8.5.min.js',
        'js/jquery-3.6.0.min.js',
        'js/jquery-ui-1.12.1.min.js',
        'js/jquery.tablesorter-2.31.3.min.js',
        'js/viewer.js'
    ]
    for js in js_replace:
        js_string = open(ifolder + '/' + js, 'rb').read()
        ori = b'src="' + js.encode() + b'">'
        rep = b'>' + js_string
        html_string = html_string.replace(ori, rep)
    # open and read style.css and replace it in the HTML
    css_replace = [
        'css/jquery.tablesorte.theme.default-2.31.2.min.css',
        'css/viewer.css'
    ]
    for css_file in css_replace:
        css_bytes = open(ifolder + '/' + css_file, 'rb').read()
        ori = b'<link href="' + css_file.encode() + b'" rel="stylesheet" type="text/css"/>'
        rep = b'<style type="text/css">' + css_bytes + b'</style>'
        html_string = html_string.replace(ori, rep)
    # replace the network
    net_string = open(ifolder + '/network.json', 'rb').read()
    ori = b'src="' + 'network.json'.encode() + b'">'
    rep = b'>' + net_string
    html_string = html_string.replace(ori, rep)
    return html_string


if __name__ == '__main__':
    pass
