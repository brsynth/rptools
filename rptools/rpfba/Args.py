from argparse import ArgumentParser
from typing import (
    List,
)
from rptools._version import __version__
from rptools.rpfba.medium import (
    __MEDIUM_DEFAULT_ID,
    __MEDIUM_PATH,
    read_medium_ids
)

def add_arguments(
    parser: ArgumentParser):
    parser.add_argument(
        'pathway_file',
        type=str,
        help='SBML file that contains an heterologous pathway'
    )
    parser.add_argument(
        'model_file',
        type=str,
        help='GEM model file (SBML)'
    )
    parser.add_argument(
        'compartment_id',
        type=str,
        help='model compartment id to consider (e.g. \'c\' or \'MNXC3\')'
    )
    parser.add_argument(
        'outfile',
        type=str,
        help='output file'
    )
    parser.add_argument(
        '--objective_rxn_id',
        type=str,
        default='rxn_target',
        help='reaction ID to optimise (default: rxn_target)'
    )
    parser.add_argument(
        '--biomass_rxn_id',
        type=str,
        default='biomass',
        help='biomass reaction ID (default: biomass). Note: Only for \'fraction\' simulation'
    )
    parser.add_argument(
        '--sim',
        type=str,
        choices=['fba', 'pfba', 'fraction'],
        default='fraction',
        help='type of simulation to use (default: fraction)'
    )
    parser.add_argument(
        '--fraction_of',
        type=float,
        default=0.75,
        help='fraction of the optimum (default: 0.75). Note: this value is ignored is \'fba\' is used'
    )
    parser.add_argument(
        '--merge',
        action='store_true',
        default=False,
        help='output the full merged model in addition of heterologous pathway only (default: False)'
    )
    parser.add_argument(
        '--ignore_orphan_species',
        action='store_true',
        default=True,
        help='ignore metabolites that are only consumed or produced (default: True)'
    )

    parser_medium = parser.add_argument_group('Medium', 'Medium modifications')
    parser_medium.add_argument('--medium_compartment_id',
        type=str,
        default='MNXC2',
        help='Model compartiment id corresponding to the extra-cellular compartment'
    )
    parser_medium.add_argument(
        '--medium_file',
        type=str,
        help='Provide a csv file with an header as <coumpond_id>,<upper_bound>. \
                This file provides information about metabolites (Metanetx Id) to add or remove (giving upper bound value)'
    )
    parser_medium.add_argument(
        '--medium_id',
        type=str,
        default=__MEDIUM_DEFAULT_ID,
        choices=[__MEDIUM_DEFAULT_ID] + read_medium_ids(__MEDIUM_PATH),
        help='Use a base medium composition. Data can be add with the option --medium_file'
    )
    parser_medium.add_argument(
        '--minimal_medium_file',
        type=str,
        help='Provide a path of a CSV file to output information about minimal medium'
    )
 
    return parser
