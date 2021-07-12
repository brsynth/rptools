from argparse  import ArgumentParser
from rptools._version import __version__


def add_arguments(parser):
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
        help='model compartment id to consider (e.g. \'c\' or \'MNXC3\' or \'c|MNXC3|cytosol|cytoplasm\')'
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
        help='output the full merged model instead of heterologous pathway only (default: False)'
    )
    parser.add_argument(
        '--ignore_orphan_species',
        action='store_true',
        default=True,
        help='ignore metabolites that are only consumed or produced (default: True)'
    )
    return parser
