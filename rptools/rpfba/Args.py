from argparse  import ArgumentParser
from rptools._version import __version__

default_upper_flux_bound = 10000
default_lower_flux_bound = 0

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
        help='model compartment id to consider (e.g. \'c\' or \'MNXC3\')'
    )
    parser.add_argument(
        'outfile',
        type=str,
        help='output file'
    )
    parser.add_argument(
        '--pathway_id',
        type=str,
        default='rp_pathway',
        help='id of the heterologous pathway (default: rp_pathway)'
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
        help='Biomass reaction ID (default: biomass)'
    )
    parser.add_argument(
        '--sim',
        type=str,
        choices=['fba', 'pfba', 'fraction'],
        default='fraction',
        help='type of simulation to use (default: fraction)'
    )
    # parser.add_argument(
    #     '--source_reaction',
    #     type=str,
    #     default='biomass',
    #     help='reaction id of the source reaction (default: biomass)'
    # )
    # parser.add_argument(
    #     '--target_reaction',
    #     type=str,
    #     default='rxn_target',
    #     help='reaction id of the target reaction (default: rxn_target). Note: if \'fba\' or \'rpfba\' options are used, then these are ignored'
    # )
    parser.add_argument(
        '--biomass_coeff',
        type=float,
        default=1.0,
        help='source coefficient (default: 1.0)'
    )
    parser.add_argument(
        '--objective_coeff',
        type=float,
        default=1.0,
        help='target coefficient (default: 1.0)'
    )
    # parser.add_argument('--num_workers',
    #                     type=int,
    #                     default=10,
    #                     help='number of workers (multi-threads)')
    parser.add_argument(
        '--is_max',
        action='store_true',
        default=True,
        help='maximise the objective'
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
    parser.add_argument(
        '--upper_flux_bound',
        type=float,
        default=default_upper_flux_bound,
        help='flux constaints (upper bound) for FBA (default: 10000)'
    )
    parser.add_argument(
        '--lower_flux_bound',
        type=float,
        default=default_lower_flux_bound,
        help='Flux constaints (lower bound) for FBA (default: 0)'
    )
    return parser
