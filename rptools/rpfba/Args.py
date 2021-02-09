from argparse  import ArgumentParser
from rptools._version import __version__


def build_args_parser():
    parser = ArgumentParser(
        prog='rpfba',
        description='Calculate FBA to generate rpFBA collection'
    )
    parser = _add_arguments(parser)
    return parser


def _add_arguments(parser):
    parser.add_argument(
        'pathway',
        type=str,
        help='SBML file that contains an heterologous pathway'
    )
    parser.add_argument(
        'model',
        type=str,
        help='GEM model file')
    parser.add_argument(
        'outfile',
        type=str,
        help='output file')
    parser.add_argument(
        '--pathway_id',
        type=str,
        default='rp_pathway',
        help='id of the heterologous pathway (default: rp_pathway)'
    )
    parser.add_argument(
        '--sink_species_group_id',
        type=str,
        default='rp_sink_species',
        help='id of the central species (default: central_species)'
    )
    parser.add_argument(
        '--species_group_id',
        type=str,
        default='central_species',
        help='id of the sink species (default: rp_sink_species)'
    )
    parser.add_argument(
        '--objective_id',
        type=str,
        default=None,
        help='overwrite the auto-generated id of the results (default: None)'
    )
    parser.add_argument(
        '--compartment_id',
        type=str,
        default='MNXC3',
        help='SBML compartment id (default: MNXC3)'
    )
    parser.add_argument(
        '--sim',
        type=str,
        choices=['fba', 'pfba', 'fraction'],
        default='fraction',
        help='type of simulation to use (default: fraction)'
    )
    parser.add_argument(
        '--source_reaction',
        type=str,
        default='biomass',
        help='reaction id of the source reaction (default: biomass)'
    )
    parser.add_argument(
        '--target_reaction',
        type=str,
        default='rxn_target',
        help='reaction id of the target reaction (default: rxn_target). Note: if \'fba\' or \'rpfba\' options are used, then these are ignored'
    )
    parser.add_argument(
        '--source_coefficient',
        type=float,
        default=1.0,
        help='source coefficient (default: 1.0)'
    )
    parser.add_argument(
        '--target_coefficient',
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
        '--dont_merge',
        action='store_true',
        help='output the merged model (default)'
    )
    parser.add_argument(
        '--log',
        metavar='ARG',
        type=str,
        choices=[
            'debug', 'info', 'warning', 'error', 'critical',
            'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'
        ],
        default='def_info',
        help='Adds a console logger for the specified level (default: error)'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )
    return parser
