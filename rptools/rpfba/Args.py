from argparse  import ArgumentParser
from rptools._version import __version__


def add_arguments(parser):
    parser.add_argument(
        'pathway',
        type=str,
        help='SBML file that contains an heterologous pathway'
    )
    parser.add_argument(
        'model',
        type=str,
        help='GEM model file (SBML)'
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
        default='fraction',
        help='overwrite the auto-generated id of the results (default: fraction)'
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
        '--merge',
        action='store_true',
        default=False,
        help='output the full merged model instead of heterologous pathway only (default: False)'
    )
    parser.add_argument(
        '--dont_ignore_orphan_species',
        action='store_true',
        default=False,
        help='Ignore metabolites that are only consumed or produced (default: False)'
    )
    return parser
