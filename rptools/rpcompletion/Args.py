from argparse import ArgumentParser
from os import path as os_path


default_upper_flux_bound = 10000
default_lower_flux_bound = -default_upper_flux_bound
default_max_subpaths_filter = 10


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument(
        'rp2_metnet',
        type=str,
        help='Retrosynthesis network provided by RetroPath2.0'
    )
    # parser.add_argument(
    #     '--rp2_metnet',
    #     nargs='?',
    #     type=str,
    #     help='Retrosynthesis network provided by RetroPath2.0'
    # )
    parser.add_argument(
        'sink',
        type=str,
        help='List of compounds in the sink'
    )
    # parser.add_argument(
    #     '--rp2_sink',
    #     nargs='?',
    #     type=str,
    #     help='List of compounds in the sink'
    # )
    parser.add_argument('rp2paths_compounds', type=str)
    # parser.add_argument('--rp2paths_compounds', nargs='?', type=str)
    parser.add_argument('rp2paths_pathways', type=str)
    # parser.add_argument('--rp2paths_pathways', nargs='?', type=str)
    parser.add_argument('outdir', type=str)
    parser.add_argument(
        '--cache-dir',
        default='',
        type=str,
        help='Path to the cache to generate or read from'
    )
    # parser.add_argument('--outdir', nargs='?', type=str)
    # parser.add_argument(
    #     '--out_format',
    #     type=str,
    #     default='JSON',
    #     choices=['RPSBML', 'JSON', 'rpsbml', 'json']
    # )
    parser.add_argument('--upper_flux_bound', type=int, default=default_upper_flux_bound)
    parser.add_argument('--lower_flux_bound', type=int, default=default_lower_flux_bound)
    parser.add_argument(
        '--max_subpaths_filter',
        type=int,
        default=default_max_subpaths_filter,
        help=f'Define the topX pathways to keep (default: {default_max_subpaths_filter}, 0 = no filtering)')
    # parser.add_argument('--pathway_id', type=str, default='rp_pathway')
    # parser.add_argument('--compartment_id', type=str, default='MNXC3')
    # parser.add_argument('--species_group_id', type=str, default='rp_trunk_species')
    # parser.add_argument('--sink_species_group_id', type=str, default='rp_sink_species')
    # parser.add_argument('--pubchem_search', type=str, default='False')
    return parser
