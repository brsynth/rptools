from argparse import ArgumentParser
from os import path as os_path


default_upper_flux_bound = 10000
default_lower_flux_bound = -default_upper_flux_bound
default_maxsubpaths = 10
default_cofactors = None
default_cspace = 'mnx3.1'


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
        '--chemical-space',
        dest='cspace',
        default=default_cspace,
        type=str,
        help='Chemical space to use (e.g. mnx3.1, mnx4.0...). Determines which configuration files and folders to use both the cache and the input cache (default: %(default)s).'
    )
    parser.add_argument(
        '--cache-dir',
        dest='cache_dir',
        default=None,
        type=str,
        help='Directory where the rrCache is stored. If not provided, the default directory will be used.'
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
        '--maxsubpaths',
        type=int,
        default=default_maxsubpaths,
        help=f'Define the topX subpathways per master pathway to keep after completion (default: {default_maxsubpaths}, 0 = no filtering)')
    parser.add_argument(
        '--cofactors',
        type=str,
        help='Name of the file containing the list of cofactors to ignore (default: None)',
        default=None
    )
    parser.add_argument(
        '--forward', dest='forward',
        help='Consider reactions in the forward direction',
        required=False, action='store_true'
    )
    # parser.add_argument('--pathway_id', type=str, default='rp_pathway')
    # parser.add_argument('--compartment_id', type=str, default='MNXC3')
    # parser.add_argument('--species_group_id', type=str, default='rp_trunk_species')
    # parser.add_argument('--sink_species_group_id', type=str, default='rp_sink_species')
    # parser.add_argument('--pubchem_search', type=str, default='False')
    return parser
