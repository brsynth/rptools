from argparse  import ArgumentParser

def build_args_parser():
    parser = ArgumentParser(prog='rpcompletion', description='Python wrapper to parse RP2 to generate rpSBML collection of unique and complete (cofactors) pathways')
    parser = _add_arguments(parser)

    return parser

def _add_arguments(parser):
    parser.add_argument('rp2_pathways', type=str)
    parser.add_argument('rp2paths_compounds', type=str)
    parser.add_argument('rp2paths_pathways', type=str)
    parser.add_argument('outdir', type=str)
    parser.add_argument('--upper_flux_bound', type=int, default=999999)
    parser.add_argument('--lower_flux_bound', type=int, default=0)
    parser.add_argument('--max_subpaths_filter', type=int, default=10)
    parser.add_argument('--pathway_id', type=str, default='rp_pathway')
    parser.add_argument('--compartment_id', type=str, default='MNXC3')
    parser.add_argument('--species_group_id', type=str, default='central_species')
    parser.add_argument('--sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('--pubchem_search', type=str, default='False')
    parser.add_argument('--log', metavar='ARG',
                        type=str, choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='error',
                        help='Adds a console logger for the specified level (default: error)')
    return parser
