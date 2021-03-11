from argparse  import ArgumentParser


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument('rp2_pathways', type=str)
    parser.add_argument('rp2paths_compounds', type=str)
    parser.add_argument('rp2paths_pathways', type=str)
    parser.add_argument('outdir', type=str)
    parser.add_argument('--upper_flux_bound', type=int, default=999999)
    parser.add_argument('--lower_flux_bound', type=int, default=0)
    parser.add_argument('--max_subpaths_filter', type=int, default=0)
    parser.add_argument('--pathway_id', type=str, default='rp_pathway')
    parser.add_argument('--compartment_id', type=str, default='MNXC3')
    parser.add_argument('--species_group_id', type=str, default='central_species')
    parser.add_argument('--sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('--pubchem_search', type=str, default='False')
    return parser
