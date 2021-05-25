from argparse import ArgumentParser


def add_arguments(parser):

    # positional arguments
    parser.add_argument('infile',    type=str, help='pathway as JSON file')
    parser.add_argument('outfile',   type=str, help='updated pathway as JSON file')

    # optional arguments
    parser.add_argument('--pathway_id'     , type=str   , default='rp_pathway')
    parser.add_argument('--ph'             , type=float)
    parser.add_argument('--ionic_strength' , type=float)
    parser.add_argument('--pMg'            , type=float)
    parser.add_argument('--temp_k'         , type=float)

    return parser
