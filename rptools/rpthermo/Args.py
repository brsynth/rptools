from argparse import ArgumentParser


def add_arguments(parser: ArgumentParser) -> ArgumentParser:

    # positional arguments
    parser.add_argument('infile',    type=str, help='pathway as rpSBML file')
    parser.add_argument('outfile',   type=str, help='updated pathway as rpSBML file')

    # optional arguments
    parser.add_argument('--ph'             , type=float)
    parser.add_argument('--ionic_strength' , type=float)
    parser.add_argument('--pMg'            , type=float)
    parser.add_argument('--temp_k'         , type=float)

    return parser
