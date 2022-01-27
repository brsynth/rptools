from argparse import ArgumentParser

DEFAULT_pH = 7.5
DEFAULT_pMg = 3.0
DEFAULT_ionic_strength = 0.25

def add_arguments(parser: ArgumentParser) -> ArgumentParser:

    # positional arguments
    parser.add_argument('infile',    type=str, help='pathway as rpSBML file')
    parser.add_argument('outfile',   type=str, help='updated pathway as rpSBML file')

    # optional arguments
    parser.add_argument(
        '--pH',
        type=float,
        default=DEFAULT_pH
    )
    parser.add_argument(
        '--ionic_strength',
        type=float,
        default=DEFAULT_ionic_strength
    )
    parser.add_argument(
        '--pMg',
        type=float,
        default=DEFAULT_pMg
    )
    # parser.add_argument('--temp_k'         , type=float)

    return parser
