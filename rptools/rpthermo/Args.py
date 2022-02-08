from argparse import ArgumentParser

DEFAULT_pH = 7.5
DEFAULT_pMg = 3.0
DEFAULT_ionic_strength = 0.25
MIN_pH = 0
MAX_pH = 14
MIN_ionic_strength = 0
MAX_ionic_strength = 500


def add_arguments(parser: ArgumentParser) -> ArgumentParser:

    # positional arguments
    parser.add_argument('infile',    type=str, help='pathway as rpSBML file')
    parser.add_argument('outfile',   type=str, help='updated pathway as rpSBML file')

    # optional arguments
    parser.add_argument(
        '--pH',
        type=float,
        choices=range(MIN_pH, MAX_pH+1),
        metavar=f'[{MIN_pH}-{MAX_pH}]',
        default=DEFAULT_pH
    )
    parser.add_argument(
        '--ionic_strength',
        type=float,
        choices=range(MIN_ionic_strength, MAX_ionic_strength+1),
        metavar=f'[{MIN_ionic_strength}-{MAX_ionic_strength}]',
        default=DEFAULT_ionic_strength
    )
    parser.add_argument(
        '--pMg',
        type=float,
        default=DEFAULT_pMg
    )
    # parser.add_argument('--temp_k'         , type=float)

    return parser
