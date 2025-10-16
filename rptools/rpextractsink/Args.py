from argparse import ArgumentParser

default_comp = 'c'
default_cspace = 'mnx3.1'


def add_arguments(parser: ArgumentParser) -> ArgumentParser:
    parser.add_argument(
        'input_sbml',
        type=str,
        help="input SBML file"
    )
    parser.add_argument(
        'output_sink',
        type=str,
        help="output sink file"
    )
    parser.add_argument(
        '--compartment-id',
        type=str,
        default=default_comp,
        help=(f"SBML compartment id from which to extract "
              f"the chemical species (default: {default_comp})")
    )
    parser.add_argument(
        '--remove-dead-end',
        action='store_true',
        help=('upon FVA evaluation, ignore chemical'
             'species that do not have any flux')
    )
    parser.add_argument(
        '--chemical-space',
        dest='cspace',
        type=str,
        default=default_cspace,
        help=(f'Type of cache data to use (default: {default_cspace})')
    )
    parser.add_argument(
        '--standalone',
        action='store_true',
        help=('do not get the InChI from Internet')
    )
    return parser
