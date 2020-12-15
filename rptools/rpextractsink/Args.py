from argparse import ArgumentParser


def build_args_parser():
    parser = ArgumentParser(prog='rpextractsink', description='Generate the sink from a model SBML by specifying the compartment')
    parser = _add_arguments(parser)

    return parser


def _add_arguments(parser):
    parser.add_argument('input_sbml',
                        type=str,
                        help="input SBML file")
    parser.add_argument('output_sink',
                        type=str,
                        help="output sink file")
    parser.add_argument('--compartment_id',
                        type=str,
                        default='MNXC3',
                        help='SBML compartment id from which to extract the chemical species')
    parser.add_argument('--remove_dead_end',
                        action='store_true',
                        help='upon FVA evaluation, ignore chemical species that do not have any flux')
    parser.add_argument('--log', metavar='ARG',
                        type=str, choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='error',
                        help='Adds a console logger for the specified level (default: error)')
    return parser
