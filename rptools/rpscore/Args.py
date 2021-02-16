from argparse  import ArgumentParser
from rptools._version import __version__


def build_args_parser():
    parser = ArgumentParser(
        prog='rpscore',
        description='Calculate global score by combining all scores (rules, FBA, Thermo)'
    )
    parser = _add_arguments(parser)
    return parser


def _add_arguments(parser):
    parser.add_argument(
        'pathway_file',
        type = str,
        help = 'Input rpSBML filename with scores (rules, FBA, Thermo...)'
    )
    parser.add_argument(
        '--outfile',
        type = str,
        default = '',
        help = 'Path to write rpSBML with global score'
    )
    parser.add_argument(
        '--weight_rp_steps',
        type = float,
        default = '0.10002239003499142',
        help = 'The weight associated with the number of steps (Default: 0.10002239003499142)'
    )
    parser.add_argument(
        '--weight_rule_score',
        type = float,
        default = '0.13346271414277305',
        help = 'The weight associated with the mean of reaction rule scores (Default: 0.13346271414277305)'
    )
    parser.add_argument(
        '--weight_fba',
        type = float,
        default = '0.6348436269211155',
        help = 'The weight associated with the flux of the target (Default: 0.6348436269211155)'
    )
    parser.add_argument(
        '--weight_thermo',
        type = float,
        default = '0.13167126890112002',
        help = 'The weight associated with the sum of reaction Gibbs free energy (Default: 0.13167126890112002)'
    )
    parser.add_argument(
        '--max_rp_steps',
        type = int,
        default = '15',
        help = 'The maximal number of steps are run in RP2 (Default: 15)'
    )
    parser.add_argument(
        '--thermo_ceil',
        type = float,
        default = '5000.0',
        help = 'The upper limit of Gibbs free energy for each reaction (Default: 5000.0)'
    )
    parser.add_argument(
        '--thermo_floor',
        type = float,
        default = '-5000.0',
        help = 'The lower limit of Gibbs free energy for each reaction (Default: -5000.0)'
    )
    parser.add_argument(
        '--fba_ceil',
        type = float,
        default = '5.0',
        help = 'The upper flux limit of the heterologous pathway (Default: 5.0)'
    )
    parser.add_argument(
        '--fba_floor',
        type = float,
        default = '0.0',
        help = 'The lower flux limit of the heterologous pathway (Default: 5.0)'
    )
    parser.add_argument(
        '--pathway_id',
        type = str,
        default = 'rp_pathway',
        help = 'The ID of the heterologous pathway (Default: rp_pathway)'
    )
    parser.add_argument(
        '--objective_id',
        type = str,
        default = 'obj_fraction',
        help = 'The ID of the FBA objective (Default: obj_fraction)'
    )
    parser.add_argument(
        '--thermo_id',
        type = str,
        default = 'dfG_prime_m',
        help = 'The ID of the Gibbs free energy that may be used. May be either dfG_prime_m or dfG_prime_o (Default: dfG_prime_m)'
    )
    parser.add_argument(
        '--log',
        metavar='ARG',
        type=str,
        choices=[
            'debug', 'info', 'warning', 'error', 'critical', 'silent', 'quiet'
            'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'SILENT', 'QUIET'
        ],
        default='def_info',
        help='Adds a console logger for the specified level (default: error)'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(__version__),
        help='show the version number and exit'
    )
    return parser
