from os import path as os_path
from argparse  import ArgumentParser
from rptools._version import __version__

__CURRENT_PATH = os_path.dirname(
    os_path.abspath(__file__)
)
__MODELS_PATH = os_path.join(
    __CURRENT_PATH,
    'models'
)
__DATA_TRAIN_FILE = os_path.join(
    __MODELS_PATH,
    'data_train.h5'
)

def add_arguments(parser):
    parser.add_argument(
        'pathways',
        type = str,
        nargs="+",
        help='Pathway file(s) (rpSBML) with scores (rules, FBA, Thermo...)'
    )
    parser.add_argument(
        '--outdir',
        type=str,
        help='Path to write pathway files (rpSBML) with global score (default: out)',
        default='out'
    )
    parser.add_argument(
        '--outfile',
        type=str,
        help='Path to write pathway file (rpSBML) with global score',
    )
    parser.add_argument(
        '--no_of_rxns_thres',
        type=int,
        help='Number of reactions above which pathway are not scored (too long) (default: 10)',
        default=10
    )
    # parser.add_argument('-ttdf', '--test_data_file', required=True, type=str)
    # parser.add_argument('-ttsf', '--test_score_file', required=True, type=str)
    parser.add_argument(
        '--data_train_file',
        type=str,
        default=__DATA_TRAIN_FILE
    )

    return parser
