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

def add_arguments(parser):
    # parser.add_argument(
    #     'infile',
    #     type = str,
    #     help = 'Pathway file (rpSBML) with scores (rules, FBA, Thermo...)'
    # )
    # parser.add_argument(
    #     'outfile',
    #     type = str,
    #     help = 'Path to write pathway file (rpSBML) with global score'
    # )
    parser.add_argument('-ttdf', '--test_data_file', required=True, type=str)
    parser.add_argument('-ttsf', '--test_score_file', required=True, type=str)
    parser.add_argument(
        '--data_train_file',
        type=str,
        default=os_path.join(
            __MODELS_PATH,
            'data_train.h5'
        )
    )
    parser.add_argument(
        '--data_predict_file',
        type=str,
        default=os_path.join(
            __MODELS_PATH,
            'data_predict.h5'
        )
    )

    return parser
