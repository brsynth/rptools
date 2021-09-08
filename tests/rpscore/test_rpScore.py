from unittest import TestCase
from os import path as os_path
from rptools.rplibs import rpPathway
from brs_utils import (
    create_logger,
)

from rptools.rpscore import predict_score
from rptools.rpscore.Args import (
    __DATA_TRAIN_FILE as DATA_TRAIN_FILE,
    __MODELS_PATH as MODELS_PATH
)



class Test_rpScore(TestCase):

    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    rpsbml_path = os_path.join(
        data_path,
        'pathway.xml'
    )


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        self.pathway = rpPathway.from_rpSBML(
            infile=self.rpsbml_path,
            logger=self.logger
        )

    def test_score(self):
        self.assertEqual(
            0.1922311931848526,
            predict_score(
                pathways=[self.pathway],
                data_train_file=DATA_TRAIN_FILE,
                models_path=MODELS_PATH,
                no_of_rxns_thres=10
            )[0]
        )