from unittest import TestCase
from os import path as os_path
from rptools.rplibs import rpPathway
from brs_utils import (
    create_logger,
)

from rptools.rpscore import (
    predict_score,
    ThermoError,
    FBAError
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
            # 0.1922311931848526,
            0.9923127055168152,
            predict_score(
                pathway=self.pathway,
                no_of_rxns_thres=10
            )
        )

    def test_score_woThermo(self):
        pathway_path = os_path.join(
            self.data_path,
            'pathway_woThermo.xml'
        )
        self.pathway = rpPathway.from_rpSBML(
            infile=pathway_path,
            logger=self.logger
        )
        with self.assertRaises(ThermoError):
            predict_score(
                pathway=self.pathway,
                no_of_rxns_thres=10
            )

    def test_score_woFBA(self):
        pathway_path = os_path.join(
            self.data_path,
            'pathway_woFBA.xml'
        )
        self.pathway = rpPathway.from_rpSBML(
            infile=pathway_path,
            logger=self.logger
        )
        with self.assertRaises(FBAError):
            predict_score(
                pathway=self.pathway,
                no_of_rxns_thres=10
            )
