from unittest import TestCase
from numpy import array as np_array
# from os import path as os_path
from rptools.rpthermo.rpThermo import (
    minimize,
)
# from tempfile import (
#     TemporaryDirectory,
#     mkdtemp
# )
# from shutil    import rmtree
from brs_utils import (
    create_logger,
)


class Test_rpThermo(TestCase):

    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

    # Reactions
    # |- rxn_1: 1.0 MNXM188 + 1.0 MNXM4 + 1.0 MNXM6 + 3.0 MNXM1 --> 1.0 CMPD_0000000004 + 1.0 CMPD_0000000003 + 1.0 MNXM13 + 1.0 MNXM15 + 1.0 MNXM5
    # |- rxn_2: 1.0 MNXM4 + 2.0 CMPD_0000000003 --> 2.0 MNXM1 + 1.0 TARGET_0000000001
    # |- rxn_3: 1.0 MNXM4 + 1.0 MNXM6 + 3.0 CMPD_0000000004 --> 1.0 MNXM13 + 1.0 MNXM5
    # Compounds ordered: CMPD_0000000003, CMPD_0000000004
    # Reactions ordered: rxn_1, rxn_2*, rxn_3
    # * to be optimised
    def test_minimize_1(self):
        sto_mat = np_array(
            [
                [1, -2, 0],
                [1, 0, -3]
            ]
        )
        rxn_tgt_idx = 1
        coeffs = minimize(
            sto_mat,
            rxn_tgt_idx,
            self.logger
        )
        self.assertSequenceEqual(
            coeffs.tolist(),
            [1 , 1/2 , 1/3]
        )

    # Reactions
    # |- rxn_1: 1.0 MNXM4 + 1.0 MNXM421 + 1.0 MNXM6 + 1.0 MNXM1 --> 1.0 CMPD_0000000015 + 1.0 MNXM2 + 1.0 MNXM5
    # |- rxn_2: 1.0 MNXM1 + 1.0 CMPD_0000000015 + 1.0 MNXM2 --> 1.0 CMPD_0000000010 + 1.0 MNXM15
    # |- rxn_3: 1.0 MNXM1 + 1.0 CMPD_0000000010 --> 1.0 CMPD_0000000003 + 1.0 MNXM13
    # |- rxn_4: 1.0 MNXM4 + 1.0 CMPD_0000000003 --> 2.0 MNXM1 + 1.0 TARGET_0000000001
    # Compounds ordered: CMPD_0000000003, CMPD_0000000010, CMPD_0000000015
    # Reactions ordered: rxn_1, rxn_2, rxn_3, rxn_4*
    # * to be optimised
    def test_minimize_2(self):
        sto_mat = np_array(
            [
                [0, 0, 1, -1],

                [0, 1, -1, 0],

                [1, -1, 0, 0]
            ]
        )
        rxn_tgt_idx = 3
        coeffs = minimize(
            sto_mat,
            rxn_tgt_idx,
            self.logger
        )
        self.assertSequenceEqual(
            coeffs.tolist(),
            [1 , 1 , 1, 1]
        )

    def test_minimize_1cmpd(self):
        sto_mat = np_array(
            [
                [ 1, -1,  0,  0]
            ]
        )
        rxn_tgt_idx = 2
        coeffs = minimize(
            sto_mat,
            rxn_tgt_idx,
            self.logger
        )
        self.assertSequenceEqual(
            coeffs.tolist(),
            [1, 1, 1, 1]
        )