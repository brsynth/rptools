from unittest import TestCase
from numpy import array as np_array
from rptools.rpthermo.rpThermo import (
    minimize,
    build_stoichio_matrix,
    remove_unknown_compounds,
    get_target_rxn_idx,
    minimize
)
from brs_utils import (
    create_logger,
)
from chemlite import Reaction


class Test_rpThermo(TestCase):

    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')
        self.rxn_1 = Reaction(
            id='rxn_1',
            reactants={'MNXM188': 1, 'MNXM4': 1, 'MNXM6': 1, 'MNXM1': 3},
            products={'CMPD_0000000004': 1, 'CMPD_0000000003': 1, 'MNXM13': 1, 'MNXM15': 3, 'MNXM5': 1},
        )
        self.rxn_2 = Reaction(
            id='rxn_2',
            reactants={'MNXM4': 1, 'CMPD_0000000003': 2},
            products={'MNXM1': 1, 'TARGET_0000000001': 1},
        )
        self.rxn_3 = Reaction(
            id='rxn_3',
            reactants={'CMPD_0000000004': 3, 'MNXM4': 1, 'MNXM6': 1},
            products={'MNXM13': 1, 'MNXM5': 1},
        )
        self.reactions = [self.rxn_1, self.rxn_2, self.rxn_3]
        self.sto_mat_1 = [
            [-3.0, 1.0, 0.0],
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 1.0],
            [3.0, 0.0, 0.0],
            [1.0, 0.0, -3.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, -2.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-1.0, 0.0, -1.0]
        ]

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

    def test_build_stoichio_matrix(self):
        # Ignore the order of matrix lines because
        # it is not relevant for our resolution system
        self.assertCountEqual(
            build_stoichio_matrix(self.reactions).tolist(),
            self.sto_mat_1
        )

    def test_build_stoichio_matrix_w_sel_cmpds(self):
        # Ignore the order of matrix lines because
        # it is not relevant for our resolution system
        self.assertCountEqual(
            build_stoichio_matrix(
                reactions=self.reactions,
                compounds=['CMPD_0000000003']
            ).tolist(),
            [self.sto_mat_1[7]]
        )

    def test_get_target_rxn_idx(self):
        self.assertEqual(
            get_target_rxn_idx(
                reactions=self.reactions,
                rxn_target_id=self.rxn_2.get_id(),
            ),
            self.reactions.index(self.rxn_2)
        )

    def test_remove_unknown_compounds(self):
        compd_id = 'UNK_CMPD_FOOBAR'
        self.rxn_1.add_product(stoichio=1, compound_id=compd_id)
        self.rxn_3.add_reactant(stoichio=1, compound_id=compd_id)
        sto_mat = build_stoichio_matrix(self.reactions)
        print(sto_mat)
        reactions = remove_unknown_compounds(
            unk_compounds=[compd_id],
            reactions=self.reactions,
            rxn_target_id=self.rxn_2.get_id(),
        )
        for rxn in reactions:
            print(rxn)
        exit()
