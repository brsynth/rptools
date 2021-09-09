from unittest import TestCase
from copy import deepcopy
from numpy import array as np_array
from rptools.rpthermo.rpThermo import (
    build_stoichio_matrix,
    get_target_rxn_idx,
    minimize,
    remove_compounds,
    # eQuilibrator,
    # initThermo,
    # get_compounds_from_cache
)
from brs_utils import (
    create_logger,
    Cache
)
from chemlite import (
    Compound,
    Reaction,
    Pathway
)


species = {
    "TARGET_0000000001": Compound(
        id="TARGET_0000000001",
        smiles="[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]",
        inchi="InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)",
        inchikey="TXXHDPDFNKHHGW-UHFFFAOYSA-N"
    ),
    "CMPD_0000000010": Compound(
        id="CMPD_0000000010",
        smiles="[H]OC(=O)c1c([H])c([H])c(O[H])c(O[H])c1[H]",
        inchi="InChI=1S/C7H6O4/c8-5-2-1-4(7(10)11)3-6(5)9/h1-3,8-9H,(H,10,11)",
        inchikey="YQUVCSBJEUQKSH-UHFFFAOYSA-N"
    ),
    "MNXM23": Compound(
        id="MNXM23",
        formula="C3H3O3",
        smiles="CC(=O)C(=O)O]",
        inchi="InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)",
        inchikey="LCTONWCANYUPML-UHFFFAOYSA-N",
        name="pyruvate"
    ),
    "CMPD_0000000025": Compound(
        id="CMPD_0000000025",
        smiles="[H]OC(=O)c1c([H])c([H])c([H])c(O[H])c1[H]",
        inchi="InChI=1S/C7H6O3/c8-6-3-1-2-5(4-6)7(9)10/h1-4,8H,(H,9,10)",
        inchikey="IJFXRHURBJZNAO-UHFFFAOYSA-N"
    ),
    "CMPD_0000000003": Compound(
        id="CMPD_0000000003",
        smiles="[H]Oc1c([H])c([H])c([H])c([H])c1O[H]",
        inchi="InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H",
        inchikey="YCIMNLLNPGFGHC-UHFFFAOYSA-N"
    ),
    "CMPD_0000000003_wo_smiles": Compound(
        id="CMPD_0000000003_wo_smiles",
        inchi="InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H",
        inchikey="YCIMNLLNPGFGHC-UHFFFAOYSA-N"
    ),
    "CMPD_0000000004_wo_smiles": Compound(
        id="CMPD_0000000003_wo_smiles",
        inchi="InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H",
        inchikey="YCIMNLLNPGFGHC-UHFFFAOYSA-N"
    ),
    "CMPD_0000000003_w_smiles_None": Compound(
        id="CMPD_0000000003_wo_smiles",
        inchi="InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H",
        inchikey="YCIMNLLNPGFGHC-UHFFFAOYSA-N",
        smiles=None
    ),
    "MNXM337": Compound(
        id="MNXM337",
        smiles="[H]OC(=O)C(OC1([H])C([H])=C(C(=O)O[H])C([H])=C([H])C1([H])O[H])=C([H])[H]",
        inchi="InChI=1S/C10H10O6/c1-5(9(12)13)16-8-4-6(10(14)15)2-3-7(8)11/h2-4,7-8,11H,1H2,(H,12,13)(H,14,15)",
        inchikey="WTFXTQVDAKGDEY-UHFFFAOYSA-N"
    ),
    "MNXM2": Compound(
        id="MNXM2",
        smiles="[H]O[H]",
        inchi="InChI=1S/H2O/h1H2",
        inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    ),
    "MNXM13": Compound(
        id="MNXM13",
        smiles="O=C=O",
        inchi="InChI=1S/CO2/c2-1-3",
        inchikey="CURLTUGMZLYLDI-UHFFFAOYSA-N",
        formula="CO2",
        name="CO2"
    ),
    "MNXM5": Compound(
        id="MNXM5",
        smiles="N=C(O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
        inchi="InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)/p+1",
        inchikey="XJLXINKUBYWONI-UHFFFAOYSA-O",
        formula="C21H25N7O17P3",
        name="NADP(+)"
    ),
    "MNXM4": Compound(
        id="MNXM4",
        smiles="O=O",
        inchi="InChI=1S/O2/c1-2",
        inchikey="MYMOFIZGZYHOMD-UHFFFAOYSA-N"
    ),
    "MNXM1": Compound(
        id="MNXM1",
        smiles="[H+]",
        inchi="InChI=1S/p+1",
        inchikey="GPRLSGONYQIRFK-UHFFFAOYSA-N"
    ),
    "MNXM6": Compound(
        id="MNXM6",
        smiles="[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]",
        inchi="InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)",
        inchikey="ACFIXJIJDZMPPO-UHFFFAOYSA-N"
    )
}

class Test_rpThermo(TestCase):

    __test__ = False

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
            [3.0, 1.5, 1.0]
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
        _coeffs = deepcopy(coeffs)
        for coeff_idx in range(len(_coeffs)):
            if (
                _coeffs[coeff_idx] == 0
                or _coeffs[coeff_idx] == abs(float("inf"))
            ):
                _coeffs[coeff_idx] = 1.
        self.assertSequenceEqual(
            list(_coeffs),
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

    def test_remove_compounds(self):
        pathway = Pathway(id='thermo')
        for rxn in self.reactions:
            pathway.add_reaction(rxn)
        compd_id1 = 'UNK_CMPD_FOOBAR'
        compd_id2 = 'UNK_CMPD_FOOBAR_2'
        self.rxn_1.add_product(stoichio=2, compound_id=compd_id1)
        self.rxn_1.add_product(stoichio=3, compound_id=compd_id2)
        self.rxn_2.add_reactant(stoichio=2, compound_id=compd_id2)
        self.rxn_3.add_reactant(stoichio=1, compound_id=compd_id1)
        reactions = remove_compounds(
            compounds=[compd_id1, compd_id2],
            reactions=pathway.get_list_of_reactions(),
            rxn_target_id=self.rxn_2.get_id(),
        )
        self.assertDictEqual(
            Reaction.sum_stoichio(reactions),
            {'MNXM1': -1.5, 'MNXM188': -1.0, 'MNXM4': -4.5, 'MNXM6': -3.0, 'CMPD_0000000003': -2.0, 'CMPD_0000000004': -5.0, 'MNXM13': 3.0, 'MNXM15': 3.0, 'MNXM5': 3.0, 'TARGET_0000000001': 1.5}
        )
        # cc = initThermo()
        # species, unk_compounds = get_compounds_from_cache(
        #     compounds=pathway.get_species(),
        #     cc=cc
        # )

        # results = eQuilibrator(
        #     species_stoichio=pathway.net_reaction(),
        #     species=species,
        #     cc=cc
        # )
