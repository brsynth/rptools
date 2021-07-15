"""
Created on May 28 2021

@author: Joan Hérisson
"""

from unittest import TestCase
from copy import deepcopy
from rr_cache import rrCache
from rptools.rplibs import (
    rpPathway,
    rpReaction,
    rpCompound
)


class Test_rpPathway(TestCase):

    species = {
        "TARGET_0000000001": rpCompound(
            id="TARGET_0000000001",
            smiles="[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]",
            inchi="InChI=1S/C6H6O4/c7-5(8)3-1-2-4-6(9)10/h1-4H,(H,7,8)(H,9,10)",
            inchikey="TXXHDPDFNKHHGW-UHFFFAOYSA-N"
        ),
        "CMPD_0000000010": rpCompound(
            id="CMPD_0000000010",
            smiles="[H]OC(=O)c1c([H])c([H])c(O[H])c(O[H])c1[H]",
            inchi="InChI=1S/C7H6O4/c8-5-2-1-4(7(10)11)3-6(5)9/h1-3,8-9H,(H,10,11)",
            inchikey="YQUVCSBJEUQKSH-UHFFFAOYSA-N"
        ),
        "MNXM23": rpCompound(
            id="MNXM23",
            formula="C3H3O3",
            smiles="CC(=O)C(=O)O]",
            inchi="InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)",
            inchikey="LCTONWCANYUPML-UHFFFAOYSA-N",
            name="pyruvate"
        ),
        "CMPD_0000000025": rpCompound(
            id="CMPD_0000000025",
            smiles="[H]OC(=O)c1c([H])c([H])c([H])c(O[H])c1[H]",
            inchi="InChI=1S/C7H6O3/c8-6-3-1-2-5(4-6)7(9)10/h1-4,8H,(H,9,10)",
            inchikey="IJFXRHURBJZNAO-UHFFFAOYSA-N"
        ),
        "CMPD_0000000003": rpCompound(
            id="CMPD_0000000003",
            smiles="[H]Oc1c([H])c([H])c([H])c([H])c1O[H]",
            inchi="InChI=1S/C6H6O2/c7-5-3-1-2-4-6(5)8/h1-4,7-8H",
            inchikey="YCIMNLLNPGFGHC-UHFFFAOYSA-N"
        ),
        "MNXM337": rpCompound(
            id="MNXM337",
            smiles="[H]OC(=O)C(OC1([H])C([H])=C(C(=O)O[H])C([H])=C([H])C1([H])O[H])=C([H])[H]",
            inchi="InChI=1S/C10H10O6/c1-5(9(12)13)16-8-4-6(10(14)15)2-3-7(8)11/h2-4,7-8,11H,1H2,(H,12,13)(H,14,15)",
            inchikey="WTFXTQVDAKGDEY-UHFFFAOYSA-N"
        ),
        "MNXM2": rpCompound(
            id="MNXM2",
            smiles="[H]O[H]",
            inchi="InChI=1S/H2O/h1H2",
            inchikey="XLYOFNOQVPJJNP-UHFFFAOYSA-N"
        ),
        "MNXM13": rpCompound(
            id="MNXM13",
            smiles="O=C=O",
            inchi="InChI=1S/CO2/c2-1-3",
            inchikey="CURLTUGMZLYLDI-UHFFFAOYSA-N",
            formula="CO2",
            name="CO2"
        ),
        "MNXM5": rpCompound(
            id="MNXM5",
            smiles="N=C(O)c1ccc[n+](C2OC(COP(=O)(O)OP(=O)(O)OCC3OC(n4cnc5c(N)ncnc54)C(OP(=O)(O)O)C3O)C(O)C2O)c1",
            inchi="InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11,13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)/p+1",
            inchikey="XJLXINKUBYWONI-UHFFFAOYSA-O",
            formula="C21H25N7O17P3",
            name="NADP(+)"
        ),
        "MNXM4": rpCompound(
            id="MNXM4",
            smiles="O=O",
            inchi="InChI=1S/O2/c1-2",
            inchikey="MYMOFIZGZYHOMD-UHFFFAOYSA-N"
        ),
        "MNXM1": rpCompound(
            id="MNXM1",
            smiles="[H+]",
            inchi="InChI=1S/p+1",
            inchikey="GPRLSGONYQIRFK-UHFFFAOYSA-N"
        ),
        "MNXM6": rpCompound(
            id="MNXM6",
            smiles="[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]",
            inchi="InChI=1S/C21H30N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1,3-4,7-8,10-11,13-16,20-21,29-31H,2,5-6H2,(H2,23,32)(H,36,37)(H,38,39)(H2,22,24,25)(H2,33,34,35)",
            inchikey="ACFIXJIJDZMPPO-UHFFFAOYSA-N"
        )
    }

    def setUp(self):
        self.reactants = {
            "CMPD_0000000003": 1,
            "MNXM4": 1
        }
        self.products = {
            "TARGET_0000000001": 1,
            "MNXM1": 2
        }
        self.rxn = rpReaction(
            id="rxn_4",
            ec_numbers=[
                "1.13.11.1"
            ],
            reactants=self.reactants,
            products=self.products
        )
        self.reactions = [
            self.rxn,
            rpReaction(
                id="rxn_3",
                ec_numbers=[
                    "4.1.1.63"
                ],
                reactants={
                    "CMPD_0000000010": 1,
                    "MNXM1": 1
                },
                products={
                    "CMPD_0000000003": 1,
                    "MNXM13": 1
                }
            ),
            rpReaction(
                id="rxn_2",
                ec_numbers=[
                    "1.14.13.23"
                ],
                reactants={
                    "CMPD_0000000025": 1,
                    "MNXM4": 1,
                    "MNXM6": 1,
                    "MNXM1": 1
                },
                products={
                    "CMPD_0000000010": 1,
                    "MNXM2": 1,
                    "MNXM5": 1
                }
            ),
            rpReaction(
                id="rxn_1",
                ec_numbers=[
                    "4.1.3.45"
                ],
                reactants={
                    "MNXM337": 1
                },
                products={
                    "CMPD_0000000025": 1,
                    "MNXM23": 1
                }
            )
        ]
        self.fba = 0.57290585662576
        self.thermo = {
            "thermo_dG0_prime": {
                "value": -884.6296371353768,
                "error": 9.819227446307337,
                "units": "kilojoule / mole"
            },
            "thermo_dGm_prime": {
                "value": -884.6296371353768,
                "error": 9.819227446307337,
                "units": "kilojoule / mole"
            },
            "thermo_dG_prime": {
                "value": -884.6296371353768,
                "error": 9.819227446307337,
                "units": "kilojoule / mole"
            }
        }
        self.rp2_transfo_id = 'TRS_0_0_0'
        self.rule_id = 'RR-02-a0cc0be463ff412f-16-F'
        self.tmpl_rxn_id = 'MNXR96458'
        self.idx_in_path = 1
        self.rule_score = 0.5982208769718989
        self.id = 'pathway'
        cache = rrCache(
            db='file',
            attrs=[
                'comp_xref',
                'deprecatedCompID_compid',
            ]
        )
        self.compartments = [
            {
                'id': 'MNXC3',
                'name': 'cytosol',
                'annot': cache.get('comp_xref')[
                    cache.get('deprecatedCompID_compid')[
                        'MNXC3'
                    ]
                ]
            }
        ]

        self.parameters = {
            "upper_flux_bound": {
                "value": 999999.0,
                "units": "mmol_per_gDW_per_hr"
            },
            "lower_flux_bound": {
                "value": 0.0,
                "units": "mmol_per_gDW_per_hr"
            }
        }
        self.unit_def = {
            "mmol_per_gDW_per_hr": [
                {
                    "kind": 23,
                    "exponent": 1,
                    "scale": -3,
                    "multiplier": 1.0
                },
                {
                    "kind": 8,
                    "exponent": 1,
                    "scale": 0,
                    "multiplier": 1.0
                },
                {
                    "kind": 28,
                    "exponent": 1,
                    "scale": 0,
                    "multiplier": 3600.0
                }
            ],
            "kj_per_mol": [
                {
                    "kind": 13,
                    "exponent": 1,
                    "scale": 3,
                    "multiplier": 1.0
                },
                {
                    "kind": 13,
                    "exponent": -1,
                    "scale": 1,
                    "multiplier": 1.0
                }
            ]
        }
        self.pathway = rpPathway(
            id=self.id,
        )
        self.rxn.set_rp2_transfo_id(self.rp2_transfo_id)
        self.rxn.set_rule_id(self.rule_id)
        self.rxn.set_tmpl_rxn_id(self.tmpl_rxn_id)
        self.rxn.set_idx_in_path(self.idx_in_path)
        self.rxn.set_rule_score(self.rule_score)
        self.pathway.add_reaction(rxn=self.rxn, target_id='TARGET_0000000001')
        for rxn in self.reactions[1:]:
            self.pathway.add_reaction(rxn)
        self.sink = ['MNXM23', 'MNXM6', 'MNXM13']
        self.pathway.set_sink(self.sink)
        for key, value in self.thermo.items():
            self.pathway.set_thermo_info(key, value)
        self.pathway.set_fba_fraction(self.fba)

    ## READ METHODS
    def test_get_sink(self):
        self.assertListEqual(
            sorted(self.pathway.get_sink()),
            sorted(self.sink)
        )

    def test_get_fba(self):
        self.assertEqual(
            self.pathway.get_fba_fraction(),
            self.fba
        )

    def test_get_thermo(self):
        self.assertDictEqual(
            self.pathway.get_thermo(),
            self.thermo
        )

    def test_rpSBML(self):
        self.assertEqual(
            self.pathway,
            rpPathway.from_rpSBML(rpsbml=self.pathway.to_rpSBML())
        )

    # def test_get_compartments(self):
    #     self.assertDictEqual(
    #         self.pathway.get_compartments(),
    #         self.compartments
    #     )

    # def test_get_compartment(self):
    #     self.assertDictEqual(
    #         self.pathway.get_compartment('MNXC3'),
    #         self.compartments['MNXC3']
    #     )

    # def test_get_compartment_wrong_id(self):
    #     self.assertEqual(
    #         self.pathway.get_compartment('WRONG_ID'),
    #         None
    #     )

    # def test_get_parameters(self):
    #     self.assertDictEqual(
    #         self.pathway.get_parameters(),
    #         self.parameters
    #     )

    # def test_get_parameter(self):
    #     self.assertDictEqual(
    #         self.pathway.get_parameter('upper_flux_bound'),
    #         self.parameters['upper_flux_bound']
    #     )

    # def test_get_parameter_wrong_ID(self):
    #     self.assertEqual(
    #         self.pathway.get_parameter('WRONG_ID'),
    #         None
    #     )

    # def test_get_unit_def(self):
    #     self.assertDictEqual(
    #         self.pathway.get_unit_def(),
    #         self.unit_def
    #     )

    # def test_get_unit_def_id(self):
    #     id = 'mmol_per_gDW_per_hr'
    #     self.assertListEqual(
    #         self.pathway.get_unit_def(id),
    #         self.unit_def[id]
    #     )

    # def test_get_unit_def_wrong_id(self):
    #     id = 'WRONG_ID'
    #     self.assertEqual(
    #         self.pathway.get_unit_def(id),
    #         None
    #     )

    def test_set_id(self):
        name = 'test pathway'
        self.pathway.set_id(name)
        self.assertEqual(
            self.pathway.get_id(),
            name
        )

    def test_add_reaction(self):
        rxn = rpReaction(id='rxn')
        self.pathway.add_reaction(rxn)
        self.assertListEqual(
            self.pathway.get_list_of_reactions(),
            self.reactions + [rxn]
        )
    
    # def test_add_compartment(self):
    #     name = 'MNXC4'
    #     compartment = {
    #         "name": name,
    #         "annot": name
    #     }
    #     self.pathway.add_compartment(name, compartment)
    #     self.assertDictEqual(
    #         self.pathway.get_compartment(name),
    #         compartment
    #     )

    # def test_add_parameter(self):
    #     name = 'new_param'
    #     param = {
    #         "value": 999999.0,
    #         "units": "mmol_per_gDW_per_hr"
    #     }
    #     self.pathway.add_parameter(name, param)
    #     # print(self.pathway.get_parameter(name))
    #     self.assertDictEqual(
    #         self.pathway.get_parameter(name),
    #         param
    #     )

    # def test_add_unit_def(self):
    #     name = 'mmol_per_gDW_per_hr'
    #     unit_def = {
    #         "kind": 23,
    #         "exponent": 1,
    #         "scale": -3,
    #         "multiplier": 1.0
    #     }
    #     self.pathway.add_unit_def(name, unit_def)
    #     self.assertDictEqual(
    #         self.pathway.get_unit_def(name),
    #         unit_def
    #     )

    # def test__to_dict(self):
    #     reactions = {rxn.get_id():rxn._to_dict() for rxn in self.reactions}
    #     species = {spe.get_id():spe._to_dict() for rxn in self.reactions for spe in rxn.get_species_compounds()}
    #     pathway_dict = {
    #             **{
    #                 'id': self.id,
    #                 'reactions': reactions,
    #                 'species': species,
    #                 'sink': self.sink,
    #                 'target': self.species[list(self.species.keys())[0]].get_id(),
    #                 'parameters': self.parameters,
    #                 'thermo_dG_prime': 
    #                 'rpsbml_infos': {
    #                     'compartments': self.compartments,
    #                     'unit_def': self.unit_def,
                        
    #                 }
    #             },
    #             **self.pathway.get_fba(),
    #             **self.thermo
    #         }
    #     for key in self.pathway._to_dict().keys():
    #         if self.pathway._to_dict().get(key) != pathway_dict.get(key):
    #             if isinstance(self.pathway._to_dict().get(key), dict):
    #                 for sub_key in self.pathway._to_dict().get(key):
    #                     print(key, sub_key)
    #                     print(self.pathway._to_dict())
    #                     if self.pathway._to_dict().get(key).get(sub_key) != pathway_dict.get(key).get(sub_key):
    #                         print(sub_key)
    #                         print(self.pathway._to_dict().get(key).get(sub_key))
    #                         print(pathway_dict.get(key).get(sub_key))
    #             else:
    #                 print(self.pathway._to_dict().get(key))
    #                 print(pathway_dict.get(key))
    #             exit()
    #     self.assertDictEqual(
    #         self.pathway._to_dict(),
    #         pathway_dict
    #     )

    def test_eq(self):
        pathway = rpPathway(id='test_pathway')
        pathway.add_reaction(rxn=self.rxn, target_id='TARGET_0000000001')
        for rxn in self.reactions[1:]:
            pathway.add_reaction(rxn)
        self.assertEqual(
            self.pathway,
            pathway
        )

    def test_eq_wrong_type(self):
        self.assertNotEqual(
            self.pathway,
            0
        )

