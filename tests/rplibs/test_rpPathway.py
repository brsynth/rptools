"""
Created on May 28 2021

@author: Joan Hérisson
"""

from os import path as os_path
from tempfile import NamedTemporaryFile
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
        self.target = rpCompound(
            id='TARGET_0000000001',
            smiles='[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]'
        )
        self.reactants = {
            "CMPD_0000000003": 1,
            "MNXM4": 1
        }
        self.products = {
            self.target.get_id(): 1,
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
            "dG0_prime": {
                "value": -884.6296371353768,
                "error": 9.819227446307337,
                "units": "kilojoule / mole"
            },
            "dGm_prime": {
                "value": -884.6296371353768,
                "error": 9.819227446307337,
                "units": "kilojoule / mole"
            },
            "dG_prime": {
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
        self.pathway.set_parameters(self.parameters)
        self.pathway.set_unit_defs(self.unit_def)
        self.rxn.set_rp2_transfo_id(self.rp2_transfo_id)
        self.rxn.set_rule_id(self.rule_id)
        self.rxn.set_tmpl_rxn_id(self.tmpl_rxn_id)
        self.rxn.set_idx_in_path(self.idx_in_path)
        self.rxn.set_rule_score(self.rule_score)
        self.pathway.add_reaction(rxn=self.rxn, target_id=self.target.get_id())
        for rxn in self.reactions[1:]:
            self.pathway.add_reaction(rxn)
        self.sink = ['MNXM23', 'MNXM6', 'MNXM13']
        self.pathway.set_sink(self.sink)
        for key, value in self.thermo.items():
            self.pathway.set_thermo_info(key, value)
        self.pathway.set_fba_fraction(self.fba)

    ## READ METHODS
    def test__to_dict(self):
        self.assertDictEqual(
            self.pathway._to_dict(specific=True),
            {
                'sink': self.pathway.get_sink(),
                'target': self.pathway.get_target_id(),
                'parameters': self.pathway.get_parameters(),
                'unit_defs': self.pathway.get_unit_defs(),
                'compartments': self.pathway.get_compartments(),
                'fba_fraction': self.pathway.get_fba_fraction(),
                'thermo_dG0_prime': self.pathway.get_thermo_dG0_prime(),
                'thermo_dGm_prime': self.pathway.get_thermo_dGm_prime(),
                'thermo_dG_prime': self.pathway.get_thermo_dG_prime(),
            }
        )

    def test_completed_species(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.set_completed_species([species[0]])
        self.pathway.add_completed_species(species)
        self.assertCountEqual(
            self.pathway.get_completed_species(),
            species
        )

    def test_add_species_group_new(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.add_species_group(
            'completed',
            species
        )
        self.assertCountEqual(
            self.pathway.get_completed_species(),
            species
        )

    def test_add_species_group(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.set_completed_species([species[0]])
        self.pathway.add_species_group(
            'completed',
            species
        )
        self.assertCountEqual(
            self.pathway.get_completed_species(),
            species
        )

    def test___set_species_group_wrong_type(self):
        species = {'SPE_1': 1, 'SPE_2': 2}
        self.pathway.add_species_group(
            'completed',
            species
        )
        self.assertListEqual(
            sorted(self.pathway.get_completed_species()),
            list(species.keys())
        )

    def test_fba_ignored_species(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.set_fba_ignored_species(species)
        self.assertCountEqual(
            self.pathway.get_fba_ignored_species(),
            species
        )

    def test_thermo_substituted_species(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.set_thermo_substituted_species(species)
        self.assertCountEqual(
            self.pathway.get_thermo_substituted_species(),
            species
        )

    def test_intermediate_species(self):
        species = ['SPE_1', 'SPE_2']
        self.pathway.set_trunk_species(
            species
            + self.sink[:1]
        )
        self.pathway.add_trunk_species([self.pathway.get_target_id()])
        self.assertCountEqual(
            self.pathway.get_intermediate_species(),
            species
        )

    def test_get_sink(self):
        self.assertCountEqual(
            self.pathway.get_sink(),
            self.sink
        )

    def test_get_target_rxn_id(self):
        self.assertEqual(
            self.pathway.get_target_rxn_id(),
            self.rxn.get_id()
        )

    def test_get_rxn_target(self):
        self.assertEqual(
            self.pathway.get_rxn_target(),
            self.rxn
        )

    def test_get_target(self):
        self.assertEqual(
            self.pathway.get_target(),
            self.target
        )

    def test_get_parameters(self):
        self.assertDictEqual(
            self.pathway.get_parameters(),
            self.parameters
        )

    def test_get_unit_def(self):
        unit = 'mmol_per_gDW_per_hr'
        self.assertListEqual(
            self.pathway.get_unit_def(unit),
            self.unit_def[unit]
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

    def test_rpSBML_rpsbml(self):
        # print(self.pathway.get_species_ids())
        # print(self.pathway)
        self.assertEqual(
            self.pathway,
            rpPathway.from_rpSBML(rpsbml=self.pathway.to_rpSBML())
        )

    def test_rpSBML_file(self):
        with NamedTemporaryFile() as tempf:
            self.pathway.to_rpSBML().write_to_file(tempf.name)
            self.assertEqual(
                self.pathway,
                rpPathway.from_rpSBML(
                    infile=tempf.name
                )
            )

    def test_rpSBML_file_rpsbml(self):
        with NamedTemporaryFile() as tempf:
            self.pathway.to_rpSBML().write_to_file(tempf.name)
            self.assertEqual(
                self.pathway,
                rpPathway.from_rpSBML(
                    infile=tempf.name,
                    rpsbml=None
                )
            )

    def test_rename_compound(self):
        for spe_id in ['CMPD_0000000003', 'TARGET_0000000001']:
            with self.subTest("Message for this subtest", spe_id=spe_id):
                species = deepcopy(self.pathway.get_species_ids())
                new_spe_id = 'NEW'+spe_id
                self.pathway.rename_compound(
                    spe_id,
                    new_spe_id
                )
                self.assertCountEqual(
                    list(set(self.pathway.get_species_ids())),
                    list(
                        (
                            set(species)
                            - {spe_id}
                        )
                        | {new_spe_id}
                    )
                )

    def test_cobraize_uncobraize(self):
        from rptools.rpfba.cobra_format import at_pattern
        compartment = 'COMP'
        species = deepcopy(self.pathway.get_species_ids())
        self.pathway.cobraize(compartment)
        self.assertCountEqual(
            self.pathway.get_species_ids(),
            [spe_id+at_pattern+compartment for spe_id in species]
        )
        self.pathway.uncobraize()
        self.assertCountEqual(
            self.pathway.get_species_ids(),
            species
        )

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
    
    def test_add_parameter(self):
        name = 'upper_flux_bound'
        params = {
            "value": 0,
            "units": "other"
        }
        self.pathway.add_parameter(name, params)
        self.assertDictEqual(
            self.pathway.get_parameter(name),
            self.parameters[name]
        )

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

