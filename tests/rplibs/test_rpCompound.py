"""
Created on July 22 2021

@author: Joan Hérisson
"""

from unittest import TestCase
from copy import deepcopy
from chemlite import Compound
from rptools.rplibs.rpCompound import rpCompound


class Test_rpCompound(TestCase):

    __id = 'test'
    __dict = {
        'id': __id,
        'name': '',
        'smiles': '',
        'inchi': '',
        'inchikey': '',
        'formula': ''
    }
    __dg = 1.0
    __biomass = 2.0
    __fraction = 3.0
    __fba = 4.0
    __pfba = 5.0

    def setUp(self):
        self.rpcompound_empty = rpCompound(id=f'{self.__id}_empty')
        self.rpcompound = rpCompound(id=self.__id)
        self.rpcompound.set_thermo_standard_dg_formation(self.__dg)
        self.rpcompound.set_fba_biomass_shadow_price(self.__biomass)
        self.rpcompound.set_fba_fraction_shadow_price(self.__fraction)
        self.rpcompound.set_fba_fba_shadow_price(self.__fba)
        self.rpcompound.set_fba_pfba_shadow_price(self.__pfba)

    def test__to_dict_empty(self):
        self.__dict['id'] = f'{self.__id}_empty'
        self.assertDictEqual(
            self.rpcompound_empty._to_dict(),
            self.__dict
        )

    def test__to_dict(self):
        self.assertDictEqual(
            self.rpcompound._to_dict(),
            {
                **self.__dict,
                **self.rpcompound._to_dict(full=False)
            }
        )

    def test_get_thermo_attr(self):
        attrs = [
            'standard_dg_formation',
        ]
        for attr in attrs:
            with self.subTest(f'Testing get/set thermo_{attr}()', attr=attr):
                getattr(self.rpcompound_empty, f'set_thermo_{attr}')(self.__dg)
                self.assertEqual(
                    getattr(self.rpcompound_empty, f'get_thermo_{attr}')(),
                    self.__dg
                )

    def test_get_fba_attr(self):
        attrs = [
            'biomass',
            'fraction',
            'fba',
            'pfba',
        ]
        for attr in attrs:
            with self.subTest(f'Testing get/set fba_{attr}_shadow_price()', attr=attr):
                getattr(self.rpcompound_empty, f'set_fba_{attr}_shadow_price')(self.__biomass)
                self.assertEqual(
                    getattr(self.rpcompound_empty, f'get_fba_{attr}_shadow_price')(),
                    self.__biomass
                )

    def test_compartment(self):
        comp_id = 'COMP_ID'
        self.rpcompound_empty.set_compartment(comp_id)
        self.assertEqual(
            self.rpcompound_empty.get_compartment(),
            comp_id
        )    
    def test_from_compound(self):
        # Load.
        cid = "MNXM23"
        smiles = "CC(=O)C(=O)O]"
        inchi = "InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)"
        inchikey = "LCTONWCANYUPML-UHFFFAOYSA-N"
        name = "target"
        formula = "C3H3O3"
        compound = Compound(
            id=cid,
            smiles=smiles,
            inchi=inchi,
            inchikey=inchikey,
            name=name,
            formula=formula
        )
        rp_compound = rpCompound.from_compound(
            compound=compound
        )
        # Return type.
        self.assertIsInstance(
            rp_compound,
            rpCompound
        )
        # Challenge - 1
        self.assertEqual(
            compound._Compound__to_dict(),
            rp_compound._Compound__to_dict()
        )
        # Challenge - 2
        rp_compound = rpCompound.from_compound(
            compound=compound,
            compartment_id='e'
        )
        self.assertEqual(
            rp_compound.get_compartment(),
            'e'
        )
