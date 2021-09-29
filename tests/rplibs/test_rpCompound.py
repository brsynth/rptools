"""
Created on July 22 2021

@author: Joan Hérisson
"""

from unittest import TestCase
from copy import deepcopy
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