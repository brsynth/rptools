"""
Created on July 22 2021

@author: Joan Hérisson
"""

from unittest import TestCase
from copy import deepcopy
from rptools.rplibs.rpObject import rpObject


class Test_rpSBML(TestCase):

    __fba = {'value': 1, 'error': 0.1}
    __thermo = {'value': 2, 'error': 0.01}
    __key = 'test'

    def setUp(self):
        self.rpobject_empty = rpObject()
        self.rpobject = rpObject()
        self.rpobject.set_fba_info(self.__key, self.__fba)
        self.rpobject.set_thermo_info(self.__key, self.__thermo)

    def test__to_dict_empty(self):
        self.assertDictEqual(
            self.rpobject_empty._to_dict(),
            {}
        )

    def test__to_dict(self):
        self.assertDictEqual(
            self.rpobject._to_dict(),
            {
                f'{self.rpobject.get_fba_prefix()}_{self.__key}': self.__fba,
                f'{self.rpobject.get_thermo_prefix()}_{self.__key}': self.__thermo,
            }
        )

    def test___eq__(self):
        rpobject = deepcopy(self.rpobject)
        self.assertEqual(
            self.rpobject,
            rpobject
        )

    def test_not__eq__(self):
        # test type
        self.assertNotEqual(
            self.rpobject,
            1
        )
        # test with empty rpObject
        self.assertNotEqual(
            self.rpobject,
            self.rpobject_empty
        )
        # test by modifying a value
        rpobject = deepcopy(self.rpobject)
        fba = deepcopy(self.__fba)
        fba['value'] += 1
        rpobject.set_fba_info(self.__key, fba)
        self.assertNotEqual(
            self.rpobject,
            rpobject
        )
        # test by modifying a key
        rpobject = deepcopy(self.rpobject)
        rpobject.set_fba_info(self.__key*2, fba)
        self.assertNotEqual(
            self.rpobject,
            rpobject
        )

    def test_get_thermo(self):
        self.assertDictEqual(
            self.rpobject.get_thermo(),
            {self.__key: self.__thermo}
        )

    def test_get_thermo_info(self):
        self.assertDictEqual(
            self.rpobject.get_thermo_info(self.__key),
            self.__thermo
        )

    def test_get_fba_info(self):
        self.assertDictEqual(
            self.rpobject.get_fba_info(self.__key),
            self.__fba
        )

    def test_get_fba(self):
        self.assertDictEqual(
            self.rpobject.get_fba(),
            {self.__key: self.__fba}
        )

    def test_get_thermo_attr(self):
        value = 1.2
        attrs = [
            'dG0_prime',
            'dGm_prime',
            'dG_prime',
            'dG',
        ]
        for attr in attrs:
            with self.subTest(f'Testing get/set thermo_{attr}()', attr=attr):
                getattr(self.rpobject_empty, f'set_thermo_{attr}')(value)
                self.assertEqual(
                    getattr(self.rpobject_empty, f'get_thermo_{attr}')(),
                    value
                )

    def test_get_fba_attr(self):
        value = 1.2
        attrs = [
            'biomass',
            'fraction',
            'fba',
            'pfba',
        ]
        for attr in attrs:
            with self.subTest(f'Testing get/set fba_{attr}()', attr=attr):
                getattr(self.rpobject_empty, f'set_fba_{attr}')(value)
                self.assertEqual(
                    getattr(self.rpobject_empty, f'get_fba_{attr}')(),
                    value
                )
