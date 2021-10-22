"""
Created on June 17 2020

@author: Joan Hérisson
"""

from unittest       import TestCase
from rptools.rplibs import rpSBML, rpGraph
from os             import path as os_path
from pathlib        import Path
from main_rplibs import Main_rplibs


class Test_rpGraph(Main_rplibs):

    __test__ = False

    def setUp(self):
        super().setUp()
        self.rpsbml  = rpSBML(
            inFile = self.rpsbml_path,
            logger = self.logger
        )
        self.rpgraph = rpGraph(
            rpsbml = self.rpsbml,
            logger = self.logger
        )


    def test_onlyConsumedSpecies(self):
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, True),
                              ['MNXM89557__64__MNXC3',
                               'MNXM1__64__MNXC3',
                               'MNXM6__64__MNXC3',
                               'MNXM3__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, False),
                              ['MNXM89557__64__MNXC3', 'MNXM1__64__MNXC3'])


    #onlyProducedSpecies
    def test_onlyProducedSpecies(self):
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, True),
                              ['TARGET_0000000001__64__MNXC3',
                              'MNXM9__64__MNXC3',
                              'MNXM5__64__MNXC3',
                              'MNXM7__64__MNXC3',
                              'MNXM20__64__MNXC3',
                              'MNXM13__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, False),
                              ['TARGET_0000000001__64__MNXC3'])
