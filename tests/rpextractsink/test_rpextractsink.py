"""
Created on June 17 2020

@author: Joan Hérisson
"""

# Generic for test process
from unittest import TestCase

# Specific for tool
from rptools.rpextractsink import genSink
from rr_cache import rrCache

# Specific for tests themselves
from tempfile import NamedTemporaryFile
from re import findall as re_findall
from os import (
    path as os_path,
    remove
)
from brs_utils import (
    create_logger,
    extract_gz
)
from shutil    import rmtree
from tempfile  import mkdtemp


# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
# 'Test_' prefix is mandatory
class Test_rpExtractSink(TestCase):


    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    e_coli_model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )

    cache = rrCache(
        ['cid_strc']
    )


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        # Create persistent temp folder
        # to deflate compressed data file so that
        # it remains reachable outside of this method.
        # Has to remove manually it in tearDown() method 
        self.temp_d = mkdtemp()

        self.e_coli_model_path = extract_gz(
            self.e_coli_model_path_gz,
            self.temp_d
        )


    def tearDown(self):
        rmtree(self.temp_d)


    def test_genSink(self):
        # Test with dead ends
        self._test_genSink(
            cache = self.cache,
            input_sbml = self.e_coli_model_path,
            remove_dead_end = False,
            compartment_id = 'MNXC3',
            ref_file = 'output_sink.csv'
        )
        # Test without dead ends
        self._test_genSink(
            cache = self.cache,
            input_sbml = self.e_coli_model_path,
            remove_dead_end = True,
            compartment_id = 'MNXC3',
            ref_file = 'output_sink_woDE.csv'
        )


    def _test_genSink(
        self,
        cache: rrCache,
        input_sbml: str,
        remove_dead_end: bool,
        compartment_id: str,
        ref_file: str
    ):
        test_sink = genSink(
            self.cache,
            input_sbml = input_sbml,
            remove_dead_end = remove_dead_end,
            compartment_id = compartment_id,
            logger = self.logger
        )
        ref_sink = {}
        with open(
            os_path.join(self.data_path, ref_file),
            'r'
        ) as ref_f:
            ref_content = ref_f.readlines()
            for line in ref_content[1:]: # skip header
                id, inchi = re_findall(r'"([^"]+)"', line)
                ref_sink[id] = inchi
        self.assertDictEqual(test_sink, ref_sink)
