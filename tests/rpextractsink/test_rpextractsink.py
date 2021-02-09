"""
Created on June 17 2020

@author: Joan Hérisson
"""

# Generic for test process
from tests.main import Main

# Specific for tool
from rptools.rpextractsink import genSink
from rptools.rplibs        import rpCache

# Specific for tests themselves
from pathlib  import Path
from tempfile import NamedTemporaryFile
from filecmp  import cmp
from os       import path as os_path
import logging


# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
# 'Test_' prefix is mandatory
class Test_rpExtractSink(Main):


    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )

    rpcache = rpCache(
        'file',
        ['cid_strc']
    )


    def setUp(self):
        super().setUp()


    def test_genSink(self):
        outfile = NamedTemporaryFile(delete=True)
        genSink(
            self.rpcache,
            input_sbml = self.e_coli_model_path,
            output_sink = outfile.name,
            remove_dead_end = False,
            logger = self.logger
        )
        self.assertTrue(
            cmp(
                Path(outfile.name),
                os_path.join(
                    self.data_path,
                    'output_sink.csv'
                )
            )
        )
        outfile.close()


    def test_genSink_rmDE(self):
        outfile = NamedTemporaryFile(delete=True)
        genSink(
            self.rpcache,
            input_sbml = self.e_coli_model_path,
            output_sink = outfile.name,
            remove_dead_end = True
        )
        self.assertTrue(
            cmp(
                Path(outfile.name),
                os_path.join(
                    self.data_path,
                    'output_sink_woDE.csv'
                )
            )
        )
        outfile.close()
