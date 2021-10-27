from unittest  import TestCase
from os        import path as os_path
from shutil    import rmtree
from tempfile  import mkdtemp

from brs_utils import (
    create_logger,
)
from rptools.rplibs.rpSBML import rpSBML

class Main_rpfba(TestCase):

    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    e_coli_model_path_gz = os_path.join(
        data_path,
        'e_coli_model.sbml.gz'
    )
    e_coli_model_path = os_path.join(
        data_path,
        'e_coli_iML1515.sbml'
    )
 
    pathway_path = os_path.join(
        data_path,
        'pathway.json'
    )
    medium_path = os_path.join(
        data_path,
        'medium'
    )

    rpsbml = None

    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')
        # Create persistent temp folder
        # to deflate compressed data file so that
        # it remains reachable outside of this method.
        # Has to remove manually it in tearDown() method
        self.temp_d = mkdtemp()
        self.rpsbml = rpSBML(
            inFile=self.e_coli_model_path,
            logger=self.logger
        )
 
    def tearDown(self):
        rmtree(self.temp_d)

