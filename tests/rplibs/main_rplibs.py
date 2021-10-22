from unittest  import TestCase
from os        import path as os_path
from shutil    import rmtree
from brs_utils import (
    create_logger,
)
from tempfile  import mkdtemp

class Main_rplibs(TestCase):

    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    #e_coli_model_path_gz = os_path.join(
    #    data_path,
    #    'e_coli_model.sbml.gz'
    #)
    rpsbml_lycopene_path = os_path.join(
        data_path,
        'lycopene.xml'
    )
    #merged_path_gz = os_path.join(
    #    data_path,
    #    'merged_sbml.xml.gz'
    #)

    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')

        # Create persistent temp folder
        # to deflate compressed data file so that
        # it remains reachable outside of this method.
        # Has to remove manually it in tearDown() method 

        #self.temp_d = mkdtemp()
        #self.e_coli_model_path = extract_gz(
        #    self.e_coli_model_path_gz,
        #    self.temp_d
        #)
        #self.merged_path = extract_gz(
        #    self.merged_path_gz,
        #    self.temp_d
        #)

    def tearDown(self):
        pass
        #rmtree(self.temp_d)
