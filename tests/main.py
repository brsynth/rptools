from unittest  import TestCase
from os        import path as os_path
from shutil    import rmtree
from brs_utils import (
    create_logger,
    extract_gz
)
from tempfile  import mkdtemp


class Main(TestCase):

    common_data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    e_coli_model_path_gz = os_path.join(
        common_data_path,
        'e_coli_model.sbml.gz'
    )


    def setUp(self):
        self.logger = create_logger(__name__, 'ERROR')
        self.temp_d = mkdtemp()
        self.e_coli_model_path = extract_gz(
            self.e_coli_model_path_gz,
            self.temp_d
        )


    def tearDown(self):
        rmtree(self.temp_d)