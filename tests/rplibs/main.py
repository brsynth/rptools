from os       import path as os_path
from tests.main     import Main
from brs_utils import extract_gz


class Main_rplibs(Main):

    data_path = os_path.join(
        os_path.dirname(__file__),
        'data'
    )
    rpsbml_path = os_path.join(
        data_path,
        'rpsbml.xml'
    )
    merged_path_gz = os_path.join(
        data_path,
        'merged_sbml.xml.gz'
    )

    def setUp(self):
        super().setUp()
        self.merged_path = extract_gz(
            self.merged_path_gz,
            self.temp_d
        )

