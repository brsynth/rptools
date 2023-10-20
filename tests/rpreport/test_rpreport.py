from rptools.rpreport import run_report
import os
import fnmatch
import tempfile
import filecmp

__test__ = True

data_path = os.path.join(
        os.path.dirname(__file__),
        'data'
    )
data_input_dir_path = os.path.join(
        data_path,
        'input'
    )
data_input_tar_file = os.path.join(
        data_input_dir_path,
        'input_rpSBML.tar'
    )
data_output_dir_path = os.path.join(
        data_path,
        'output'
    )
data_output_standalone_file = os.path.join(
        data_path,
        'standalone_output',
        'index.html'
    )
data_output_js_file = os.path.join(
        data_path,
        'output',
        'js',
        'data.js'
    ) 

#files = fnmatch.filter(os.listdir(data_input_dir_path), "*.xml")
#files.sort()

def test_run_report_standalone_from_dir():
    # testing standalone file output from files into a directory
    with tempfile.TemporaryDirectory() as tmp_folder:
        run_report(True, data_input_dir_path, tmp_folder, False, False, True)
        tested_output_single_file_html = os.path.join(
            tmp_folder,
            'index.html'
        )
        with open(tested_output_single_file_html, 'r') as test_f:
            test_content = test_f.read()
            with open(data_output_standalone_file, 'r') as ref_f:
                ref_content = ref_f.read()
                assert test_content == ref_content
        # assert filecmp.cmp(tested_output_single_file_html, data_output_standalone_file, shallow=False)


def test_run_report_standalone():
    # testing standalone file output from tar file
    with tempfile.TemporaryDirectory() as tmp_folder:
        run_report(False, data_input_tar_file, tmp_folder, False, False, True)
        tested_output_single_file_html = os.path.join(
            tmp_folder,
            'index.html'
        )    
        with open(tested_output_single_file_html, 'r') as test_f:
            test_content = test_f.read()
            with open(data_output_standalone_file, 'r') as ref_f:
                ref_content = ref_f.read()
                assert test_content == ref_content
        # assert filecmp.cmp(tested_output_single_file_html, data_output_standalone_file, shallow=False)

def test_run_report():
    # testing files output from tar file
    with tempfile.TemporaryDirectory() as tmp_folder:
        run_report(False, data_input_tar_file, tmp_folder, False, False, False)
        tested_output_js_file = os.path.join(
            tmp_folder,
            'js',
            'data.js'
        ) 
        assert os.path.exists(os.path.join(data_output_dir_path, 'css', 'ag-grid.css'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'css', 'ag-grid.min.css'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'css', 'ag-theme-alpine.min.css'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'css', 'bootstrap.min.css'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'css', 'bootstrap.min.css.map'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'ag-charts-community.min.js'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'ag-grid-community.min.noStyle.js'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'bootstrap.bundle.min.js'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'bootstrap.bundle.min.js.map'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'jquery-3.6.0.min.js'))
        assert os.path.exists(os.path.join(data_output_dir_path, 'js', 'main.js'))
        assert filecmp.cmp(tested_output_js_file, data_output_js_file, shallow=False)