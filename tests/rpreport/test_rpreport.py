from rptools.rpreport.rp_report import (
    to_data_js,
    run_report
)
import os
import fnmatch
import tempfile
import filecmp

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

files = fnmatch.filter(os.listdir(data_input_dir_path), "*.xml")
files.sort()

rp_list = [{'pathway_name': 'rp_002_0001', 'dfG_prime_m': None, 'global_score': 0.6006580897054999, 'fba_obj_fraction': 0.5493228602383657, 'norm_rule_score': 0.6346154014079132, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.5869134225065102}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0021', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0022', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0023', 'dfG_prime_m': 235.0, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.31415937648615916, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': 140.0, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': 300.0, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': 235.0, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0136', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0137', 'dfG_prime_m': -35.0, 'global_score': None, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0138', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0251', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': None, 'norm_rule_score': 0.6885137648615917, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': 0.5683242688496836}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': 0.7486085128675456}}}, {'pathway_name': 'rp_003_0252', 'dfG_prime_m': None, 'global_score': 0.6242847593495533, 'fba_obj_fraction': 0.4051338393927223, 'norm_rule_score': None, 'nb_reactions': 3, 'reactions': {'rxn_1': {'rxn_idx': 1, 'ec_code': ['2.5.1.85', '2.5.1.84', '2.5.1.31', '2.5.1.81', '2.5.1.88', '2.5.1.-', '2.5.1', '2.5.1.29'], 'dfG_prime_m': None, 'rule_score': None}, 'rxn_2': {'rxn_idx': 2, 'ec_code': ['2.5.1.32', '2.5.1.99', '2.5.1.96'], 'dfG_prime_m': None, 'rule_score': None}, 'rxn_3': {'rxn_idx': 3, 'ec_code': ['1.3.99.31'], 'dfG_prime_m': None, 'rule_score': None}}}]

def test_to_data_js():
    assert to_data_js(files, data_input_dir_path, '') == rp_list

def test_run_report_standalone_from_dir():
    # testing standalone file output from files into a directory
    with tempfile.TemporaryDirectory() as tmp_folder:
        run_report(True, data_input_dir_path, tmp_folder, False, False, True)
        tested_output_single_file_html = os.path.join(
            tmp_folder,
            'index.html'
        )    
        assert filecmp.cmp(tested_output_single_file_html, data_output_standalone_file, shallow=False)


def test_run_report_standalone():
    # testing standalone file output from tar file
    with tempfile.TemporaryDirectory() as tmp_folder:
        run_report(False, data_input_tar_file, tmp_folder, False, False, True)
        tested_output_single_file_html = os.path.join(
            tmp_folder,
            'index.html'
        )    
        assert filecmp.cmp(tested_output_single_file_html, data_output_standalone_file, shallow=False)

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