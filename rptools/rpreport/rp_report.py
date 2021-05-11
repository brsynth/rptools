#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""rp_report.py: Parse SBML files and generate a web page with relevant pathways data."""

__author__ = 'Olivier TELLE - INRAE'
__license__ = 'MIT'

import os
import sys
import json
import logging
import tarfile
import tempfile
import fnmatch
from shutil import copyfile
from pathlib import Path
from rptools.rplibs import rpSBML
from rptools.rpreport.dictor.dictor import dictor


def get_reactions_data(rxn_dict):
    """Extract, sort and return a dictionary of reactions data

    :param rxn_dict: data of 'reactions' element of a pathway dictionary
    :type rxn_dict: dict
    :return: return relevant and sorted reactions data
    :rtype: dict
    """

    if not isinstance(rxn_dict, dict):
        raise AttributeError('get_reactions_data() expects dict as argument.')

    # init
    reaction_dict = {}
    reaction = {}

    for rxn_name in rxn_dict.keys():
        # We store step number of the reaction
        reaction_dict['rxn_idx'] = dictor(rxn_dict, f"{rxn_name}.brsynth.rxn_idx")

        # We store all the ec-codes
        reaction_dict['ec_code'] = dictor(rxn_dict, f"{rxn_name}.miriam.ec-code")

        # We store the dfg_prime
        reaction_dict['dfG_prime_m'] = dictor(rxn_dict, f"{rxn_name}.brsynth.dfG_prime_m.value")

        # We store the rule_score
        reaction_dict['rule_score'] = dictor(rxn_dict, f"{rxn_name}.brsynth.rule_score")

        reaction[rxn_name] = reaction_dict
        reaction_dict = {}  # emptying is useless?

    # sorting dict by rxn_idx by reinserting values
    reaction = dict(sorted(reaction.items(), key=lambda item: item[1]['rxn_idx']))

    return reaction

def to_data_js(sbml_files: list, source_path: str, output_folder: str, verbose=False, dev=False):
    """
    Return a list of dictionaries parsed from sbml files
    """
    # creating list where all necessary elements will be compiled
    rp_list = []

    # loop and operations for each sbml.xml files found
    for name in sbml_files:
        if verbose:
            print("Parsing", name)

        rpsbml = rpSBML(inFile=os.path.join(source_path, name))

        pathway_dict = rpsbml.toDict()

        # if pathway name is found
        if dictor(pathway_dict, 'pathway.brsynth.path_id', default=False):
            rp_name = dictor(pathway_dict, "pathway.brsynth.path_id")
            if verbose:
                print("Path_id found:", rp_name)

            # adding necessary values to the list
            rp_list.append({
                'pathway_name': rp_name,
                'dfG_prime_m': dictor(pathway_dict, "pathway.brsynth.dfG_prime_m.value"),
                'global_score': dictor(pathway_dict, "pathway.brsynth.global_score"),
                'fba_obj_fraction': dictor(pathway_dict, "pathway.brsynth.fba_obj_fraction.value"),
                'norm_rule_score': dictor(pathway_dict, "pathway.brsynth.norm_rule_score"),
                'nb_reactions': dictor(pathway_dict, "pathway.brsynth.nb_reactions"),
                'reactions': get_reactions_data(pathway_dict['reactions'])
            })

            # sorting list by pathway_name (the 1st element)
            rp_list = sorted(rp_list, key=lambda k: k['pathway_name'])

            if dev:
                #  Saving pathway_dict into separate json file
                with open(output_folder + '/dev/' + rp_name + '.json', "w") as f:
                    json.dump(pathway_dict, f, indent=4)
        elif verbose:
            print("No path_id found, file ignored!")

    return rp_list


def write_to_one_html(templates_dir, data):
    """
            Write js, css and data files into a standalone html file
    """

    # Load the html template file into memory
    with open(templates_dir + '/index.html', 'r') as file:
        standalone_html_file = file.read()
    with open(templates_dir + '/js/ag-grid-community.min.noStyle.js', 'r') as file:
        ag_grid_community_min_nostyle_js = file.read()
    with open(templates_dir + '/js/ag-charts-community.min.js', 'r') as file:
        ag_charts_community_min_js = file.read()
    with open(templates_dir + '/js/main.js', 'r') as file:
        main_js = file.read()
    with open(templates_dir + '/js/bootstrap.bundle.min.js', 'r') as file:
        bootstrap_bundle_min_js = file.read()
    with open(templates_dir + '/css/bootstrap.min.css', 'r') as file:
        bootstrap_min_css = file.read()
    with open(templates_dir + '/css/ag-grid.min.css', 'r') as file:
        ag_grid_min_css = file.read()
    with open(templates_dir + '/css/ag-theme-alpine.min.css', 'r') as file:
        ag_theme_alpine_min_css = file.read()

    # Replace the link to data.js with its content itself
    standalone_html_file = standalone_html_file.replace('<script src="./js/data.js">',
                                                        '<script>' + data)
    standalone_html_file = standalone_html_file.replace('<script src="./js/ag-grid-community.min.noStyle.js">',
                                                        '<script>' + ag_grid_community_min_nostyle_js)
    standalone_html_file = standalone_html_file.replace('<script src="./js/ag-charts-community.min.js">',
                                                        '<script>' + ag_charts_community_min_js)
    standalone_html_file = standalone_html_file.replace('<script src="./js/main.js">',
                                                        '<script>' + main_js)
    standalone_html_file = standalone_html_file.replace('<script src="./js/bootstrap.bundle.min.js">',
                                                        '<script>' + bootstrap_bundle_min_js)

    # Replace css > inline
    standalone_html_file = standalone_html_file.replace('<link rel="stylesheet" href="./css/bootstrap.min.css">',
                                                        '<style>' + bootstrap_min_css + '</style>')
    standalone_html_file = standalone_html_file.replace('<link rel="stylesheet" href="./css/ag-grid.min.css">',
                                                        '<style>' + ag_grid_min_css + '</style>')
    standalone_html_file = standalone_html_file.replace('<link rel="stylesheet" href="./css/ag-theme-alpine.min.css">',
                                                        '<style>' + ag_theme_alpine_min_css + '</style>')

    return standalone_html_file


def run_report(input_dir:bool, source_path:str, output_folder:str, dev:bool, verbose:bool, standalone:bool):
    """
        Converting SBML RP files into web report.
    """

    # Warning & errors logging
    logging.basicConfig(stream=sys.stderr,
                        level=logging.WARNING,
                        datefmt='%d/%m/%Y %H:%M:%S',
                        format='%(asctime)s -- %(levelname)s -- %(message)s')

    # Make output folder(s) if needed
    if not os.path.isfile(output_folder):
        try:
            os.makedirs(output_folder, exist_ok=True)
        except IOError as e:
            raise e
    if not standalone:  # no need if it's a standalone html file
        if not os.path.isfile(output_folder + '/js'):
            try:
                os.makedirs(output_folder + '/js', exist_ok=True)
            except IOError as e:
                raise e
        if not os.path.isfile(output_folder + '/css'):
            try:
                os.makedirs(output_folder + '/css', exist_ok=True)
            except IOError as e:
                raise e
    if dev:
        if not os.path.isfile(output_folder + '/dev'):
            try:
                os.makedirs(output_folder + '/dev', exist_ok=True)
            except IOError as e:
                raise e

    # if -d option exists then parse files in the directory
    if input_dir:
        files = fnmatch.filter(os.listdir(source_path), "*.xml")
        rp_list = to_data_js(files, source_path, output_folder, verbose, dev)
    else:
        if not os.path.isfile(source_path):
            logging.error('File "{}" not found, exit'.format(source_path))
            sys.exit(1)
        if tarfile.is_tarfile(source_path):
            with tempfile.TemporaryDirectory() as tmp_folder:
                tar = tarfile.open(source_path, mode='r')
                tar.extractall(path=tmp_folder)
                tar.close()
                files = os.listdir(tmp_folder)

                rp_list = to_data_js(files, tmp_folder, output_folder, verbose, dev)
        else:
            rp_list = ''

    if dev:
        # Saving rp_list into json file just for dev purpose.
        with open(os.path.join(output_folder, 'dev/digest.json'), "w") as f:
            json.dump(rp_list, f, indent=4)

    # Saving data into variable
    data = json.dumps(rp_list)
    data = "const dataJSON = " + data

    templates_dir = str(Path(__file__).resolve().parent / 'templates')

    # In case its not a standalone file, make dirs and copy css and js files
    if standalone:
        standalone_file = write_to_one_html(templates_dir, data)
        with open(output_folder + '/index.html', "w") as f:
            f.write(standalone_file)
            f.close()
    else:
        copyfile(templates_dir + '/index.html', output_folder + '/index.html')
        if not os.path.isfile(output_folder + '/js/ag-charts-community.min.js'):
            copyfile(templates_dir + '/js/ag-charts-community.min.js',
                            output_folder + '/js/ag-charts-community.min.js')
        if not os.path.isfile(output_folder + '/js/ag-grid-community.min.noStyle.js'):
            copyfile(templates_dir + '/js/ag-grid-community.min.noStyle.js',
                            output_folder + '/js/ag-grid-community.min.noStyle.js')
        if not os.path.isfile(output_folder + '/js/bootstrap.bundle.min.js'):
            copyfile(templates_dir + '/js/bootstrap.bundle.min.js', output_folder + '/js/bootstrap.bundle.min.js')
        if not os.path.isfile(output_folder + '/js/bootstrap.bundle.min.js.map'):
            copyfile(templates_dir + '/js/bootstrap.bundle.min.js.map',
                            output_folder + '/js/bootstrap.bundle.min.js.map')
        if not os.path.isfile(output_folder + '/js/jquery-3.6.0.min.js'):
            copyfile(templates_dir + '/js/jquery-3.6.0.min.js', output_folder + '/js/jquery-3.6.0.min.js')
        if not os.path.isfile(output_folder + '/js/main.js'):
            copyfile(templates_dir + '/js/main.js', output_folder + '/js/main.js')
        if not os.path.isfile(output_folder + '/css/ag-grid.css'):
            copyfile(templates_dir + '/css/ag-grid.css', output_folder + '/css/ag-grid.css')
        if not os.path.isfile(output_folder + '/css/ag-grid.min.css'):
            copyfile(templates_dir + '/css/ag-grid.min.css', output_folder + '/css/ag-grid.min.css')
        if not os.path.isfile(output_folder + '/css/ag-theme-alpine.min.css'):
            copyfile(templates_dir + '/css/ag-theme-alpine.min.css',
                            output_folder + '/css/ag-theme-alpine.min.css')
        if not os.path.isfile(output_folder + '/css/bootstrap.min.css'):
            copyfile(templates_dir + '/css/bootstrap.min.css', output_folder + '/css/bootstrap.min.css')
        if not os.path.isfile(output_folder + '/css/bootstrap.min.css.map'):
            copyfile(templates_dir + '/css/bootstrap.min.css.map', output_folder + '/css/bootstrap.min.css.map')
        # Saving data as constant in js file (can't load json files without webserver)
        with open(output_folder + '/js/data.js', "w") as f:
            f.write(data)
            f.close()
