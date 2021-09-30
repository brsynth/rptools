#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle html generation."""
 
__author__ = 'Thomas Duigou'
__license__ = 'MIT'


import os
import pandas
import distutils.dir_util
from pathlib import Path

from bs4 import BeautifulSoup


class Viewer(object):
    """Viewer factory."""

    def __init__(self, out_folder='viewer'):
        """Initialising."""
        # Default configuration
        self.out_folder = Path(out_folder).resolve()
        self.template_folder = Path(__file__).resolve().parent / 'templates'
        self.template_html_file = self.template_folder / 'template.html'

        # These are not expected to be changed by the user
        self.json_file = self.out_folder / 'network_elements.js'
        self.html_file = self.out_folder / 'index.html'

    def copy_templates(self):
        """Copy the complete template tree

        :return: None
        """
        distutils.dir_util.copy_tree(
            str(self.template_folder),
            str(self.out_folder)
            )

    def write_json_deprecated(self, dict_paths, scores):
        """
        Write the html file, according the html template.

        Not needed anymore as HTML edits are made by JS

        :param dict_paths: dict of pathways
        :param scores: {'score_type_1': {'pathway_1': int, 'pathway_2': int, ...}, ...}
        :return: None
        """
        # Get template content
        with open(self.template_html_file) as ifh:
            soup = BeautifulSoup(ifh, 'html.parser')

        # Update selectbox with scores
        select_script = soup.find(id="selectbox")
        for i in scores:
            new_tag = soup.new_tag("option")
            new_tag["value"] = i
            new_tag.append(i)
            select_script.append(new_tag)

        # Pathway table
        df = pandas.DataFrame()
        pathway_name_list = []
        select_input_list = []
        score_list = []
        for pathway_name in dict_paths:
            pathway_name_list.append(pathway_name)
            select_input_list.append('')
            score_list.append('')
        df['Pathway'] = pathway_name_list
        df['Select'] = select_input_list
        df['Score'] = score_list

        # Append the pathway table
        html_str = df.to_html(index=False)
        html_str = html_str.replace('&lt;', '<').replace('&gt;', '>')
        table = soup.find(id='table_path')
        table.append(BeautifulSoup(html_str, 'html.parser'))

        # Write
        html = soup.prettify('utf-8')
        with open(self.html_file, 'wb') as ofh:
            ofh.write(html)
