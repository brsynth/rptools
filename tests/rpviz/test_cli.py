"""Test cases for the rpviz CLI module."""
import json
from typing import Dict
from pathlib import Path

import pytest
import deepdiff

from rptools.rpviz.__main__ import (
    __build_arg_parser,
    __run,
)


REF_IN_DIR = Path(__file__).resolve().parent / 'inputs' / 'as_dir'
REF_IN_TAR = Path(__file__).resolve().parent / 'inputs' / 'as_tar.tgz'
REF_OUT_DIR = Path(__file__).resolve().parent / 'outputs'
COF_FILE = (
    Path(__file__).resolve().parent.parent.parent
    / "rptools"
    / "rpviz"
    / "data"
    / "cofactors_mnx_202507.tsv"
)


def __dump_file(path: str, data: Dict) -> None:
    """Dump a dictionary to a file."""
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)


def __read_and_sort(path: str, skip=[]) -> list:
    with open(path, encoding='utf-8') as fh:
        lines = []
        for line in fh.readlines():
            if not any(
                skip_item in line
                for skip_item in skip
            ):
                lines.append(line.strip())
    return sorted(lines)


def __read_multi_object_json(path: str) -> Dict:
    """Read a JSON-like file containing multiple JSON objects."""
    objects = {}
    with open(path, encoding='utf-8') as fh:
        obj_name = ''
        obj_content = ''
        curly_count = 0
        for line in fh:
            for char in line.strip():
                if char == '{':
                    if curly_count == 0:
                        # Start of a new object
                        obj_name = line.strip().split()[0]
                        obj_content = ''
                    obj_content += char
                    curly_count += 1
                elif char == '}':
                    obj_content += char
                    if curly_count == 1:
                        # End of the current object
                        objects[obj_name] = json.loads(obj_content)
                        obj_name = ''
                        obj_content = ''
                    curly_count -= 1
                else:
                    obj_content += char
        # Last object
        if curly_count == 0 and obj_content:
            objects[obj_name] = json.loads(obj_content)
    return objects


def test_build_arg_parser(mocker):
    """Test the argument parser."""
    args = ['prog']  # no arguments
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args()


def test_dir_input(mocker, tmpdir):
    """Test the CLI with a directory input."""
    args = ['prog', str(REF_IN_DIR), str(tmpdir)]
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    ref_file = REF_OUT_DIR / 'network.json'
    test_file = tmpdir / 'network.json'
    ref_objects = __read_multi_object_json(ref_file)
    test_objects = __read_multi_object_json(test_file)
    assert not deepdiff.DeepDiff(ref_objects, test_objects, ignore_order=True)


def test_tar_input(mocker, tmpdir):
    """Test the CLI with a tar input."""
    args = ['prog', str(REF_IN_TAR), str(tmpdir)]
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    ref_file = REF_OUT_DIR / 'network.json'
    test_file = tmpdir / 'network.json'
    ref_objects = __read_multi_object_json(ref_file)
    test_objects = __read_multi_object_json(test_file)
    assert not deepdiff.DeepDiff(ref_objects, test_objects, ignore_order=True)


def test_with_cofactors(mocker, tmpdir):
    """Test the CLI with no cofactor file."""
    args = ['prog', str(REF_IN_TAR), str(tmpdir), '--cofactor-file', str(COF_FILE)]
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    test_file = tmpdir / 'network.json'
    test_objects = __read_multi_object_json(test_file)
    assert any(
        node['data']['cofactor']
        for node in test_objects['network']['elements']['nodes']
        if node['data']['type'] == 'chemical'
    )


def test_without_cofactors(mocker, tmpdir):
    """Test the CLI without a cofactor file."""

    # Case 1: with cofactor file set to None
    args = ['prog', str(REF_IN_TAR), str(tmpdir), '--cofactor-file', 'None']
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    test_file = tmpdir / 'network.json'
    test_objects = __read_multi_object_json(test_file)
    assert not any(
        node['data']['cofactor']
        for node in test_objects['network']['elements']['nodes']
        if node['data']['type'] == 'chemical'
    )

    # Case 2: with no cofactor detection switch
    args = ['prog', str(REF_IN_TAR), str(tmpdir), '--no-cofactor-detection']
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    test_file = tmpdir / 'network.json'
    test_objects = __read_multi_object_json(test_file)
    assert not any(
        node['data']['cofactor']
        for node in test_objects['network']['elements']['nodes']
        if node['data']['type'] == 'chemical'
    )
