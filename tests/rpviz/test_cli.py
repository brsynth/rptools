import pytest
from pathlib import Path
from rptools.rpviz.__main__ import (
    __build_arg_parser,
    __run,
)


REF_IN_DIR = Path(__file__).resolve().parent / 'inputs' / 'as_dir'
REF_IN_TAR = Path(__file__).resolve().parent / 'inputs' / 'as_tar.tgz'
REF_OUT_DIR = Path(__file__).resolve().parent / 'outputs'


def __read_and_sort(
    path: str
    ) -> list:
    with open(path) as fh:
        lines = [_ for _ in fh.readline()]
    return lines.sort()


def test_build_arg_parser(mocker):
    # No args
    args = ['prog']
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    with pytest.raises(SystemExit):
        args = parser.parse_args()


def test_dir_input(mocker, tmpdir):
    out_dir = tmpdir / 'odir'
    args = ['prog', str(REF_IN_DIR), str(out_dir)]
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    ref_file = REF_OUT_DIR / 'network.json'
    test_file = out_dir / 'network.json'
    assert __read_and_sort(test_file) == __read_and_sort(ref_file)


def test_tar_input(mocker, tmpdir):
    out_dir = tmpdir / 'odir'
    args = ['prog', str(REF_IN_TAR), str(out_dir)]
    mocker.patch('sys.argv', args)
    parser = __build_arg_parser()
    args = parser.parse_args()
    __run(args)
    ref_file = REF_OUT_DIR / 'network.json'
    test_file = out_dir / 'network.json'
    assert __read_and_sort(test_file) == __read_and_sort(ref_file)