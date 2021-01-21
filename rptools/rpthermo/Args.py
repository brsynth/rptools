from argparse import ArgumentParser


def build_args_parser():
    parser = ArgumentParser(prog='rpthermo', description='Python wrapper to add cofactors to generate rpSBML collection')
    parser = _add_arguments(parser)

    return parser


def _add_arguments(parser):

    # positional arguments
    parser.add_argument('input',            type=str)
    parser.add_argument('output',           type=str)

    # optional arguments
    parser.add_argument('--pathway_id',     type=str,   default='rp_pathway')
    parser.add_argument('--ph',             type=float, default=7.5)
    parser.add_argument('--ionic_strength', type=float, default=200.0)
    parser.add_argument('--pMg',            type=float, default=10.0)
    parser.add_argument('--temp_k',         type=float, default=298.15)
    params = parser.parse_args()

    return parser
