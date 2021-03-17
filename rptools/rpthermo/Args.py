from argparse import ArgumentParser


def add_arguments(parser):

    # positional arguments
    parser.add_argument('rpsbml_infile',    type=str, help='input SBML file')
    parser.add_argument('rpsbml_outfile',   type=str, help='output SBML file')

    # optional arguments
    parser.add_argument('--pathway_id'     , type=str   , default='rp_pathway')
    parser.add_argument('--ph'             , type=float , default=7.5         )
    parser.add_argument('--ionic_strength' , type=float , default=200.0       )
    parser.add_argument('--pMg'            , type=float , default=10.0        )
    parser.add_argument('--temp_k'         , type=float , default=298.15      )

    return parser
