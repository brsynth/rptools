#!/usr/bin/env python

from rptools.rpthermo import runThermo
from rptools.rpthermo.Args import add_arguments
from rptools import build_args_parser


def _cli():
    parser = build_args_parser(
        prog = 'rpthermo',
        description = 'Calculate score by processing thermodynamics',
        m_add_args = add_arguments
    )
    args   = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    runThermo(args.rpsbml_infile, args.rpsbml_outfile,
              args.pathway_id,
              args.ph, args.ionic_strength, args.pMg, args.temp_k,
              logger=logger)


if __name__ == '__main__':
    _cli()
