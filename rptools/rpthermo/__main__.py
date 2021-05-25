#!/usr/bin/env python

from json import (
    dump as json_dump,
    load as json_load
)
from logging import (
    Logger,
    getLogger
)
from typing import (
    Dict,
    List,
    Tuple
)
from colored import fg, bg, attr
from rptools.rpthermo import thermo
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

    ## READ PATHWAY FROM FILE
    with open(args.infile, 'r') as fp:
        pathway = json_load(fp)

    ## RUN THERMO
    pathway = thermo(
        pathway,
        args.pathway_id,
        args.ph,
        args.ionic_strength,
        args.pMg,
        args.temp_k,
        logger=logger
    )

    ## WRITE PATHWAY TO FILE
    with open(args.outfile, 'w') as fp:
        json_dump(pathway, fp, indent=4)

    ## PRINT OUT RESULTS
    print_results(
        pathway['measures']['thermo'],
        logger
    )


def print_results(
    results: Dict,
    logger: Logger=getLogger(__name__)
) -> None:
    logger.info(
        "{color}{typo}Results{rst}".format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    logger.info(
        "   {color}{typo}|- ΔG'° = {value} {rst}+/- {error} {units}".format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset'),
            value=results['dG0_prime']['value'],
            error=results['dG0_prime']['error'],
            units=results['dG0_prime']['units']
        )
    )
    logger.info(
        "   {color}{typo}|- ΔG'm = {value} {rst}+/- {error} {units}".format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset'),
            value=results['dGm_prime']['value'],
            error=results['dGm_prime']['error'],
            units=results['dGm_prime']['units']
        )
    )
    logger.info(
        "   {color}{typo}|- ΔG' = {value} {rst}+/- {error} {units}".format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset'),
            value=results['dG_prime']['value'],
            error=results['dG_prime']['error'],
            units=results['dG_prime']['units']
        )
    )


if __name__ == '__main__':
    _cli()
