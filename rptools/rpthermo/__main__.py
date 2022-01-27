#!/usr/bin/env python

from logging import (
    Logger,
    getLogger
)
from typing import Dict
from colored import fg, bg, attr
from brs_utils import (
    print_OK_adv as print_OK,
    print_title_adv as print_title
)
from rptools.rpthermo import runThermo
from rptools.rpthermo.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import rpPathway


def _cli():
    parser = build_args_parser(
        prog = 'rpthermo',
        description = 'Calculate score by processing thermodynamics',
        m_add_args = add_arguments
    )
    args   = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)

    msg = f'Parameters\n----------\n'
    for param in ['pH', 'ionic_strength', 'pMg']:
        value = getattr(args, param)
        msg += f'- {param}: {value}\n'
    logger.info(
        '{color}{msg}{rst}'.format(
            color=fg('light_cyan'),
            msg=msg,
            rst=attr('reset')
        )
    )

    ## READ PATHWAY FROM FILE
    pathway = rpPathway.from_rpSBML(
      infile=args.infile,
      logger=logger
    )

    # RUN THERMO
    results = runThermo(
        pathway=pathway,
        ph=args.pH,
        ionic_strength=args.ionic_strength,
        pMg=args.pMg,
        logger=logger
    )

    # Print results
    print_results(pathway, results, logger)
    # Write pathway into file
    pathway.to_rpSBML().write_to_file(args.outfile)
    logger.info(
        "{color}{typo}Written into file: {file}{rst}".format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset'),
            file=args.outfile
        )
    )


def print_results(
    pathway: rpPathway,
    results: Dict,
    logger: Logger=getLogger(__name__)
) -> None:
    logger.info(
        "{color}{typo}Results {net_rxn}{rst}".format(
            net_rxn=results['optimized_net_reaction'],
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    print_thermo_results(
        results['net_reaction'],
        logger
    )
    for rxn_id in pathway.get_reactions_ids():
        logger.info(
            "{color}{typo}Results {rxn}{rst}".format(
                rxn=pathway.get_reaction(rxn_id),
                color=fg('white'),
                typo=attr('bold'),
                rst=attr('reset')
            )
        )
        print_thermo_results(
            results['reactions'][rxn_id],
            logger
        )

def print_thermo_results(
    results: Dict,
    logger: Logger=getLogger(__name__)
) -> None:
    for key, value in results.items():
        logger.info(
            "   {color}{typo}|- {key} = {value} {rst}+/- {error} {units}".format(
                color=fg('white'),
                typo=attr('bold'),
                rst=attr('reset'),
                key=key,
                value=results[key]['value'],
                error=results[key]['error'],
                units=results[key]['units']
            )
        )


if __name__ == '__main__':
    _cli()
