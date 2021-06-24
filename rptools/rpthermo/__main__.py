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
from brs_utils import (
    print_OK_adv as print_OK,
    print_title_adv as print_title
)
from rptools.rpthermo import runThermo
from rptools.rpthermo.Args import add_arguments
from rptools import build_args_parser
from rptools.rplibs import (
    rpSBML,
    rpPathway
)


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
    pathway = rpSBML(
        inFile=args.infile,
        logger=logger
    ).to_Pathway()

    # RUN THERMO
    results = runThermo(
        pathway,
        args.ph,
        args.ionic_strength,
        args.pMg,
        args.temp_k,
        logger=logger
    )

    # Print results
    print_results(pathway, results, logger)
    # Write pathway into file
    rpSBML.from_Pathway(pathway).write_to_file(args.outfile)
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
            net_rxn=pathway.net_reaction(),
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    print_thermo_results(
        results['net_reaction'],
        # pathway.get_infos()['thermo'],
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
            # pathway.get_reaction(rxn_id).get_info('thermo'),
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
