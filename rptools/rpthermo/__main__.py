#!/usr/bin/env python

from logging import Logger, getLogger
from os import makedirs as os_makedirs
from os import path as os_path

from colored import attr, fg

from rptools import build_args_parser
from rptools.rplibs import rpPathway
from rptools.rpthermo import runThermo
from rptools.rpthermo.args import add_arguments


def _cli():
    parser = build_args_parser(
        prog="rpthermo",
        description="Calculate score by processing thermodynamics",
        m_add_args=add_arguments,
    )
    args = parser.parse_args()

    from rptools.__main__ import init

    logger = init(parser, args)

    msg = "Parameters\n----------\n"
    for param in ["pH", "ionic_strength", "pMg"]:
        value = getattr(args, param)
        msg += f"- {param}: {value}\n"
    logger.info(
        "{color}{msg}{rst}".format(color=fg("light_cyan"), msg=msg, rst=attr("reset"))
    )

    ## READ PATHWAY FROM FILE
    pathway = rpPathway(infile=args.infile, logger=logger)

    # RUN THERMO
    results = runThermo(
        pathway=pathway,
        ph=args.pH,
        ionic_strength=args.ionic_strength,
        pMg=args.pMg,
        logger=logger,
    )

    # Print results
    print_results(pathway, results, logger)
    # Write pathway into file
    # Create the output directory if not exists
    os_makedirs(os_path.dirname(args.outfile), exist_ok=True)
    pathway.to_rpSBML().write_to_file(args.outfile)
    logger.info(
        "{color}{typo}Written into file: {file}{rst}".format(
            color=fg("white"), typo=attr("bold"), rst=attr("reset"), file=args.outfile
        )
    )


def print_results(
    pathway: rpPathway, results: dict, logger: Logger | None = None
) -> None:
    logger = logger or getLogger(__name__)
    logger.info(
        "{color}{typo}Results {net_rxn}{rst}".format(
            net_rxn=results["optimized_net_reaction"],
            color=fg("white"),
            typo=attr("bold"),
            rst=attr("reset"),
        )
    )
    print_thermo_results(results["net_reaction"], logger)
    for rxn_id in pathway.get_reactions_ids():
        logger.info(
            "{color}{typo}Results {rxn}{rst}".format(
                rxn=pathway.get_reaction(rxn_id),
                color=fg("white"),
                typo=attr("bold"),
                rst=attr("reset"),
            )
        )
        print_thermo_results(results["reactions"][rxn_id], logger)


def print_thermo_results(results: dict, logger: Logger | None = None) -> None:
    logger = logger or getLogger(__name__)
    for key, value in results.items():
        logger.info(
            "   {color}{typo}|- {key} = {value} {rst}+/- {error} {units}".format(
                color=fg("white"),
                typo=attr("bold"),
                rst=attr("reset"),
                key=key,
                value=value["value"],
                error=value["error"],
                units=value["units"],
            )
        )


if __name__ == "__main__":
    _cli()
