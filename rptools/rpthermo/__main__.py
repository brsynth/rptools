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
from rptools.rpthermo import thermo
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
    results = thermo(
        pathway,
        args.ph,
        args.ionic_strength,
        args.pMg,
        args.temp_k,
        logger=logger
    )

    # results = {'net_reaction': {'dG0_prime': {'value': -789.8933819367826, 'error': 7.7248778091405255, 'units': 'kilojoule / mole'}, 'dGm_prime': {'value': -789.8933819367826, 'error': 7.7248778091405255, 'units': 'kilojoule / mole'}, 'dG_prime': {'value': -789.8933819367826, 'error': 7.7248778091405255, 'units': 'kilojoule / mole'}}, 'reactions': {'rxn_2': {'dG0_prime': {'value': -324.1942486194258, 'error': 5.946918851751192, 'units': 'kilojoule / mole'}, 'dGm_prime': {'value': -307.0794110847048, 'error': 5.946918851751192, 'units': 'kilojoule / mole'}, 'dG_prime': {'value': -324.1942486194258, 'error': 5.946918851751192, 'units': 'kilojoule / mole'}}, 'rxn_1': {'dG0_prime': {'value': -465.6991333173571, 'error': 5.39929292564583, 'units': 'kilojoule / mole'}, 'dGm_prime': {'value': -482.8139708520781, 'error': 5.39929292564583, 'units': 'kilojoule / mole'}, 'dG_prime': {'value': -465.6991333173571, 'error': 5.39929292564583, 'units': 'kilojoule / mole'}}}, 'species': {'MNXM188': {'standard_dg_formation': {'value': -201.33117647305943, 'units': 'kilojoule / mole'}}, 'MNXM15': {'standard_dg_formation': {'value': -81.91776181003092, 'units': 'kilojoule / mole'}}, 'MNXM4': {'standard_dg_formation': {'value': 16.399999999865035, 'units': 'kilojoule / mole'}}, 'MNXM6': {'standard_dg_formation': {'value': -3070.5260288612576, 'units': 'kilojoule / mole'}}, 'CMPD_0000000003': {'standard_dg_formation': {'value': -290.7486822708706, 'units': 'kilojoule / mole'}}, 'TARGET_0000000001': {'standard_dg_formation': {'value': -508.17794736445586, 'units': 'kilojoule / mole'}}, 'MNXM13': {'standard_dg_formation': {'value': -386.0000000000019, 'units': 'kilojoule / mole'}}, 'MNXM5': {'standard_dg_formation': {'value': -3098.931269475794, 'units': 'kilojoule / mole'}}, 'MNXM1': {'standard_dg_formation': {'value': None, 'units': 'kilojoule / mole'}}}}

    # Print results
    print_results(pathway, results, logger)

    # Write results into pathway
    write_results(pathway, results, logger)

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


def write_results(
  pathway: rpPathway,
  results: Dict,
  logger: Logger = getLogger(__name__)
) -> None:
  # Write species results
  for spe_id, score in results['species'].items():
    for k, v in score.items():
      pathway.get_specie(spe_id).add_info(
        key='thermo_'+k,
        value=v
      )
  # Write reactions results
  for rxn_id, score in results['reactions'].items():
    for k, v in score.items():
      pathway.get_reaction(rxn_id).add_info(
        key='thermo_'+k,
        value=v
      )
  # Write pathway result
  for k, v in results['net_reaction'].items():
    pathway.add_info(
      key='thermo_'+k,
      value=v
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
