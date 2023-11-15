from os import path as os_path, makedirs as os_makedirs
from sys import exit as sys_exit
from tempfile import NamedTemporaryFile
from errno import EEXIST as errno_EEXIST
from rptools import build_args_parser
from rptools.__main__ import init
from .Args import add_arguments
from .rpfba import (
    preprocess,
    runFBA,
    build_results,
    write_results_to_pathway
)


def _make_dir(filename):
    dirname = os_path.dirname(filename)
    if dirname != "" and not os_path.exists(dirname):
        try:
            os_makedirs(dirname)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno_EEXIST:
                raise


def entry_point():
    parser = build_args_parser(
        prog="rpfba",
        description="Process to Flux Balance Analysis",
        m_add_args=add_arguments,
    )
    args = parser.parse_args()

    logger = init(parser, args)

    # PREPROCESSING
    (
        merged_model,
        pathway,
        ids
    ) = preprocess(
        args=args,
        logger=logger
    )

    # FBA
    results = runFBA(
        model=merged_model,
        compartment_id=ids['comp_id'],
        biomass_rxn_id=ids['biomass_rxn_id'],
        objective_rxn_id=ids['obj_rxn_id'],
        sim_type=args.sim,
        fraction_coeff=args.fraction_of,
        logger=logger,
    )
    # with NamedTemporaryFile() as tmpfile:
    #     merged_model.write_to_file(tmpfile.name)
    #     results = runFBA(
    #         model_file=tmpfile.name,
    #         compartment_id=ids['comp_id'],
    #         biomass_rxn_id=ids['biomass_rxn_id'],
    #         objective_rxn_id=ids['obj_rxn_id'],
    #         sim_type=args.sim,
    #         fraction_coeff=args.fraction_of,
    #         hidden_species=hidden_species,
    #         logger=logger,
    #     )

    # RESULTS
    hidden_species = merged_model.get_isolated_species()
    results = build_results(
        results=results,
        pathway=pathway,
        compartment_id=ids['comp_id'],
        hidden_species=hidden_species,
        logger=logger,
    )

    # # Remove the Cobra standard ('compound@compartment') from all compounds
    # pathway.uncobraize()

    # results = uncobraize_results(results, cobra_suffix(compartment_id))

    # Write results into the pathway
    write_results_to_pathway(pathway, results, logger)

    if not results:
        logger.info("No results written. Exiting...")
    else:
        logger.info("Writing into file...")
        _make_dir(args.outfile)
        pathway.write_to_file(args.outfile)
        logger.info("   |--> written in " + args.outfile)

    return 0

if __name__ == "__main__":
    sys_exit(entry_point())
