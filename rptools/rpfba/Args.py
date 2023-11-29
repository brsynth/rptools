from argparse import ArgumentParser
from typing import (
    List,
)
from rptools._version import __version__


DEFAULT_ARGS = {
    "pathway_file": "",
    "model_file": "",
    "compartment_id": "c",
    "outfile": "",
    "objective_rxn_id": "rxn_target",
    "biomass_rxn_id": "biomass",
    "sim": "fraction",
    "fraction_coeff": 0.75,
    "merge": "",
    "with_orphan_species": False,
}

def add_arguments(parser: ArgumentParser):
    parser.add_argument(
        "pathway_file", type=str, help="SBML file that contains an heterologous pathway"
    )
    parser.add_argument("model_file", type=str, help="GEM model file (SBML)")
    parser.add_argument(
        "compartment_id",
        type=str,
        help="model compartment id to consider (e.g. 'c' or 'MNXC3')",
    )
    parser.add_argument("outfile", type=str, help="output file")
    parser.add_argument(
        "--objective_rxn_id",
        type=str,
        default=DEFAULT_ARGS["objective_rxn_id"],
        help="reaction ID to optimise (default: rxn_target)",
    )
    parser.add_argument(
        "--biomass_rxn_id",
        type=str,
        default=DEFAULT_ARGS["biomass_rxn_id"],
        help="biomass reaction ID (default: biomass). Note: Only for 'fraction' simulation",
    )
    parser.add_argument(
        "--sim",
        type=str,
        choices=["fba", "pfba", "fraction"],
        default=DEFAULT_ARGS["sim"],
        help="type of simulation to use (default: fraction)",
    )
    parser.add_argument(
        "--fraction_of",
        type=float,
        default=DEFAULT_ARGS["fraction_coeff"],
        help="fraction of the optimum (default: 0.75). Note: this value is ignored is 'fba' is used",
    )
    parser.add_argument(
        "--merge",
        type=str,
        default=DEFAULT_ARGS["merge"],
        help="output the full merged model in addition of heterologous pathway only (default: False)",
    )
    parser.add_argument(
        "--with_orphan_species",
        action="store_true",
        default=DEFAULT_ARGS["with_orphan_species"],
        help="Take metabolites that are only consumed (default: False)",
    )

    return parser
