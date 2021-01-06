# rpFBA

Perform FBA on a single or collection of SBML files containing heterologous pathways, as tar.xz archives. The package performs the following steps: 1) it merges a user defined GEM SBML model to a given heterologous pathway. 2) it performs FBA using the [cobrapy](https://opencobra.github.io/cobrapy/) package using a user defined mathod that include, FBA, parsimonious FBA or fraction of optimum of another reaction. For the first two, the user must know the reaction name that the model will optimise to, while the latter the use must provide the target reaction but also another reaction that will be restricted. The first step involves performing FBA using the "source" reaction as the objective. Then the flux of that reaction has its upper and lower bounds set to the same value, determined as a fraction of its FBA flux value. Thereafter the objective is set to the initial target reaction and FBA is performed once again. The tool uses the [FBC](https://co.mbine.org/specifications/sbml.level-3.version-1.fbc.version-2.release-1) package to manage the objective and flux bounds.

## Input

Required:
* **input_sbml**: (string) Path to the input file
* **gem_sbml**: (string) Path to the GEM SBML model

Advanced options:
* **--pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway
* **--compartment_id**: (string, default=MNXC3 (i.e. cytoplasm)) ID of the compartment ID that contains the heterologous pathway
* **--sim_type**: (string, default=fraction) Valid options include: fraction, fba, pfba. The type of constraint based modelling method
* **--source_reaction**: (string, default=biomass) Name of the source reaction that will be restricted in the "fraction" simulation type. This parameter is ignored for "fba" and "pfba"
* **--target_reaction**: (string, default=RP1_sink) Heterologous pathway flux sink reaction. This parameters is required in all simulation type
* **--source_coefficient**: (float, default=1.0) Objective coefficient for the source reaction. This parameter is ignored for "fba" and "pfba"
* **--target_coefficient**: (float, default=1.0) Objective coefficient for the target reaction.
* **--is_max**: (boolean, default=True) Maximise or minimise the objective function
* **--fraction_of**: (float, default=0.75) Portion of the maximal flux used to set the maximal and minimal bounds for the source reaction of the "fraction" simulation type
* **--dont_merge**: (boolean, default=True) Return the merged GEM+heterologous pathway SBML or only the heterologous pathway SBML files
* **--log**: (string, default=error) Set the log level, choices are 'debug', 'info', 'warning', 'error', 'critical'

## Output

* **output**: (string) Path to the output file


## Install
rpFBA is part of rpTools suite:
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda rptools
```

## Run

### rpFBA process
**From Python code**
```python
from rptools.rpfba import runFBA, build_args_parser

parser = build_args_parser()
args  = parser.parse_args()

logger = logging.getLogger(__name__)
logger.setLevel(getattr(logging, args.log.upper()))
logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

result = runFBA(args.input_sbml, args.gem_sbml, args.outfile,
                args.sim_type,
                args.source_reaction, args.target_reaction,
                args.source_coefficient, args.target_coefficient,
                args.is_max,
                args.fraction_of,
                args.dont_merge,
                args.pathway_id,
                args.objective_id,
                args.compartment_id,
                args.species_group_id,
                args.sink_species_group_id,
                logger)
```
**From CLI**
```sh
python -m rptools.rpfba <input_sbml> <gem_sbml> <outfile>
```

## Tests
Test can be run with the following commands:

### Natively
```bash
cd tests
pytest -v
```

# CI/CD
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).


## Authors

* **Melchior du Lac**
* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou


## Licence
rpFBA is released under the MIT licence. See the LICENCE file for details.
