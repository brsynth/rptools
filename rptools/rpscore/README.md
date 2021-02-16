# rpGlobalScore

Tool that reads a collection of rpSBML files (in a tar.xz) and calculates a global score based on a series of weights given by the user

## Getting Started

This is a docker galaxy tools, and thus, the docker needs to be built locally where Galaxy is installed. 

## Input

Required information:
* **-input**: (string) Path to either tar.xz input collection of rpSBML files or a single rpSBML file.
* **-input_format**: (string) Format of the input

Advanced options:
* **-fba_ceil**: (float, default=3.0) FBA ceiling
* **-fba_floor**: (float, default=0.0) FBA floor
* **-thermo_ceil**: (float, default=8901.2) Thermodynamics ceiling
* **-thermo_floor**: (float, default=-7570.2) Thermodynamics floor
* **-weight_rp_steps**: (float, default=0.0) Number of steps weight
* **-max_rp_steps**: (integer, default=15) Maximal number of steps (as run in RetroPath2.0)
* **-weight_rule_score**: (float, default=0.0) Reation rule score weight
* **-weight_fba**: (float, default=0.699707) FBA weight
* **-weight_thermo**: (float, default=0.8334961) Pathway thermodynamics weight
* **-pathway_id**: (string, default=rp_pathway) Name of the heterologous pathway
* **-objective_id**: (string, default=obj_RP1_sink__restricted_biomass) Name of the heterologous pathway objective function
* **-thermo_id**: (string, default=dfG_prime_m) Name of the thermodynamics

## Output

* **output**: (string) Path to the output file

## Install
rpScore is part of rpTools suite:
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda rptools
```

## Run

### rpScore process
**From Python code**
```python
from rptools.score  import compute_globalscore
from rptools.rplibs import rpSBML

rpsbml = rpSBML(inFile = args.pathway_file)

rpsbml_dict = compute_globalscore(rpsbml = rpsbml)

global_score = rpsbml_dict['pathway']['brsynth']['global_score']
```
**From CLI**
```sh
python -m rptools.rpscore <input_sbml>
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

* **Joan Hérisson**
* **Melchior du Lac**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou

### How to cite rpGlobalScore?