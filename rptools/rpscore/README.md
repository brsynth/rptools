# rpScore

Computes a global score for a pathway from the following:


<p><img src="https://render.githubusercontent.com/render/math?math=GS = \frac{1}{N} \sum_{i=1}^{N} Scores[i] \cdot Weights[i]"></p>

<p>with <img src="https://render.githubusercontent.com/render/math?math=N = \overline{Scores}">, and <img src="https://render.githubusercontent.com/render/math?math=Scores = [RuleScore, FBAScore, ThermoScore, NbReactions]"></p>


The global score (*GS*) is the weighted arithmetic mean of the following scores:

* *RuleScore*: arithmetic mean of each reaction's rule score in the pathway

* *FBAScore*: score of the flux balance analysis of the pathway

* *ThermoScore*: score of the thermodynamics of the pathway

* *NbReactions*: number of reactions in the pathway

Weights are given by the user and have the following default values:

* *weight_NbReactions* = 0.10002239003499142
* *weight_Rule* = 0.13346271414277305
* *weight_FBA* = 0.6348436269211155
* *weight_Thermo* = 0.13167126890112002


## Input

Positional argument:
* **pathway_file**: Path to a pathway (rpSBML) file

Optional arguments:
* **--outfile**: path to the file where the input file augmented with the global score will be written
* **--fba_ceil**: (float, default=3.0) FBA ceiling
* **--fba_floor**: (float, default=0.0) FBA floor
* **--thermo_ceil**: (float, default=8901.2) Thermodynamics ceiling
* **--thermo_floor**: (float, default=-7570.2) Thermodynamics floor
* **--weight_rp_steps**: (float, default=0.10002239003499142) Number of steps weight
* **--max_rp_steps**: (integer, default=15) Maximal number of steps (as run in RetroPath2.0)
* **--weight_rule_score**: (float, default=0.13346271414277305) Reation rule score weight
* **--weight_fba**: (float, default=0.6348436269211155) FBA weight
* **--weight_thermo**: (float, default=0.13167126890112002) Pathway thermodynamics weight
* **--pathway_id**: (string, default=rp_pathway) Name of the heterologous pathway
* **--objective_id**: (string, default=obj_RP1_sink__restricted_biomass) Name of the heterologous pathway objective function
* **--thermo_id**: (string, default=dfG_prime_m) Name of the thermodynamics
* **--log**: (string, default=error) Set the log level, choices are 'debug', 'info', 'warning', 'error', 'critical'
* **--silent**: (default=False) Runs the program silently
* **--version**: Display program's version and exit

## Output

Depending of which option(s) is(are) set:
**No option**: print the global score of the pathway
**--outfile**: write a new rpSBMl file with global score added in it


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