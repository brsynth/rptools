# rpScore

Computes a global score for a heterologous pathway. The score is calculated from a learning process based on reaction rules score, flux balance analysis and thermodynamics metrics, and the number of reactions in the pathway.


## Input

Positional argument:
* **infile**: (string) Path to pathway (rpSBML) file
* **outfile**: (string, required if one single input file is passed) path to store the scored pathway

Optional arguments:
* **--no_of_rxns_thres**: (int, default=10) number of reactions above which a pathway is not scored (too long)
* **--data_train_file**: (string) path to the trained data
* **--log**: (string, default=error) Set the log level, choices are 'debug', 'info', 'warning', 'error', 'critical'
* **--version**: Display program's version and exit


## Install
Please see `rptool` documentation.

## Run

<!-- ### rpScore process -->
**From Python code**
```python
from rptools.rpscore import predict_score
from rptools.rplibs import rpPathway

pathway = rpPathway(
    infile='tests/rpscore/data/pathway.xml'
)

global_score = predict_score(pathway)
```
**From CLI**
```sh
python -m rptools.rpscore <input_rpsbml> <output_rpsbml>
```

## Tests
Test can be run with the following commands:

### Natively
```bash
cd tests
pytest -v
```

## CI/CD
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).

## Authors

* **Jean-Loup Faulon**
* **Joan Hérisson**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thomas Duigou

<!-- ### How to cite rpScore? -->