# rpRank

Rank a list of pathways according their global score.


## Input

Positional arguments:
* **pathways**: (list of string) Paths to pathway (rpSBML) files

## Install
Please see `rptool` documentation.

## Run

<!-- ### rpScore process -->
**From Python code**
```python
from rptools.rpscore import predict_score
from rptools.rplibs import rpPathway

pathway = rpPathway.from_rpSBML(
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