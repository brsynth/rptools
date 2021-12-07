# rpExtractSink


RetroPath2 sink generator

## Input

Required:
* **input_sbml**: (string) Path to the input SBML file

Optional:
* **--remove_dead_end**: (boolean, default: True) Perform FVA evaluation to remove dead end metabolites
* **--compartment_id**: (string, default: MNXC3) Specify the compartment from which to extract the sink molecules. The default are for MetaNetX files

## Output

* **output_sbml**: (string) Path to the output csv file


## Install
Please see `rptool` documentation.

## Use

### Function call from Python code
```python
from rptools.rpextractsink import rpextractsink

sink = rpExtractSink(input_sbml, output_sink)
sink.genSink()
```

If parameters from CLI have to be parsed, the function `build_args_parser` is available:
```python
from rptools.pextractsink import build_args_parser

parser = buildparser()
params = parser.parse_args()
```

### Run from CLI
```sh
python -m rptools.rpextractsink \
    <input_sbml> \
    <output_sink> \
    [--compartment_id COMPARTMENT_ID] \
    [--remove_dead_end REMOVE_DEAD_END]
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

* **Melchior du Lac**
* Thomas Duigou, Joan Hérisson

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
