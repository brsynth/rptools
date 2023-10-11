# rpExtractSink

RetroPath2 sink generator. From a given SBML file, the tool will extract all the sink molecules and generate a csv file with the InChI structures. For each molecule, the InChI will be get from:

1. The local cache (rrCache), if available
2. The MetaNetX database from MIRIAM URLs (MetaNetX first)

## Input

Required:
* **input_sbml**: (string) Path to the input SBML file

Optional:
* **--remove-dead-end**: (boolean, default: True) Perform FVA evaluation to remove dead end metabolites
* **--compartment-id**: (string, default: 'c') Specify the compartment from which to extract the sink molecules. The default are for MetaNetX files
* **--standalone**: (boolean, default: False) If True, do not retrieve InChI from Internet
* **--cache-dir**: (string, default: None) Path to the cache directory

## Output

* **output_sbml**: (string) Path to the output csv file


## Install
Please see `rptool` documentation.

## Use

### Function call from Python code
```python
from rr_cache import rrCache
from rptools.rpextractsink import genSink

cache = rrCache(
    attrs=['cid_strc'],
    cache_dir=args.cache_dir,
    logger=logger
)
sink = genSink(
    cache,
    args.input_sbml
)
```


### Run from CLI
```sh
python -m rptools.rpextractsink --help
```

## Tests
Test can be run with the following commands:

### Natively
```bash
python -m pytest tests/rpextractsink
```

## CI/CD
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).

## Authors

* **Joan Hérisson**
* Thomas Duigou, Melchior du Lac

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
