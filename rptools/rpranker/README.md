# rpRanker

Rank pathways according to their global score (computed by [rpScore](https://github.com/brsynth/rptools/tree/master/rptools/rpscore) tool).


## Input

Positional arguments:
* **pathways**: Pathway (rpSBML) filenames

Optional arguments:
* **--top**: (default=10) Numbers of pathways to be selected
* **--rank_outfile**: Path to store the ranking result
* **--data_outfile**: Path to store selected pathways as a .tar.gz archive
* **--data_outdir**: Path to store selected pathways within a folder
* **--rename**: (default=False) Rename files when --data_outdir or --data_outfile is set
* **--log**: (string, default=error) Set the log level, choices are 'debug', 'info', 'warning', 'error', 'critical'
* **--silent**: (default=False) Runs the program silently
* **--version**: Display program's version and exit

## Output

Depending of which option(s) is(are) set:
**No option**: print the list of sorted pathways
**--rank_outfile**: store the list of sorted pathways into a file
**--data_outfile**: store ranked pathways into a compressed file
**--data_outdir**: write to disk the ranked pathways into a folder


## Install
rpRanker is part of rpTools suite:
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda rptools
```

## Run

### rpRanker process
**From Python code**
```python
from rptools.rpranker import rank, build_args_parser

parser = build_args_parser()
args  = parser.parse_args()

from rptools.__main__ import init
logger = init(parser, args)

ranked_pathways = rank(
    pathways = args.pathways,
    logger = logger
)
```
**From CLI**
```sh
python -m rptools.ranker pathway_1 [pathway_2 ...]
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

## Acknowledgments

* Thomas Duigou


## Licence
rpRanker is released under the MIT licence. See the LICENCE file for details.
