# rpCompletion

Completes mono-component reactions output by RetroPath2.0 with the appropriate cofactors. Creates sub-paths when multiple reaction rules are associated with a single reaction. Input is a single pathways file produced by RP2Paths. It stands on rrCache which store pre-computed data.

All metabolic pathways will be built as the following:

- In each "master" pathway,
  - each chemical transformation could be produced by multiple reaction rules, and
    - each reaction rule could have been produced by multiple template chemical reactions.

Thus, each different template reaction for each different reaction rule for each different chemical transformation provides one single possible pathway. The algorithm explores the combinatorics of all possible pathways and for each "master pathway" (the one from chemical transformations), keeps only top (defined by `max_subpaths_filter` CLI option, default: 10) ones.

## Input

Required:
* **rp2_metnet**: (string) Path to the metabolic network file built by RetroPath2.0
* **sink**: (string) Path to the rpextractsink file containing infos on molecules in the sink
* **rp2paths_compounds**: (string) Path to the rp2paths compounds file
* **rp2paths_pathways**: (string) Path to the rp2paths pathways file
* **outdir**: (string) Path to the folder where result files are written

Advanced options:
* **--upper_flux_bound**: (integer, default=10000) Upper flux bound value for all new reactions created
* **--lower_flux_bound**: (integer, default=-10000) Lower flux bound value for all new reactions created
* **--max_subpaths_filter**: (integer, default=10) Number of subpaths per master pathway



<!-- ## Memory management

### File mode
This is the default mode. All cache data are stored into files on disk and loaded in memory each time the tool is used. In this mode, fingerprint in memory is equal to the size of cache files loaded in memory multiplied by the number of processes which are running at the same time. Option can be specified by `--store-mode file`.

### DB mode
In order to save memory space, cache data can be loaded once in a database (redis) so that the memory space taken is equal to one instance of the cache, whatever the number of processes whic are running. Option can be specified by `--store-mode <db_host>`, where `db_host` is the hostname on which redis server is running. -->


## Install
Please see `rptool` documentation.

## Run

### rpCompletion process
**From Python code**
```python
from rptools.rpcompletion import rp_completion

pathways = rp_completion(
    rp2_metnet_filename,
    sink_filename,
    rp2paths_compounds_filename,
    rp2paths_pathways_filename,
)
```
**From CLI**
```sh
python -m rptools.rpcompletion \
  rp2_metnet.csv \
  sink.csv \
  rp2paths_compounds.csv \
  rp2paths_pathways.csv \
  <outdir>
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
* Melchior du Lac

## Licence
rpCompletion is released under the MIT licence. See the LICENCE file for details.
