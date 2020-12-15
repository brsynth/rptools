# rpCompletion

[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpcompletion/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rpcompletion)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/rpcompletion/badges/version.svg)](https://anaconda.org/brsynth/rpcompletion)
![Test suite](https://github.com/brsynth/rpCompletion/workflows/Test%20suite/badge.svg)

Completes mono-component reactions output by RetroPath2.0 with the appropriate cofactors. Creates sub-paths when multiple reaction rules are associated with a single reaction. Input is a single pathways file produced by RP2Paths. It stands on rpCache which store pre-computed data.

## Input

Required:
* **rp2_pathways**: (string) Path to the RetroPath2.0 pathways file
* **rp2paths_compounds**: (string) Path to the rp2paths compounds file
* **rp2paths_pathways**: (string) Path to the rp2paths pathways file
* **outdir**: (string) Path to the folder where result files are written

Advanced options:
* **-upper_flux_bound**: (integer, default=9999) Upper flux bound value
* **-lower_flux_bound**: (integer, default=0) Lower flux bound value
* **-maxSubPaths_filter**: (integer, default=10) Number of subpaths per path
* **-pathway_id**: (string, default=rp_pathway) ID of the heterologous pathway
* **-compartment_id**: (string, default=MNXC3 (i.e. cytoplasm)) Heterologous pathway compartment ID
* **-species_group_id**: (string, default=central_species) ID of the central species, i.e. not cofactors, in the heterologous reactions
* **--store-mode, -sm**: (optional, string, default: file) Store mode. If 'file', rpCache is supposed to be stored in files. Else, the rpCache is supposed to be stored in a redis database which the name is the value of this input field. The redis server is considered to be up and running.



## Memory management

### File mode
This is the default mode. All cache data are stored into files on disk and loaded in memory each time the tool is used. In this mode, fingerprint in memory is equal to the size of cache files loaded in memory multiplied by the number of processes which are running at the same time. Option can be specified by `--store-mode file`.

### DB mode
In order to save memory space, cache data can be loaded once in a database (redis) so that the memory space taken is equal to one instance of the cache, whatever the number of processes whic are running. Option can be specified by `--store-mode <db_host>`, where `db_host` is the hostname on which redis server is running.


## Install
rpCompletion requires [RDKit](https://www.RDKit.org) which is not available through pip. It can be installed through Conda:
```sh
[sudo] conda install -c rdkit rdkit
```
### From pip
```sh
[sudo] python -m pip install rpcompletion
```
### From Conda
```sh
[sudo] conda install -c brsynth rpcompletion
```

## Run

### rpCompletion process
**From Python code**
```python
from rpcompletion import rpCompletion, build_args_parser

parser = build_args_parser()
args  = parser.parse_args()

rpcompletion = rpCompletion(db=args.store_mode)
rpcompletion.rp2ToSBML(args.rp2_pathways,
                       args.rp2paths_compounds,
                       args.rp2paths_pathways,
                       args.outdir)
```
**From CLI**
```sh
python -m rpcompletion \
  rp2_pathways.csv \
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

* **Melchior du Lac**
* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou


## Licence
rpCompletion is released under the MIT licence. See the LICENCE file for details.
