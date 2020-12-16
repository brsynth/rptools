# brs-libs

[![Anaconda-Server Badge](https://anaconda.org/brsynth/brs_libs/badges/latest_release_date.svg)](https://anaconda.org/brsynth/brs_libs)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/brs_libs/badges/version.svg)](https://anaconda.org/brsynth/brs_libs)
![Test suite)](https://github.com/brsynth/brs-libs/workflows/Test%20suite/badge.svg)

Libraries for rpTools:
* rpSBML
* rpCache
* inchikeyMIRIAM

## rpSBML
Defines SBML structure with additional fields relative to [RetroPath2](https://github.com/brsynth/RetroPath2-wrapper) objects.

<!-- ### Prerequisites
* Python 3 with the following modules:
    * python-libsbml
    * [RDKit](https://www.RDKit.org) -->

## rpCache

### Memory management

#### File mode
This is the default mode. All cache data are stored into files on disk and loaded in memory each time the tool is used. In this mode, fingerprint in memory is equal to the size of cache files loaded in memory multiplied by the number of processes which are running at the same time. Option can be specified by `--store-mode file`.

#### DB mode
In order to save memory space, cache data can be loaded once in a database (redis) so that the memory space taken is equal to one instance of the cache, whatever the number of processes whic are running. Option can be specified by `--store-mode <db_host>`, where `db_host` is the hostname on which redis server is running.


### Install
#### From Conda
```sh
[sudo] conda install -c brsynth -c conda-forge brs_libs
```

### Use

#### Load rpCache in memory
**Full cache into files**
```python
from brs_libs import rpCache

rpcache = rpCache(db='file')
print(rpcache.cid_src)
```

**Full cache into Redis DB**
For multiple instances of rpCache simultaneously, rpCache can be loaded into one single Redis database:
```python
from brs_libs import rpCache

rpcache = rpCache(db='localhost')
print(rpcache.cid_src)
```
`localhost` means that rpCache will look for a redis database locally. If there is not, it will start a brand new redis server. `localhost` could be replaced by any hostname that hosts the Redis database.

**A part of cache**
For less loading time and memory footprint, a part of the cache can be loaded:
```python
from brs_libs import rpCache

rpcache = rpCache(attrs='cid_strc')
print(rpcache.cid_src)
```

#### (Re-)generate the cache
**From Python code**
```python
from brs_libs import rpCache

rpCache.generate_cache(outdir)
```

**From CLI**
After having installed brs_libs Python module:
```sh
python -m brs_libs --gen_cache <folder>
```


### Test
Please follow instructions below ti run tests:
```
cd tests
pytest -v
```
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).

## inchikeyMIRIAM
Uses the rpCache to parse an SBML file to find all the chemical species, and try to recover the inchikey and add it to the MIRIAM annotation.



## Authors

* **Melchior du Lac**
* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou


## Licence
brs_libs is released under the MIT licence. See the LICENCE file for details.
