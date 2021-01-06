# rpTools

[![Anaconda-Server Badge](https://anaconda.org/brsynth/rptools/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rptools)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/rptools/badges/version.svg)](https://anaconda.org/brsynth/rptools)
![Test](https://github.com/brsynth/rpTools/workflows/Test/badge.svg)

rpTools are dedicated to work around rpSBML data structure. Tools are the following:

* rpCompletion: [README.md](rptools/rpcompletion/README.md)
* rpFBA: [README.md](rptools/rpfba/README.md)
* rpExtractSink: [README.md](rptools/rpextractsink/README.md)
* rpLibs: [README.md](rptools/rplibs/README.md)

## Install
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda rptools
```

## Run
Please see tool documentation.

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
* **Thomas Duigou**

## Licence
rpTools is released under the MIT licence. See the LICENCE file for details.
