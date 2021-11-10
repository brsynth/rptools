# rpTools

[![Anaconda-Server Badge](https://anaconda.org/brsynth/rptools/badges/latest_release_date.svg)](https://anaconda.org/brsynth/rptools)
[![Anaconda-Server Badge](https://anaconda.org/brsynth/rptools/badges/version.svg)](https://anaconda.org/brsynth/rptools)
![Tests](https://github.com/brsynth/rpTools/workflows/Tests/badge.svg)

rpTools are dedicated to work around rpSBML data structure. Tools are the following:

* rpCompletion: [README.md](rptools/rpcompletion/README.md)
* rpFBA: [README.md](rptools/rpfba/README.md)
* rpThermo: [README.md](rptools/rpthermo/README.md)
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

## CI/CD
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](./ci/README.md)).

## For developers

### Development installation

After a git clone:

```sh
cd <repository>
conda env create -f environment.yaml -n <dev_env>
conda develop -n <dev_env> .
```

Warning: if you do not specify an environment name with `-n <dev_env>`, 
then 'rptools-dev' will be used.

Test your installation with:

```sh
conda activate <dev_env>
python -m rptools
python -m rptools.rpcompletion -h
```

To uninstall:

```sh
conda deactivate
conda env remove -n <dev_env>
```

### Development installation (alternative using the ci toolkit)

After a git clone:
```sh
git clone https://github.com/breakthewall/cicd-toolkit.git
cd cicd-toolkit
make test
cd ..
conda develop -n rptools_test .
```

This will create a `rptools_test` environnement that can be activated to further develop / debug rptools.
```sh
conda activate rptools_test
python -m rptools
python -m rptools.rpcompletion -h
```

## Authors

* **Melchior du Lac**
* **Joan Hérisson**
* **Thomas Duigou**

## Licence
rpTools is released under the MIT licence. See the [LICENCE file](./LICENSE) for details.
