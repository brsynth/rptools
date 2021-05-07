# rplibs

Libraries for rpTools, contains:
* rpSBML
* inchikeyMIRIAM

## rpSBML
Defines SBML structure with additional fields relative to [RetroPath2](https://github.com/brsynth/RetroPath2-wrapper) objects.

<!-- ### Prerequisites
* Python 3 with the following modules:
    * python-libsbml
    * [RDKit](https://www.RDKit.org) -->


### Install
rpLibs is part of rpTools suite:
#### From Conda
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda rptools
```

### Test
Please follow instructions below ti run tests:
```
cd tests
pytest -v
```
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)).

## inchikeyMIRIAM
Uses the rrCache to parse an SBML file to find all the chemical species, and try to recover the inchikey and add it to the MIRIAM annotation.



## Authors

* **Melchior du Lac**
* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou


## Licence
brs_libs is released under the MIT licence. See the LICENCE file for details.
