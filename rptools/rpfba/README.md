# rpFBA

Perform FBA on a single or collection of SBML files containing heterologous pathways, as tar.xz archives. The package performs the following steps:
    1) it merges a user defined GEM SBML model to a given heterologous pathway.
    2) it performs FBA using the [cobrapy](https://opencobra.github.io/cobrapy/) package using a user defined method that includes FBA, parsimonious FBA or fraction of optimum of another reaction. For the first two, the user must know the reaction name that the model will optimise to, while the latter the user must provide the target reaction (default: 'rxn_target') but also another reaction that will be restricted (default: 'biomass'). The first step involves performing FBA using the "source" reaction as the objective. Then the flux of that reaction has its upper and lower bounds set to the same value, determined as a fraction of its FBA flux value. Thereafter the objective is set to the initial target reaction and FBA is performed once again. The tool uses the [FBC](https://co.mbine.org/specifications/sbml.level-3.version-1.fbc.version-2.release-1) package to manage the objective and flux bounds.

NOTE: In order to FBA works correctly, some of chemical species have to be ignored. These species are selected according to the following criteria:
* pathway species has not been found in the model (neither by its ID nor its InChIKey),
* the species is not the target, and
* the species is only consumed or produced with the heterologue pathway.


## Input

Required:
* **pathway_file**: (string) Path to the pathway file (rpSBML)
* **model_file**: (string) Path to the GEM SBML model
* **compartment_id**: (string, e.g. cytoplasm) ID of the compartment that contains the chemical species involved in the heterologous pathway
* **out_file**: (string) Path to the ouput upgraded pathway file

Advanced options:
* **--sim**: (string, default='fraction') Valid options include: 'fraction', 'fba', 'pfba'. The type of constraint based modelling method
* **--objective_rxn_id**: (string, default=rxn_target) Reaction ID to optimise
* **--biomass_rxn_id**: (string, default='biomass') Biomass reaction ID. Note: Only for 'fraction' simulation
* **--fraction_of**: (float, default=0.75) Portion of the maximal flux used to set the maximal and minimal bounds for the source reaction of the 'fraction' simulation type
* **--merge**: (boolean, default=False) Return the merged GEM+heterologous pathway SBML or only the heterologous pathway SBML files
* **--ignore_orphan_species**: (string, default=True) Ignore metabolites that are only consumed or produced

## Output

* **output**: (string) Path to the output file


## Install
Please see `rptool` documentation.

## Run

### rpFBA process
**From Python code**
```python
from rptools.rplibs import rpPathway
from rptools.rpfba import runFBA

pathway = rpPathway(infile='tests/rpfba/data/sets/measured_3/B.xml')

results = runFBA(
    pathway=pathway,
    gem_sbml_path='tests/rpfba/data/sets/measured_3/e_coli_iJ01366.xml',
    compartment_id='MNXC3',
    biomass_rxn_id='biomass',
    objective_rxn_id='rxn_target',
    sim_type='fraction',
    fraction_coeff=0.75,
    merge=False,
    ignore_orphan_species=True
)

pathway.get_fba_biomass().get('value')
0.7638744755010194
pathway.get_fba_biomass()
{'value': 0.7638744755010194, 'units': 'milimole / gDW / hour'}
```
**From CLI**
```sh
python -m rptools.rpfba <pathway_rpsbml> <model_gem_sbml> <compartment_id> <outfile>
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

* **Joan Hérisson**
* **Melchior du Lac**

## Acknowledgments

* Thomas Duigou


## Licence
rpFBA is released under the MIT licence. See the LICENCE file for details.
