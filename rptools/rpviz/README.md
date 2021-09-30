# rpviz -- Visualize pathways from the RetroPath Suite

## Installation

```sh
conda create --name <myenv>
conda activate <myenv>
conda install -c brsynth -c conda-forge -c bioconda rpviz
``` 
`<myenv>` has to be replaced by whatever meaningful name that will pleased the user.

## Usage

Produce the HTML files
```sh
python -m rpviz.cli <input-folder> <output-folder>
```

To view the pathways, open the `index.html` file outputted in `<out-folder>` using any web browser.

## Example

Using a folder as input:
```sh
conda activate <myenv>
python -m rpviz.cli sample/input/as_dir sample/output/as_dir
```

Using a tar file as input:
```sh
conda activate <myenv>
python -m rpviz.cli sample/input/as_tar.tgz sample/output/as_tar
```

## Command line arguments
```
positional arguments:
  input_rpSBMLs         Input file containing rpSBML files in a 
                        tar archive or a folder.
  output_folder         Output folder to be used. If it does not 
                        exist, an attempt will be made to create 
                        it.It the creation of the folder fails, 
                        IOError will be raised.

optional arguments:
  -h, --help            show this help message and exit
  --debug               Turn on debug instructions
  --template_folder TEMPLATE_FOLDER
                        Path to the folder containing templates
  --cofactor COFACTOR   File listing structures to consider as 
                        cofactors.
  --autonomous_html AUTONOMOUS_HTML
                        Optional file path, if provided will 
                        output an autonomous HTML containing
                        all dependancies.
```

## Input expected by the HTML component

Input file expected by the viewer:
- `network.json`: should contains 2 variable, namely `network` and
`pathways_info`.


## For developpers

### Install
```sh
cd <repository>
conda env create -f environment.yaml -n <dev_env>
conda activate <dev_env>
conda develop -n <dev_env> .
```

Warning: if you do not specify an environment name with -n <dev_env>, then `rpviz-dev` will be used.

### Uninstall
```sh
conda deactivate
conda env remove -n <dev_env>
```

### Test
```sh
conda activate -n <dev_env>
conda install -conda-forge pytest pytest-mock
```

### Generating local documentation
```sh
conda activate -n <dev_env>
conda install -c conda-forge sphinx sphinx_rtd_theme myst-parser 
cd docsource
make clean && make html
```


## JSON objects

### network

Below an overview of the `network` object expected by the JS viewer:

```json
{
    "elements": {
        "nodes": [
            {
                "data": {
                    "id": "string",
                    "path_ids": ["string"],
                    "type": "reaction",
                    "label": "string",
                    "all_labels": ["string"],
                    "svg": null,
                    "xlinks": [{"db_name":  "string", "entity_id": "string", "url": "string"}, ...],
                    "rsmiles": "string",
                    "rule_ids": ["string"],
                    "rxn_template_id": "string",
                    "ec_numbers": ["string"],
                    "thermo_dg_m_gibbs": float,
                    "rule_score": float,
                    "uniprot_ids": [{"UID": {"score": float, "target_ID": "string"}}, ...],
                    "smiles": null,
                    "inchi": null,
                    "inchikey": null,
                    "target_chemical": null,
                    "sink_chemical": null,
                    "thermo_dg_m_formation": null, 
                    "cofactor": null
                }
            },
            ...,
            ...,
            ...,
            {
                "data": {
                    "id": "string",
                    "path_ids": ["string"],
                    "type": "chemical",
                    "label": "string",
                    "all_labels": ["string"],
                    "svg": "string",
                    "xlinks": [{"db_name":  "string", "entity_id": "string", "url": "string"}, ...],
                    "rsmiles": null,
                    "rule_ids": null,
                    "rxn_template_id": null,
                    "ec_numbers": null,
                    "thermo_dg_m_gibbs": null,
                    "rule_score": null,
                    "uniprot_ids": null,
                    "smiles": "string",
                    "inchi": "string",
                    "inchikey": "string",
                    "target_chemical": bool,
                    "sink_chemical": bool,
                    "thermo_dg_m_formation": float, 
                    "cofactor": integer
                }
            },
            ...
            ...
            ...
        ],
        "edges": [
            {
                "data": {
                    "id": "string",
                    "path_ids": ["string"],
                    "source": "string",
                    "target": "string"
                }
            },
            ...
            ...
            ...
        ]
    }
}
```

`network` is composed of 2 types of nodes ("reaction" and "chemical"), and 1 type of edge. What ever the node type,
all the keys ('id', 'path_ids', ...) should be present in each node.


### reaction node

For reaction node, the content should be: 

- `id`, (string), __required value__ -- The canonic reaction SMILES of the reactions. It will be use as the
unique ID of the node. Example: `"id": "[H]OC(=O)C([H])=C([H])C([H])=C([H])C(=O)O[H]>>O=O.[H]Oc1c([H])c([H])c([H])c([H])c1O[H]"`
- `path_ids`: (list of strings), __required values__ -- The list of unique pathway IDs into which the reaction
is involved. It should not contains duplicates. Example: `"path_ids": ["rp_3_1", "rp_2_1", "rp_3_2", "rp_1_1"]`
- `type`, (string), __required value__ -- Should be `"reaction"` for reaction node. It is this value that define
the type of node.
- `label`, (string), __required value__ -- The label to be printed for the node.
- `all_labels`, (list of strings), __optional__ -- All possible labels for the node.
- `svg`, (string), __optional value__ -- SVG depiction of the reaction.
- `xlinks` (list of dictionaries), __optional__ -- Crosslinks to the reaction. Each individual crosslink should
be described in a dictionary having keys: "db_name", "entity_id", "url".   
- `rsmiles`, (string), __optional__ -- The canonical reaction SMILES.
- `rule_ids`, (list of strings), __optional__ -- The reaction rule IDs.
- `ec_numbers`, (list of strings), __optional__ -- The EC numbers.
- `thermo_dg_m_gibbs`, (float), __optional__ -- The dG Gibbs energy of the reaction (in mM concentration context).
- `rule_score`, (float), __optional__ -- The rule score associated to the reaction rule.
- `smiles`, (string), __not used__ -- Value should be `null` (meaningful for chemical node only).
- `inchi`, (string),  __required value__ -- Value should be `null`.
- `inchikey`, (string),  __required value__ -- Value should be `null`.
- `target_chemical`, (bool), __not used__ -- Value should be `null`.
- `sink_chemical`, (bool), __not used__ -- Value should be `null`.
- `thermo_dg_m_formation`, (string), __not used__ -- Value should be `null`.
- `cofactor`, (string), __not used__ -- Value should be `null`.


### chemical node

For chemical node, the content should be:

- `id`, (string), __required value__ -- The InChIKey of the chemical. It will be use as the unique ID of the node.
- `path_ids`: (list of strings), __required values__ -- The list of unique pathway IDs into which the chemical is
involved. It should not contains duplicates. Example: `"path_ids": ["rp_3_1", "rp_2_1", "rp_3_2", "rp_1_1"]`
- `type`, (string), __required value__ -- Should be `"chemical"` for chemical node. It is this value that define the
type of node.
- `label`, (string), __required value__ -- The label to be printed for the node.
- `all_labels`, (list of strings), __optional__ -- All possible labels for the node.
- `svg`, (string),   __optional value__ -- SVG depiction of the entity.
- `xlinks` (list of dictionaries), __optional__ -- Crosslinks to the chemical. Each individual crosslink should
be described in a dictionary having keys: "db_name", "entity_id", "url".
- `rsmiles`, (string), __not used__ -- Value should be `null` (meaningful for reaction node only).
- `rule_ids`, (list of strings), __not used__ -- Value should be `null`.
- `ec_numbers`, (list of strings), __not used__ -- Value should be `null`.
- `thermo_dg_m_gibbs`, (float), __not used__ -- Value should be `null`.
- `rule_score`, (float), __not used__ -- Value should be `null`.
- `smiles`, (string), __required value__ -- The canonic SMILES.
- `inchi`, (string),  __required value__ -- InChI.
- `inchikey`, (string),  __required value__ -- InChIKey.
- `target_chemical`, (bool), __required__ -- Flag to designate the target. Value should be either false (not the
target) or true (it is).
- `sink_chemical`, (bool), __optional value__ -- Flag to designate chemical that are available in the sink. Value
should be either false (not in the sink) or true (it is).
- `thermo_dg_m_formation`, (string), __optional value__ -- The dG of formation of the chemical (in mM concentration context).
- `cofactor`, (string), __required__ -- Flag to designate cofactor chemicals (eg: ATP, NADH, ...). Value should
be either 0 (not a cofactor) or 1 (it is).


### edge

For edge, the content should be:

- `id`, (string), __required value__ -- Some unique string to be used as edge ID. To build such ID, follow the 
convention: `A_B` where `A` is the source node ID (eg a chemical node ID) and `B` is the target node ID (eg 
a reaction node ID).
- `path_ids`: (list of strings), __required values__ -- The list of unique pathway IDs into which the edge is
involved in.
- `source`: the source node ID.
- `target`: the target node ID.


### pathways_info

The pathway_info JSON object purpose is to store information related to each pathways. Notice that order of 
pathway matter: the JS viewer will print the pathways in the order they are given in this object. Here is
the content:
```json
{
    "path_id1": {
        "path_id": "path_id1",
        "nb_steps": interger,
        "node_ids": ["node_id1", "node_id2", ...],
        "edge_ids": ["edge_id1", "edge_id2", ...],
        "scores": {
            "rule_score": float,
            "steps": int,
            "thermo_dg_m_gibbs": float,
            "fba_target_flux": float,
            "global_score": float,
            ...
        },
    },
    "path_id2": {...},
    ...
    ...
    ...
}
```

Where:
- `path_id` is the pathway ID,
- `node_ids` and `edge_ids` lists the nodes and edges involved in this pathway,
- `thermo_dg_m_gibbs` expresses the thermodynamics of the pathway
- `fba_target_flux` is the FBA flux value of the pathway (based on the artificial FBA reaction consuming the target)
- `nb_steps` is the number of reactions involved in the pathway
- `scores` gives the pathway scores.


## Implementation choices

- If a compound cannot be associated to an inchikey (typically, in the case of a generic compound like
`NAD-OR-NADP`), then the MNX ID is used as ID in the network.json file. If no MNX is available, then an
error is logged and the execution is continued.
- rpSBML file name should end by `.xml`
- Target chemical have a SBML ID starting by `TARGET`


## Known bugs and feature requests

### Future release TODO

Fix information into the JSON:
- Fix the URLs for erroneous reaction and chemical crosslinks

Add information into the JSON:
- Build a dedicated section for sequence crosslinks
- Add annotation about the rule diameter to reaction nodes
- Add crosslinks to template reactions used for the rules

## Authors

- Thomas Duigou
- Melchior du Lac 

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details