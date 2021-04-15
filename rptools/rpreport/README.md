# rpreport

## Installation
### Local conda environment setup :
We strongly recommend you to use conda package manager, and follow those steps:  
`conda create -n rpreport_env python=3`  
`conda activate rpreport_env`  
`conda install -c brsynth -c conda-forge -c bioconda rptools`  (rpSBML installation, see https://github.com/brsynth/rpTools)
`pip install dictor`  (see https://github.com/perfecto25/dictor)

### Usage
In the conda environment
`python report.py <archive.tar> <outputfolder>`
where <archive.tar> is an archive containing RPSBML files

User `-v` for verbose mode

### Pycharm local settings
Configuration de PyCharm pour utiliser conda en créant le projet :
https://www.jetbrains.com/help/pycharm/conda-support-creating-conda-virtual-environment.html
TODO : configurer git avec projet. Apparement c’est déjà fait, PyCharm lit les infos git du dépôt précédemment créé.
