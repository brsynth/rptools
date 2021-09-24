# rpreport
Generates HTML pages to explore the main characteristics (thermodynamics,
fluxes, number of metabolic steps, reaction rule score) of pathways predicted
with RetroPath suite


## Installation
Please see `rptools` documentation.

## Usage
rpreport [-h] [--log ARG] [--log_file LOG_FILE] [--silent] [--version] [-d] [--dev] [-v] [--standalone] source_path output_folder

Required:
* **source_path**: (string) Path to a tar archive (default) or folder (using '-d' option) containing rpSBML file(s).
* **output_folder**: (string) Output folder where report file(s) will be generated.

Optional:
* **-h, --help**: show this help message and exit
* **--log ARG, -l ARG**: Adds a console logger for the specified level (default:
                       error)
* **--log_file LOG_FILE**: Filename where to put logs
* **--silent, -s**: run rpreport silently
* **--version**: show the version number and exit
* **-d, --input_dir**: source_path is a folder containing standalone rpSBML
                       file(s).
* **--dev**: For dev purpose only : create supplementary files into
                       a dev folder
* **-v, --verbose**: Turn on console verbose mode.
* **--standalone**: if set will output an autonomous HTML containing all
                       css and js files dependencies.
