
Run: python3 -m SYNBIOCAD_scoring -ttdf reactions_[name].csv -ttsf pathways_[name].csv


Input: you must provide two files for reactions and pathways

Output: mean_stdev_.csv, the score of the pathway(s) listed in the pathways_[name].csv file are in the third column 'Prob1_mean'

Examples of reaction and pathway files provided in the folder:

all pathways used to train the ML model
pathways_all.csv
reactions_all.csv

all (345) pathways producing rosmarinate
pathways_rosmarinate_345.csv
reactions_rosmarinate_345.csv

one pathways producing rosmarinate
pathways_rosmarinate_1.csv
reactions_rosmarinate_1.csv

