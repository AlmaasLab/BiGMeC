![Illustration of BiGMeC pipeline](readme_header.png)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/AlmaasLab/BiGMeC?style=plastic)
![GitHub release date](https://img.shields.io/github/release-date/AlmaasLab/BiGMeC?style=plastic)
[![DOI](https://zenodo.org/badge/273448678.svg)](https://zenodo.org/badge/latestdoi/273448678)

# BiGMeC - Biosynthetic Gene cluster Metabolic pathway Constructor
The BiGMeC pipeline makes a draft reconstruction of the metabolic pathway associated with a non-ribosomal peptide synthetase (NRPS) or polyketide synthase (PKS) biosynthetic gene cluster. The pipeline takes an individual region Genbank (.gbk) file produced by [antiSMASH](https://antismash.secondarymetabolites.org) and produces a JSON-file that can be readily incoporporated into a genome-scale metabolic model using available software such as [cobrapy](cobrapy.readthedocs.io/) or [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/). The pipeline leverages the genome-scale metabolic model of _S. coelicolor_ ([Sco-GEM](https://github.com/SysBioChalmers/Sco-GEM)) as a database for reactions and metabolites.

## Requirements
- python 3 (>=3.5).
- conda or virtualenv, if you want to run the pipeline in a virtual environment (recommended).
- pip package manager. Neccessary to install required python packages. 

## Installation
The BiGMeC pipeline doesn't require any installation, you can simply clone the repository, install a few python packages and immediately run the program. 
### Clone repository
In your command line interface, change to your preferred directory and run
```
git clone https://github.com/AlmaasLab/BiGMeC.git
```
### Make virtual environment
It is recommend to run the pipeline in a virtual environment. You can create this environment with both [virtuealenv](https://virtualenv.pypa.io/en/latest/) or [conda](https://docs.anaconda.com/anaconda/install/), but here we only show how to do it with conda.
From the BiGMeC repository create the environment by running:
```
conda create -n bigmec python=3.6 
```
Activate the new environment:
```
conda activate bigmec
```
### Install required packages
Finally, use pip to install to the required packages listed in _requirements.txt_:
```
pip install -r requirements.txt
```
## Usage
To test the pipeline you can simply run:
```
python bigmec.py
```
This will run the pipeline for the MIBiG gene cluster [BGC0000001](https://mibig.secondarymetabolites.org/repository/BGC0000001/index.html#r1c1) using the corresponding antiSMASH result file [Data/mibig/1.gbk](Data/mibig/1.gkb). The result will be stored in [Data/constructed_pathways/1.json](Data/constructed_pathways/1.json).

To use your own antiSMASH result and store the results in a specified folder run:
```
python bigmec.py -f antismash_data_filename_or_folder -o output_folder
```
Note that antiSMASH provides both a complete and region-specific GenBank-files, and it is the region-specific files that should be used as input for the BiGMeC pipeline.
Further information is provided by running 
```
python bigmec.py -h
```

## Contributors
- [Snorre Sulheim](https://www.sintef.no/en/all-employees/employee/?empId=5675) ([@sulheim](https://github.com/sulheim)), SINTEF Industry, Norway / Norwegian University of Science and Technology, Norway
- [Fredrik Aunaas Fossheim](https://no.linkedin.com/in/fredrik-aunaas-fossheim)([@FredrikFossheim](https://github.com/FredrikFossheim)), Norwegian University of Science and Technology, Norway

## Contribute
Contributions are very welcome, either by raising issues or through pull requests. 
