# Emission peaks and crises

This study is investigating the relationship between emission peaks and economic crises.

This repository contains the entire scientific project. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

[![article](https://img.shields.io/badge/article-doi-blue)](https://doi.org/10.1038/s43247-023-00687-8)
[![data](https://img.shields.io/badge/data-10.5281/zenodo.7474120-blue)](https://doi.org/10.5281/zenodo.7474120)
[![workflow archive](https://img.shields.io/badge/workflow-10.5281/zenodo.7477484-blue)](https://doi.org/10.5281/zenodo.7477484)

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml --no-default-packages

## Run the analysis

    snakemake

This will run all analysis steps to reproduce results.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run the tests

    snakemake test

## Repo structure

* `build`: will contain all results (does not exist initially)
* `config`: configurations used in the study
* `envs`: contains execution environments
* `profiles`: Snakemake execution profiles
* `rules`: Snakemake rules
* `src`: contains the Python source code
* `tests`: contains the test code

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
