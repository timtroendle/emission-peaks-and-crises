# Kaya

Analysing how Kaya factors change during crises.

This repository contains the entire scientific project. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml --no-default-packages

## Run the analysis

    snakemake --profile profiles/default

This will run all analysis steps to reproduce results.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --profiles profiles/default --rulegraph | dot -Tpdf > dag.pdf

## Run the tests

    snakemake --profile profiles/default test

## Repo structure

* `src`: contains the Python source code
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)
* `profiles`: Snakemake execution profiles

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
