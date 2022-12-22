from snakemake.utils import min_version

configfile: "./config/default.yaml"
include: "./rules/download.smk"
include: "./rules/preprocess.smk"
include: "./rules/analyse.smk"
min_version("7.17")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/figures/global-emissions.png",
        "build/figures/peaker-bp.png",
        "build/figures/prepost-growth-peak-and-decline.png",
        "build/figures/prepost-growth-no-peak-and-decline.png",
        "build/analysis/prepost-growth-tests.csv",
        "build/figures/contribution-timeseries/peak-and-decline.png",
        "build/figures/contribution-timeseries/no-peak-and-decline.png",
        "build/figures/contribution-timeseries/non-crisis-peaker.png",
        "build/figures/contribution-timeseries/all.png",
        "build/logs/test-report.html"


rule push:
    message: "Package, zip, and move entire build."
    params: push_directory = config["push-directory"]
    shell:
        """
        zip -r {params.push_directory}/kaya-results-$(date -Idate).zip build
        """


rule build_data_publication:
    message: "Package data publication."
    input:
        "build/analysis/multiplicative-contributions.csv",
        "build/analysis/multiplicative-contributions.nc",
        "build/analysis/prepost-growth-rates.csv",
        "build/analysis/prepost-growth-rates.nc",
        "src/doc/README.md"
    output:
        "build/data-publication.zip"
    shell: "zip -rj {output} {input}"



rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    message: "Run tests"
    input:
        test_dir = "tests",
        tests = map(str, Path("tests").glob("**/test_*.py")),
        emissions = "build/data/emissions-in-mt-bp.csv",
        contribution_factors = rules.multiplicative_contributions.output.nc
    params:
        countries = COUNTRY_CODES
    output: "build/logs/test-report.html"
    conda: "./envs/test.yaml"
    script: "./tests/test_runner.py"
