PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"

configfile: "./config/default.yaml"
include: "./rules/download.smk"
include: "./rules/preprocess.smk"
include: "./rules/analyse.smk"


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/global-emissions.png",
        "build/peaker-bp.png",
        "build/prepost-growth-peak-and-decline.png",
        "build/prepost-growth-no-peak-and-decline.png",
        "build/contribution-timeseries/peak-and-decline.png",
        "build/contribution-timeseries/no-peak-and-decline.png",
        "build/contribution-timeseries/all.png",
        "build/logs/test-report.html"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --css=report.css --to html5"
    elif suffix == "pdf":
        return "--css=report.css --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/report.css",
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
        """


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
        emissions = "build/emissions-in-mt-bp.csv",
        contribution_factors = rules.multiplicative_contributions.output.nc
    params:
        countries = COUNTRY_CODES
    output: "build/logs/test-report.html"
    conda: "./envs/test.yaml"
    script: "./tests/test_runner.py"
