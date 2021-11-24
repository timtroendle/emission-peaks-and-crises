"""Rules to download raw data."""

rule download_gdp:
    message: "Download GDP data."
    params: url = config["worldbank"]["urls"]["gdp"]
    output: protected("data/automatic/raw-gdp.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_gdp:
    message: "Unzip GDP data."
    input:
        rules.download_gdp.output[0]
    output:
        "build/data/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_3158988.csv"
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    shell: "unzip -o {input} -d build/data"


rule download_population:
    message: "Download population data."
    params: url = config["worldbank"]["urls"]["population"]
    output: protected("data/automatic/raw-population.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_population:
    message: "Unzip population data."
    input:
        rules.download_population.output[0]
    output:
        "build/data/API_SP.POP.TOTL_DS2_en_csv_v2_3158886.csv"
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    shell: "unzip -o {input} -d build/data"


rule download_bp_stats:
    message: "Download BP statistical review of world energy."
    params: url = config["bp"]["url"]
    output: protected("data/automatic/bp-stats-review-all-data.xlsx")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"