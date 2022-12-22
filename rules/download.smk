"""Rules to download raw data."""

rule download_gdp:
    message: "Download GDP data."
    params: url = config["data-sources"]["worldbank"]["urls"]["gdp"]
    output: protected("data/automatic/raw-gdp.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_gdp:
    message: "Unzip GDP data."
    input:
        rules.download_gdp.output[0]
    output:
        "build/data/API_NY.GDP.MKTP.KD_DS2_en_csv_v2_4687500.csv"
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    shell: "unzip -o {input} -d build/data"


rule download_population:
    message: "Download population data."
    params: url = config["data-sources"]["worldbank"]["urls"]["population"]
    output: protected("data/automatic/raw-population.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_population:
    message: "Unzip population data."
    input:
        rules.download_population.output[0]
    output:
        "build/data/API_SP.POP.TOTL_DS2_en_csv_v2_4687738.csv"
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    shell: "unzip -o {input} -d build/data"


rule download_bp_stats:
    message: "Download BP statistical review of world energy."
    params: url = config["data-sources"]["bp"]["url"]
    output: protected("data/automatic/bp-stats-review-all-data.xlsx")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule download_maddison_gdp:
    message: "Download Maddison GDP data."
    params: url = config["data-sources"]["maddison"]["url"]
    output: protected("data/automatic/mpd2020.xlsx")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


def flag_url(wildcards):
    country_code = pycountry.countries.lookup(wildcards.country).alpha_3
    country_unicode = config["flags"][country_code]
    return config["data-sources"]["flags"]["url"].format(country_unicode=country_unicode)


rule download_flag:
    message: "Download flag of {wildcards.country}."
    params: url = flag_url
    output: "data/automatic/flags/{country}.png"
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"
