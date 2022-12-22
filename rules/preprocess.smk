"""Rules to preprocess Worldbank and BP data."""
import pycountry


def country_to_country_code(country_name):
    if country_name == "World":
        return "WLD"
    elif country_name == "OECD":
        return "OED"
    elif country_name == "Non-OECD":
        return "NOE"
    else:
        return pycountry.countries.lookup(country_name).alpha_3


COUNTRY_CODES = sorted([
    country_to_country_code(country_name)
    for country_name in config["countries"]
])


rule energy:
    message: "Preprocess BP energy consumption data."
    input:
        bp = rules.download_bp_stats.output[0]
    params:
        sheet_name = config["data-sources"]["bp"]["sheet_names"]["energy"],
        countries = config["countries"]
    output: "build/data/energy-in-ej.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/bp.py"


rule emissions_bp:
    message: "Preprocess BP emissions data."
    input:
        bp = rules.download_bp_stats.output[0]
    params:
        sheet_name = config["data-sources"]["bp"]["sheet_names"]["emissions"],
        countries = config["countries"]
    output: "build/data/emissions-in-mt-bp.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/bp.py"


rule emissions_iea:
    message: "Preprocess IEA emissions data."
    input:
        iea = config["data-sources"]["iea"]["path"]
    params:
        sheet_name = "GHG FC",
        countries = config["countries"]
    output: "build/data/emissions-in-mt-iea.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/iea.py"


rule population:
    message: "Preprocess Worldbank population data."
    input:
        path = rules.unzip_population.output[0]
    params:
        country_codes = COUNTRY_CODES
    output: "build/data/population.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/worldbank.py"


rule worldbank_gdp:
    message: "Preprocess Worldbank GDP data."
    input:
        path = rules.unzip_gdp.output[0]
    params:
        country_codes = COUNTRY_CODES
    output: "build/data/worldbank-gdp-in-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/worldbank.py"


rule maddison_gdp:
    message: "Preprocess Maddison GDP data."
    input:
        path = rules.download_maddison_gdp.output[0]
    params:
        country_codes = config["data-sources"]["maddison"]["countries"],
        gdp_sheet_name = config["data-sources"]["maddison"]["sheet_names"]["gdp_per_capita"],
        pop_sheet_name = config["data-sources"]["maddison"]["sheet_names"]["population"]
    output: "build/data/maddison-gdp-in-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/maddison.py"


rule gdp:
    message: "Use selected countries from Maddison GDP, otherwise use Worldbank GDP."
    input:
        worldbank = rules.worldbank_gdp.output[0],
        maddison = rules.maddison_gdp.output[0]
    output: "build/data/gdp-in-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/gdp.py"


rule energy_intensity:
    message: "Divide {input.df1} / {input.df2}"
    input:
        df1 = rules.energy.output[0],
        df2 = rules.gdp.output[0]
    output: "build/data/energy-intensity-in-ej-per-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule carbon_intensity:
    message: "Divide {input.df1} / {input.df2}"
    input:
        df1 = "build/data/emissions-in-mt-bp.csv",
        df2 = rules.energy.output[0]
    output: "build/data/carbon-intensity-in-mt-per-ej.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule energy_and_carbon_intensity:
    message: "Divide {input.df1} / {input.df2}."
    input:
        df1 = "build/data/emissions-in-mt-bp.csv",
        df2 = rules.gdp.output[0]
    output: "build/data/energy-and-carbon-intensity-in-mt-per-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule gdp_per_capita:
    message: "Divide {input.df1} / {input.df2}"
    input:
        df1 = rules.gdp.output[0],
        df2 = rules.population.output[0]
    output: "build/data/gdp-in-usd-per-capita.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule flag:
    message: "Preprocess flag of {wildcards.country}."
    input:
        flag = rules.download_flag.output[0],
    params: height = 200
    output: "build/flag/{country}.png"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/flag.py"
