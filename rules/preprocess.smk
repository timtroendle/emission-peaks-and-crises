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
        src = "src/preprocess/bp.py",
        bp = rules.download_bp_stats.output[0]
    params:
        sheet_name = config["data-sources"]["bp"]["sheetnames"]["energy"],
        countries = config["countries"]
    output: "build/energy-in-ej.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/bp.py"


rule emissions_bp:
    message: "Preprocess BP emissions data."
    input:
        src = "src/preprocess/bp.py",
        bp = rules.download_bp_stats.output[0]
    params:
        sheet_name = config["data-sources"]["bp"]["sheetnames"]["emissions"],
        countries = config["countries"]
    output: "build/emissions-in-mt-bp.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/bp.py"


rule emissions_iea:
    message: "Preprocess IEA emissions data."
    input:
        src = "src/preprocess/iea.py",
        iea = config["data-sources"]["iea"]["path"]
    params:
        sheet_name = "GHG FC",
        countries = config["countries"]
    output: "build/emissions-in-mt-iea.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/iea.py"


rule population:
    message: "Preprocess Worldbank population data."
    input:
        src = "src/preprocess/worldbank.py",
        path = rules.unzip_population.output[0]
    params:
        country_codes = COUNTRY_CODES
    output: "build/population.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/worldbank.py"


rule gdp:
    message: "Preprocess Worldbank GDP data."
    input:
        src = "src/preprocess/worldbank.py",
        path = rules.unzip_gdp.output[0]
    params:
        country_codes = COUNTRY_CODES
    output: "build/gdp-in-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/worldbank.py"


rule energy_intensity:
    message: "Divide {input.df1} / {input.df2}"
    input:
        src = "src/preprocess/divide.py",
        df1 = rules.energy.output[0],
        df2 = rules.gdp.output[0]
    output: "build/energy-intensity-in-ej-per-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule carbon_intensity:
    message: "Divide {input.df1} / {input.df2}"
    input:
        src = "src/preprocess/divide.py",
        df1 = "build/emissions-in-mt-bp.csv",
        df2 = rules.energy.output[0]
    output: "build/carbon-intensity-in-mt-per-ej.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule energy_and_carbon_intensity:
    message: "Divide {input.df1} / {input.df2}."
    input:
        src = "src/preprocess/divide.py",
        df1 = "build/emissions-in-mt-bp.csv",
        df2 = rules.gdp.output[0]
    output: "build/energy-and-carbon-intensity-in-mt-per-usd.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"


rule gdp_per_capita:
    message: "Divide {input.df1} / {input.df2}"
    input:
        src = "src/preprocess/divide.py",
        df1 = rules.gdp.output[0],
        df2 = rules.population.output[0]
    output: "build/gdp-in-usd-per-capita.csv"
    conda: "../envs/default.yaml"
    script: "../src/preprocess/divide.py"
