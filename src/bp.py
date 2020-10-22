import pandas as pd
import pycountry


COLS = "A:BD"
SKIPROWS = [0, 1]


def preprocess_bp_data(path_to_bp_data, sheet_name, countries, path_to_output):
    countries = list(map(bp_country_names, countries))
    (
        pd
        .read_excel(path_to_bp_data, sheet_name=sheet_name, index_col=0, skiprows=SKIPROWS, usecols=COLS)
        .loc[countries, :]
        .transpose()
        .rename(columns=sane_country_names)
        .rename(columns=lambda country_name: pycountry.countries.lookup(country_name).alpha_3)
        .to_csv(path_to_output, index=True, header=True)
    )


def bp_country_names(country_name):
    if country_name == "Trinidad and Tobago":
        return "Trinidad & Tobago" # BP uses this unusual way of writing the country name
    return country_name


def sane_country_names(country_name):
    if country_name == "Trinidad & Tobago":
        return "Trinidad and Tobago"
    return country_name


if __name__ == "__main__":
    preprocess_bp_data(
        path_to_bp_data=snakemake.input.bp,
        sheet_name=snakemake.params.sheet_name,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
