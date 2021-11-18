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
        .rename(columns=lambda country_name: country_to_country_code(country_name))
        .to_csv(path_to_output, index=True, header=True)
    )


def bp_country_names(country_name):
    if country_name == "Trinidad and Tobago":
        return "Trinidad & Tobago" # BP uses this unusual way of writing the country name
    elif country_name == "Hong Kong":
        return "China Hong Kong SAR"
    elif country_name == "Korea, Republic of":
        return "South Korea"
    elif country_name == "World":
        return "Total World"
    elif country_name == "OECD":
        return "of which: OECD"
    elif country_name == "Non-OECD":
        return "                 Non-OECD"
    return country_name


def sane_country_names(country_name):
    if country_name == "Trinidad & Tobago":
        return "Trinidad and Tobago"
    elif country_name == "China Hong Kong SAR":
        return "Hong Kong"
    elif country_name == "South Korea":
        return "Korea, Republic of"
    elif country_name == "Total World":
        return "World"
    elif country_name == "of which: OECD":
        return "OECD"
    elif country_name == "                 Non-OECD":
        return "Non-OECD"
    return country_name


def country_to_country_code(country_name):
    if country_name == "World":
        return "WLD"
    elif country_name == "OECD":
        return "OED"
    elif country_name == "Non-OECD":
        return "NOE"
    else:
        return pycountry.countries.lookup(country_name).alpha_3


if __name__ == "__main__":
    preprocess_bp_data(
        path_to_bp_data=snakemake.input.bp,
        sheet_name=snakemake.params.sheet_name,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
