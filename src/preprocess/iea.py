import pandas as pd
import pycountry


COLS = "A:AY"
SKIPROWS = [0, 2]


def preprocess_iea_data(path_to_data, sheet_name, countries, path_to_output):
    countries = list(map(iea_country_names, countries))
    (
        pd
        .read_excel(path_to_data, sheet_name=sheet_name, index_col=0, header=1, skiprows=SKIPROWS, usecols=COLS, na_values="..")
        .loc[countries, :]
        .transpose()
        .rename(columns=sane_country_names)
        .rename(columns=lambda country_name: country_to_country_code(country_name))
        .to_csv(path_to_output, index=True, header=True)
    )


def iea_country_names(country_name):
    if country_name == "Korea, Republic of":
        return "Korea"
    elif country_name == "Slovakia":
        return "Slovak Republic"
    elif country_name == "US":
        return "United States"
    elif country_name == "China":
        return "People's Rep. of China"
    elif country_name == "OECD":
        return "OECD Total"
    elif country_name == "Non-OECD":
        return "Non-OECD Total"
    return country_name


def sane_country_names(country_name):
    if country_name == "Korea":
        return "Korea, Republic of"
    elif country_name == "Slovak Republic":
        return "Slovakia"
    elif country_name == "United States":
        return "US"
    elif country_name == "People's Rep. of China":
        return "China"
    elif country_name == "OECD Total":
        return "OECD"
    elif country_name == "Non-OECD Total":
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
    preprocess_iea_data(
        path_to_data=snakemake.input.iea,
        sheet_name=snakemake.params.sheet_name,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
