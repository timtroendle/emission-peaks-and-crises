import numpy as np
import pandas as pd


def preprocess_worldbank_data(path_to_data, country_codes, path_to_output):
    return (
        pd
        .read_csv(path_to_data, skiprows=4)
        .iloc[:, :-1]
        .set_index("Country Code")
        .transpose()
        .drop(["Indicator Name", "Indicator Code", "Country Name"], axis="index")
        .rename(index=lambda x: int(x))
        .rename_axis(index="year", columns="country_code")
        .astype(np.float64)
        .loc[:, country_codes]
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    preprocess_worldbank_data(
        path_to_data=snakemake.input.path,
        country_codes=snakemake.params.country_codes,
        path_to_output=snakemake.output[0]
    )
