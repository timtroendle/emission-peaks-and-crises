import pandas as pd


def decoupling_index(path_to_emissions, path_to_gdp, path_to_output):
    emissions = pd.read_csv(path_to_emissions, index_col=0)
    gdp = pd.read_csv(path_to_gdp, index_col=0)
    gdp_change = gdp.diff().shift(-1) / gdp
    emissions_change = emissions.diff().shift(-1) / emissions
    index = gdp_change / emissions_change
    index.to_csv(path_to_output, index=True, header=True)


if __name__ == "__main__":
    decoupling_index(
        path_to_emissions=snakemake.input.emissions,
        path_to_gdp=snakemake.input.gdp,
        path_to_output=snakemake.output[0]
    )
