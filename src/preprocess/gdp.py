import pandas as pd


def gdp(path_to_worldbank_gdp, path_to_maddison_gdp, path_to_output):
    worldbank = pd.read_csv(path_to_worldbank_gdp, index_col=0)
    maddison = pd.read_csv(path_to_maddison_gdp, index_col=0)
    (
        maddison
        .reindex_like(worldbank)
        .fillna(worldbank)
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    gdp(
        path_to_worldbank_gdp=snakemake.input.worldbank,
        path_to_maddison_gdp=snakemake.input.maddison,
        path_to_output=snakemake.output[0]
    )
