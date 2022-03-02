import pandas as pd


def read_maddison_gdp(path_to_data, gdp_sheet_name, pop_sheet_name, country_codes, path_to_output):
    gdp = pd.read_excel(path_to_data, sheet_name=gdp_sheet_name, index_col=0, header=1, skiprows=[7])
    pop = pd.read_excel(path_to_data, sheet_name=pop_sheet_name, index_col=0, header=1, skiprows=[7])
    gdp = gdp * pop.reindex_like(gdp) * 1000
    (
        gdp
        .loc[:, country_codes]
        .to_csv(path_to_output, index=True, header=True)

    )

if __name__ == "__main__":
    read_maddison_gdp(
        path_to_data=snakemake.input.path,
        gdp_sheet_name=snakemake.params.gdp_sheet_name,
        pop_sheet_name=snakemake.params.pop_sheet_name,
        country_codes=snakemake.params.country_codes,
        path_to_output=snakemake.output[0]
    )
