import pandas as pd


def contributions(path_to_carbon_intensity, path_to_energy_intensity,
                  path_to_population, path_to_gdp, path_to_output):
    (
        pd
        .concat([
            read_contribution(path_to_carbon_intensity, "carbon-intensity"),
            read_contribution(path_to_energy_intensity, "energy-intensity"),
            read_contribution(path_to_population, "population"),
            read_contribution(path_to_gdp, "gdp"),
        ])
        .to_xarray()
        .rename("contributions")
        .to_netcdf(path_to_output)
    )


def read_contribution(path_to_contribution, name):
    return (
        pd
        .read_csv(path_to_contribution, index_col=0)
        .stack()
        .rename_axis(index=["year", "country_id"])
        .reset_index()
        .assign(factor=name)
        .set_index(["factor", "year", "country_id"])
        .loc[:, 0]
        .rename(name)
    )


if __name__ == "__main__":
    contributions(
        path_to_carbon_intensity=snakemake.input.carbon_intensity,
        path_to_energy_intensity=snakemake.input.energy_intensity,
        path_to_population=snakemake.input.population,
        path_to_gdp=snakemake.input.gdp,
        path_to_output=snakemake.output[0]
    )
