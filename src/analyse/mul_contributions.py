import pandas as pd


def contributions(path_to_carbon_intensity, path_to_energy_intensity,
                  path_to_population, path_to_gdp, path_to_output_nc, path_to_output_csv):
    contributions = (
        pd
        .concat([
            contribution(path_to_carbon_intensity, "carbon-intensity"),
            contribution(path_to_energy_intensity, "energy-intensity"),
            contribution(path_to_population, "population"),
            contribution(path_to_gdp, "gdp"),
        ])
        .rename("contributions")
    )
    contributions.to_csv(path_to_output_csv, index=True, header=True)
    (
        contributions
        .to_xarray()
        .to_netcdf(path_to_output_nc)
    )


def contribution(path_to_kaya_factor, name, n=1):
    kaya_factor = pd.read_csv(path_to_kaya_factor, index_col=0)
    contribution = kaya_factor / kaya_factor.shift(n)
    return (
        contribution
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
        path_to_output_nc=snakemake.output.nc,
        path_to_output_csv=snakemake.output.csv
    )
