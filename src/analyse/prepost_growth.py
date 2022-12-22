import pandas as pd
import xarray as xr

from crisis import Crisis


def crises_pre_post_contribution(crises: list, path_to_carbon_intensity: str, path_to_energy_intensity: str,
                                 path_to_energy_and_carbon_intensity: str, path_to_population: str, path_to_gdp: str,
                                 path_to_gdp_and_population: str, paths_to_output: list):
    kaya_factors = (
        pd
        .concat([
            kaya_factor(path_to_carbon_intensity, "carbon-intensity"),
            kaya_factor(path_to_energy_intensity, "energy-intensity"),
            kaya_factor(path_to_energy_and_carbon_intensity, "energy-and-carbon-intensity"),
            kaya_factor(path_to_population, "population"),
            kaya_factor(path_to_gdp_and_population, "gdp-and-population"),
            kaya_factor(path_to_gdp, "gdp"),
        ])
        .rename("kaya_factor")
        .to_xarray()
    )

    prepost = xr.concat([
        xr.concat([
            derive_prepost_growth_rates(kaya_factors, crisis, country_id.item())
            for country_id in kaya_factors.country_id
        ], dim='country_id')
        for crisis in crises
    ], dim='crisis')
    prepost.to_netcdf(paths_to_output.nc)
    prepost.to_series().to_csv(paths_to_output.csv, header=True, index=True)


def derive_prepost_growth_rates(ds, crisis, country_id):
    ds = ds.sel(country_id=country_id)
    period = crisis.national_period(country_id)
    pre = growth_rate(
        initial_level=ds.sel(year=period.pre_from_year),
        final_level=ds.sel(year=period.pre_to_year),
        number_years =(period.pre_to_year - period.pre_from_year)
    )
    post = growth_rate(
        initial_level=ds.sel(year=period.post_from_year),
        final_level=ds.sel(year=period.post_to_year),
        number_years=(period.post_to_year - period.post_from_year)
    )

    pre.coords['period'] = 'pre'
    post.coords['period'] = 'post'
    pre.expand_dims('period')
    post.expand_dims('period')

    prepost = xr.concat([pre, post], dim='period')
    prepost.coords['crisis'] = crisis.slug
    prepost.expand_dims('crisis')
    return prepost


def kaya_factor(path_to_kaya_factor, name):
    kaya_factor = pd.read_csv(path_to_kaya_factor, index_col=0)
    return (
        kaya_factor
        .stack()
        .rename_axis(index=["year", "country_id"])
        .reset_index()
        .assign(factor=name)
        .set_index(["factor", "year", "country_id"])
        .loc[:, 0]
        .rename(name)
    )


def growth_rate(initial_level, final_level, number_years):
    return ((final_level / initial_level) ** (1 / number_years) - 1).rename("growth_rate")


if __name__ == "__main__":
    crises_pre_post_contribution(
        crises=[Crisis.from_config(crisis_slug, snakemake.params.crises[crisis_slug])
                for crisis_slug in snakemake.params.crises],
        path_to_carbon_intensity=snakemake.input.carbon_intensity,
        path_to_energy_intensity=snakemake.input.energy_intensity,
        path_to_energy_and_carbon_intensity=snakemake.input.energy_and_carbon_intensity,
        path_to_population=snakemake.input.population,
        path_to_gdp=snakemake.input.gdp,
        path_to_gdp_and_population=snakemake.input.gdp_and_population,
        paths_to_output=snakemake.output
    )
