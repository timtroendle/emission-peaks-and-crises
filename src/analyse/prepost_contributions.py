import xarray as xr

from crisis import Crisis


def crises_pre_post_contribution(crises, path_to_multiplicative_contributions, paths_to_output):
    ds = xr.open_dataset(path_to_multiplicative_contributions)

    prepost = (
        xr
        .concat(
            [derive_prepost_contribution_factors(ds, crisis) for crisis in crises],
            dim='crisis'
        )
    )
    prepost.to_netcdf(paths_to_output.nc)
    prepost["contributions"].to_series().to_csv(paths_to_output.csv, header=True, index=True)


def derive_prepost_contribution_factors(ds, crisis):
    pre = ds.sel(year=range(crisis.pre_from_year, crisis.pre_to_year + 1)).mean(dim="year")
    post = ds.sel(year=range(crisis.post_from_year, crisis.post_to_year + 1)).mean(dim="year")

    pre.coords['period'] = 'pre'
    post.coords['period'] = 'post'
    pre.expand_dims('period')
    post.expand_dims('period')

    prepost = xr.concat([pre, post], dim='period')
    prepost.coords['crisis'] = crisis.slug
    prepost.expand_dims('crisis')
    return prepost


if __name__ == "__main__":
    crises_pre_post_contribution(
        crises=[Crisis.from_config(crisis_slug, snakemake.params.crises[crisis_slug])
                for crisis_slug in snakemake.params.crises],
        path_to_multiplicative_contributions=snakemake.input.contributions,
        paths_to_output=snakemake.output
    )
