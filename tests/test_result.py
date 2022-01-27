import pandas as pd
import xarray as xr

EPSILON = 0.01 # error is smaller than 1 percentage point
SHIFT_YEARS = 1 # contribution factors are shifted by 1 year

FACTORS = ["population", "gdp", "carbon-intensity", "energy-intensity"]


def test_positive_emissions(emissions: pd.DataFrame, country_id: str):
    assert (emissions.loc[:, country_id].dropna() > 0).all()


def test_positive_contribution_factors(contribution_factors: xr.Dataset, contribution_factor: str, country_id: str):
    assert (
        (contribution_factors.sel(country_id=country_id, factor=contribution_factor).dropna(dim="year") > 0)
        .all()
        .item()
    )


def test_decomposition(emissions: pd.DataFrame, contribution_factors: xr.Dataset, country_id: str):
    # This tests equation (2) in (Ang 2005)
    emissions = emissions.loc[:, country_id]
    contribution_factors = contribution_factors.sel(country_id=country_id, factor=FACTORS).dropna("year")

    d_tot = (emissions / emissions.shift(SHIFT_YEARS)).dropna()
    d_sum = contribution_factors.prod("factor", skipna="False").to_series().reindex_like(d_tot).dropna()

    d_tot = d_tot.reindex_like(d_sum)

    diff = (d_tot - d_sum).abs()
    assert (diff < EPSILON).all()
