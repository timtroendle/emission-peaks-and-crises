import numpy as np
import pandas as pd
import xarray as xr
import pwlf

from crisis import Crisis


def trend(path_to_emissions, path_to_gdp, path_to_population, path_to_carbon_intensity,
          path_to_energy_intensity, crisis, path_to_output):
    time_slice = slice(crisis.pre_from_year - 1, crisis.post_to_year)
    all_ts = {
        "emissions": pd.read_csv(path_to_emissions, index_col=0).loc[time_slice, crisis.country_ids],
        "gdp": pd.read_csv(path_to_gdp, index_col=0).loc[time_slice, crisis.country_ids],
        "population": pd.read_csv(path_to_population, index_col=0).loc[time_slice, crisis.country_ids],
        "carbon-intensity": pd.read_csv(path_to_carbon_intensity, index_col=0).loc[time_slice, crisis.country_ids],
        "energy-intensity": pd.read_csv(path_to_energy_intensity, index_col=0).loc[time_slice, crisis.country_ids],
    }
    t0 = [
        crisis.pre_from_year - 1, # to consider diffs, needs to start 1 year before pre-period
        crisis.from_year - 1, # to consider diffs, needs to start 1 year before crisis
        crisis.to_year, # to consider diffs, needs to start 1 year before post-period
        crisis.post_to_year
    ]
    (
        xr
        .Dataset(data_vars={
            "trend": piecewise_linear_trend(all_ts, t0=t0),
            "r_squared": r_squared(all_ts, t0=t0),
            "p_value": p_value(all_ts, t0=t0)
        })
        .to_netcdf(path_to_output)
    )


def piecewise_linear_trend(all_ts, t0):
    return pd.concat([
        ts_piecewise_linear_trend(data, name, t0)
        for name, data in all_ts.items()
    ]).to_xarray()


def r_squared(all_ts, t0):
    return pd.concat([
        ts_r_squared(data, name, t0)
        for name, data in all_ts.items()
    ]).to_xarray()


def p_value(all_ts, t0):
    return pd.concat([
        ts_p_value(data, name, t0)
        for name, data in all_ts.items()
    ]).to_xarray()


def ts_r_squared(ts, variable_name, t0):
    r_squared = pd.Series( # sum of squares of residuals
        index=ts.columns,
        data=[
            national_piecewise_linear_r_squared(ts.loc[:, country], t0)
            for country in ts.columns
        ]
    )
    return (
        r_squared
        .rename_axis(index="country")
        .to_frame()
        .assign(variable=variable_name)
        .reset_index()
        .set_index(["variable", "country"])
        .iloc[:, 0]
        .rename("r_squared")
    )


def ts_p_value(ts, variable_name, t0):
    trend = pd.DataFrame(
        index=ts.columns,
        data={
            "pre_pre_crisis": [national_piecewise_linear_p_values(ts.loc[:, country], t0)[0]
                               for country in ts.columns],
            "pre_crisis": [national_piecewise_linear_p_values(ts.loc[:, country], t0)[1]
                           for country in ts.columns],
            "crisis": [national_piecewise_linear_p_values(ts.loc[:, country], t0)[2]
                       for country in ts.columns],
            "post_crisis": [national_piecewise_linear_p_values(ts.loc[:, country], t0)[3]
                            for country in ts.columns],
        }
    )
    return (
        trend
        .rename_axis(index="country", columns="period")
        .assign(variable=variable_name)
        .reset_index()
        .set_index(["variable", "country"])
        .stack()
    )


def ts_piecewise_linear_trend(ts, variable_name, t0):
    trend = pd.DataFrame(
        index=ts.columns,
        data={
            "pre_crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[0]
                           for country in ts.columns],
            "crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[1]
                       for country in ts.columns],
            "post_crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[2]
                            for country in ts.columns],
        }
    )
    return (
        trend
        .rename_axis(index="country", columns="period")
        .assign(variable=variable_name)
        .reset_index()
        .set_index(["variable", "country"])
        .stack()
    )


def national_piecewise_linear_fit(ts, t0):
    fit = pwlf.PiecewiseLinFit(ts.index, ts.values)
    fit.fit_with_breaks(t0)
    return fit


def national_piecewise_linear_trend(ts, t0):
    if ts.isna().any():
        return np.nan, np.nan, np.nan
    fit = national_piecewise_linear_fit(ts, t0)
    pre_trend, crisis_trend, post_trend = fit.calc_slopes()
    rel_pre_trend = pre_trend / fit.predict(t0[0])[0]
    rel_crisis_trend = crisis_trend / fit.predict(t0[1])[0]
    rel_post_trend = post_trend / fit.predict(t0[2])[0]
    return rel_pre_trend, rel_crisis_trend, rel_post_trend


def national_piecewise_linear_r_squared(ts, t0):
    if ts.isna().any():
        return np.nan
    fit = national_piecewise_linear_fit(ts, t0)
    return fit.r_squared()


def national_piecewise_linear_p_values(ts, t0):
    if ts.isna().any():
        return np.nan, np.nan, np.nan, np.nan
    fit = national_piecewise_linear_fit(ts, t0)
    return fit.p_values()


if __name__ == "__main__":
    trend(
        path_to_emissions=snakemake.input.emissions,
        path_to_gdp=snakemake.input.gdp,
        path_to_population=snakemake.input.population,
        path_to_carbon_intensity=snakemake.input.carbon_intensity,
        path_to_energy_intensity=snakemake.input.energy_intensity,
        crisis=Crisis.from_config(snakemake.params.crisis, snakemake.params.countries),
        path_to_output=snakemake.output[0]
    )
