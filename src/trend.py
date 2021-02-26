import pandas as pd
import xarray as xr
import pwlf

from crisis import Crisis


def trend(path_to_emissions, path_to_gdp, path_to_population, path_to_carbon_intensity,
          path_to_energy_intensity, crisis, path_to_output):
    all_ts = {
        "emissions": pd.read_csv(path_to_emissions, index_col=0).loc[:, crisis.country_ids],
        "gdp": pd.read_csv(path_to_gdp, index_col=0).loc[:, crisis.country_ids],
        "population": pd.read_csv(path_to_population, index_col=0).loc[:, crisis.country_ids],
        "carbon-intensity": pd.read_csv(path_to_carbon_intensity, index_col=0).loc[:, crisis.country_ids],
        "energy-intensity": pd.read_csv(path_to_energy_intensity, index_col=0).loc[:, crisis.country_ids],
    }
    (
        xr
        .Dataset(data_vars={
            "trend": piecewise_linear_trend(all_ts, crisis=crisis),
            "goodness_of_fit": goodness_of_fit(all_ts, crisis=crisis)
        })
        .to_netcdf(path_to_output)
    )


def piecewise_linear_trend(all_ts, crisis):
    return pd.concat([
        ts_piecewise_linear_trend(data, name, crisis)
        for name, data in all_ts.items()
    ]).to_xarray()


def goodness_of_fit(all_ts, crisis):
    return pd.concat([
        ts_goodness_of_fit(data, name, crisis)
        for name, data in all_ts.items()
    ]).to_xarray()


def ts_goodness_of_fit(ts, variable_name, crisis):
    ts = ts.loc[crisis.pre_from_year - 1:crisis.post_to_year, :]
    t0 = [crisis.pre_from_year - 1, crisis.from_year - 1, crisis.post_from_year - 1, crisis.post_to_year]
    goodness = pd.Series(
        index=ts.columns,
        data=[
            national_piecewise_linear_trend(ts.loc[:, country], t0)[3]
            for country in ts.columns
        ]
    )
    return (
        goodness
        .rename_axis(index="country")
        .to_frame()
        .assign(variable=variable_name)
        .reset_index()
        .set_index(["variable", "country"])
        .iloc[:, 0]
        .rename("goodness_of_fit")
    )


def ts_piecewise_linear_trend(ts, variable_name, crisis):
    ts = ts.loc[crisis.pre_from_year - 1:crisis.post_to_year, :]
    t0 = [crisis.pre_from_year - 1, crisis.from_year - 1, crisis.post_from_year - 1, crisis.post_to_year]
    trend = pd.DataFrame(
        index=ts.columns,
        data={
            "pre_crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[0] for country in ts.columns],
            "crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[1] for country in ts.columns],
            "post_crisis": [national_piecewise_linear_trend(ts.loc[:, country], t0)[2] for country in ts.columns],
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


def national_piecewise_linear_trend(ts, t0):
    fit = pwlf.PiecewiseLinFit(ts.index, ts.values)
    goodness = fit.fit_with_breaks(t0)
    pre_trend = annual_relative_change(fit, t0[0])
    crisis_trend = annual_relative_change(fit, t0[1])
    post_trend = annual_relative_change(fit, t0[2])
    return (pre_trend, crisis_trend, post_trend, goodness)


def annual_relative_change(fit, year):
    return (fit.predict(year + 1)[0] - fit.predict(year)[0]) / fit.predict(year)[0]


if __name__ == "__main__":
    trend(
        path_to_emissions=snakemake.input.emissions,
        path_to_gdp=snakemake.input.gdp,
        path_to_population=snakemake.input.population,
        path_to_carbon_intensity=snakemake.input.carbon_intensity,
        path_to_energy_intensity=snakemake.input.energy_intensity,
        crisis=Crisis.from_config(snakemake.params.crisis),
        path_to_output=snakemake.output[0]
    )
