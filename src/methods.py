import numpy as np
from scipy import signal
import pandas as pd
import xarray as xr
import pwlf

from crisis import Crisis


def method_analysis(path_to_carbon_emissions, path_to_gdp, path_to_population, path_to_energy_intensity,
                    path_to_carbon_intensity, crisis, path_to_output):
    emissions = pd.read_csv(path_to_carbon_emissions, index_col=0).loc[:, crisis.country_ids]
    gdp = pd.read_csv(path_to_gdp, index_col=0).loc[:, crisis.country_ids]
    population = pd.read_csv(path_to_population, index_col=0).loc[:, crisis.country_ids]
    energy_intensity = pd.read_csv(path_to_energy_intensity, index_col=0).loc[:, crisis.country_ids]
    carbon_intensity = pd.read_csv(path_to_carbon_intensity, index_col=0).loc[:, crisis.country_ids]
    (
        xr
        .Dataset(data_vars={
            "emissions": perform_experiment(ts=emissions, crisis=crisis),
            "gdp": perform_experiment(ts=gdp, crisis=crisis),
            "population": perform_experiment(ts=population, crisis=crisis),
            "carbon_intensity": perform_experiment(ts=carbon_intensity, crisis=crisis),
            "energy_intensity": perform_experiment(ts=energy_intensity, crisis=crisis),
        })
        .to_netcdf(path_to_output)
    )


def perform_experiment(ts, crisis):
    return (
        pd
        .concat([
            calculate_value(method, ts, from_year, to_year, crisis)
            for from_year, to_year in [
                (crisis.post_from_year - 1, crisis.post_to_year - 1),
                (crisis.post_from_year - 1, crisis.post_to_year - 2),
                (crisis.post_from_year - 1, crisis.post_to_year),
                (crisis.post_from_year - 0, crisis.post_to_year)]
            for method in ["old", "linear_trend", "piecewise_linear_trend"]
        ])
        .set_index(["country", "emethod", "range"])
        .to_xarray()
        ["value"]
    )


def calculate_value(method, ts, start_year, end_year, crisis):
    if method == "old":
        return old_method(ts, start_year, end_year, crisis)
    elif method == "linear_trend":
        return linear_trend_method(ts, start_year, end_year, crisis)
    elif method == "piecewise_linear_trend":
        return piecewise_linear_trend_method(ts, start_year, end_year, crisis)
    else:
        raise ValueError(f"Method {method} unknown.")


def old_method(ts, start_year, end_year, crisis):
    ts = ts.loc[start_year:end_year, :]
    value = (ts.iloc[-1] - ts.iloc[0]) / len(ts.index) / ts.iloc[0]
    return (
        value
        .rename("value")
        .rename_axis(index="country")
        .reset_index()
        .assign(emethod="old", range=f"[{start_year}:{end_year}]")
    )


def linear_trend_method(ts, start_year, end_year, crisis):
    ts = ts.loc[start_year:end_year, :]
    trend = ts - signal.detrend(ts, axis=0)
    value = (trend.iloc[-1] - trend.iloc[0]) / len(trend.index) / trend.iloc[0]
    return (
        value
        .rename("value")
        .rename_axis(index="country")
        .reset_index()
        .assign(emethod="linear_trend", range=f"[{start_year}:{end_year}]")
    )


def piecewise_linear_trend_method(ts, start_year, end_year, crisis):
    ts = ts.loc[crisis.pre_from_year - 1:min(crisis.post_to_year, end_year), :]
    t0 = np.array(overwrite_crisis_years(crisis, start_year, end_year))
    value = pd.Series(
        index=ts.columns,
        data=[
            scalar_piecewise_linear_trend(ts.loc[:, country], t0, start_year, end_year)
            for country in ts.columns
        ]
    )
    return (
        value
        .rename("value")
        .rename_axis(index="country")
        .reset_index()
        .assign(emethod="piecewise_linear_trend", range=f"[{start_year}:{end_year}]")
    )


def overwrite_crisis_years(crisis, start_year, end_year):
    crisis_years = [crisis.pre_from_year - 1, crisis.from_year - 1, crisis.post_from_year - 1, crisis.post_to_year]
    start_idx = abs(pd.Series(crisis_years) - start_year).idxmin()
    end_idx = abs(pd.Series(crisis_years) - end_year).idxmin()
    assert abs(start_idx - end_idx) == 1, f"Index identification went wrong for years {start_year} and {end_year}."
    crisis_years[start_idx] = start_year
    crisis_years[end_idx] = end_year
    return crisis_years


def scalar_piecewise_linear_trend(ts, t0, start_year, end_year):
    my_pwlf = pwlf.PiecewiseLinFit(ts.index, ts.values)
    my_pwlf.fit_with_breaks(t0)
    return (my_pwlf.predict(end_year)[0] - my_pwlf.predict(start_year)[0]) / (end_year - start_year) / my_pwlf.predict(start_year)[0]


if __name__ == "__main__":
    method_analysis(
        path_to_carbon_emissions=snakemake.input.emissions,
        path_to_gdp=snakemake.input.gdp,
        path_to_population=snakemake.input.population,
        path_to_carbon_intensity=snakemake.input.carbon_intensity,
        path_to_energy_intensity=snakemake.input.energy_intensity,
        crisis=Crisis.from_config(snakemake.params.crisis),
        path_to_output=snakemake.output[0]
    )
