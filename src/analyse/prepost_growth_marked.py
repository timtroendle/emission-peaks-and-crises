import typing

import pandas as pd
import xarray as xr


def mark_prepost(path_to_prepost: str, peak_and_decline: dict, no_peak_and_decline: dict, peaked_before: dict,
                 paths_to_output: typing.Any):
    ds_prepost = (
        xr
        .open_dataset(path_to_prepost)
        .drop_sel(country_id=["WLD", "OED", "NOE"])
        .sel(factor=["energy-and-carbon-intensity", "gdp-and-population"], drop=True)
    )
    da_peaks = create_peaks_data_array(peak_and_decline, no_peak_and_decline, peaked_before, ds_prepost)
    ds_prepost["peaked"] = da_peaks
    ds_prepost["diffs"] = ds_prepost.sel(period="post")["growth_rate"] - ds_prepost.sel(period="pre")["growth_rate"]

    ds_prepost.to_netcdf(paths_to_output.nc)
    ds_prepost.to_dataframe().to_csv(paths_to_output.csv, header=True, index=True)


def create_peaks_data_array(peak_and_decline: dict, no_peak_and_decline: dict, peaked_before: dict,
                            prepost: xr.DataArray) -> xr.DataArray:
    crises = set(peak_and_decline.keys()) & set(no_peak_and_decline.keys())
    return (
        pd
        .concat([
            create_peaks_for_crisis(crisis, peak_and_decline[crisis], no_peak_and_decline[crisis],
                                    peaked_before[crisis], prepost)
            for crisis in crises])
        .to_xarray()
    )


def create_peaks_for_crisis(crisis_name: str, peak_and_decline: list, no_peak_and_decline: list,
                            peaked_before: list, prepost: xr.DataArray) -> pd.Series:
    if not peaked_before:
        peaked_before = list()
    assert set(peak_and_decline) & set(no_peak_and_decline) == set(), crisis_name
    assert set(no_peak_and_decline) & set(peaked_before) == set(), crisis_name
    peaks = pd.Series(
        index=prepost.country_id,
        data=pd.Categorical(
            [pd.NA] * len(prepost.country_id),
            categories=["peak and decline", "no peak and decline", "peaked before"]),
        name="peaked"
    )
    for country_id in peak_and_decline:
        peaks[country_id] = "peak and decline"
    for country_id in no_peak_and_decline:
        peaks[country_id] = "no peak and decline"
    for country_id in peaked_before:
        peaks[country_id] = "peaked before"
    peaks = (
        peaks
        .rename_axis(index="country_id")
        .to_frame()
        .assign(crisis=crisis_name)
        .reset_index()
        .set_index(["crisis", "country_id"])["peaked"]
    )
    return peaks


if __name__ == "__main__":
    mark_prepost(
        path_to_prepost=snakemake.input.prepost,
        peak_and_decline=snakemake.params.peak_and_decline,
        no_peak_and_decline=snakemake.params.no_peak_and_decline,
        peaked_before=snakemake.params.peaked_before,
        paths_to_output=snakemake.output
    )
