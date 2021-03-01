import xarray as xr

TREND_VARIABLE = "trend"
R_SQUARED_VARIABLE = "r_squared"


def csv_results(path_to_trend, path_to_trend_change, path_to_r_squared):
    ds = xr.open_dataset(path_to_trend)
    write_trend_change(ds[TREND_VARIABLE], path_to_trend_change)
    write_r_squared(ds[R_SQUARED_VARIABLE], path_to_r_squared)


def write_trend_change(da, path_to_output):
    da = da.copy() * 100
    (
        (da.sel(period="post_crisis") - da.sel(period="pre_crisis"))
        .to_series()
        .unstack()
        .T
        .to_csv(path_to_output, index=True, header=True)
    )


def write_r_squared(da, path_to_output):
    (
        da
        .to_series()
        .unstack()
        .T
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    csv_results(
        path_to_trend=snakemake.input.trend,
        path_to_trend_change=snakemake.output.trend_change,
        path_to_r_squared=snakemake.output.r_squared
    )
