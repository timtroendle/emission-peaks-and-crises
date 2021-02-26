import xarray as xr


def trend_change(path_to_trend, path_to_output):
    da = xr.open_dataset(path_to_trend)["trend"] * 100
    (
        (da.sel(period="post_crisis") - da.sel(period="pre_crisis"))
        .to_series()
        .unstack()
        .T
        .to_csv(path_to_output, index=True, header=True)
    )


if __name__ == "__main__":
    trend_change(
        path_to_trend=snakemake.input.trend,
        path_to_output=snakemake.output[0]
    )
