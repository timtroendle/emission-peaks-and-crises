import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns


COLOR_PALETTE = [ # Nature colors
    "#E64B35",
    "#4DBBD5",
    "#00A087",
    "#3C5488",
    "#F39B7F",
    "#8491B4",
    "#91D1C2",
    "#DC0000",
    "#7E6148",
    "#B09C85",
    "#00A087",
]
sns.set_palette(COLOR_PALETTE)


def timeseries(path_to_contributions, path_to_emissions, country_id, from_year, to_year, path_to_output):
    ds = (
        read_data(path_to_contributions, path_to_emissions)
        .sel(country_id=country_id, year=range(from_year, to_year + 1))
    )
    reference_year = ds.year[0].item()
    reference_value = ds.emissions.sel(year=reference_year).item()
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    ax.plot(
        ds.year,
        ds.emissions / reference_value,
        label="$\mathrm{CO_2}$ emissions",
        linewidth=2.4
    )
    for factor in ds.factor:
        series = ds.contributions.sel(factor=factor).to_series()
        ax.plot(
            ds.year,
            series.cumprod() / series.loc[reference_year],
            label=factor.item(),
            linestyle="--"
        )
    ax.set_ylabel(f"Change since {reference_year}")
    ax.set_title(country_id)
    ax.legend()

    fig.savefig(path_to_output)


def read_data(path_to_contributions, path_to_emissions):
    ds = xr.open_dataset(path_to_contributions)
    ds["emissions"] = (
        pd
        .read_csv(path_to_emissions, index_col=0)
        .unstack()
        .rename_axis(index=["country_id", "year"])
        .to_xarray()
    )
    return ds


if __name__ == "__main__":
    timeseries(
        path_to_contributions=snakemake.input.contributions,
        path_to_emissions=snakemake.input.emissions,
        country_id=snakemake.wildcards.country_id,
        from_year=int(snakemake.wildcards.from_year),
        to_year=int(snakemake.wildcards.to_year),
        path_to_output=snakemake.output[0]
    )
