import math

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import seaborn as sns
import pycountry


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

NICE_FACTOR_NAMES = {
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity",
    "gdp": "GDP",
    "population": "Population"
}


def timeseries(path_to_contributions, path_to_emissions, country_ids, from_year, to_year,
               path_to_peak_years, crises_years, path_to_output):
    nrows = math.ceil(len(country_ids) / 2)

    reference_years = pd.read_csv(path_to_peak_years, index_col=0)["peak_year"]
    ds = (
        read_data(path_to_contributions, path_to_emissions)
        .sel(country_id=country_ids, year=range(from_year, to_year + 1))
    )
    ds['emissions'] = ds.emissions

    fig = plt.figure(figsize=(8, nrows * 2))
    axes = fig.subplots(nrows, 2, sharex=True, sharey=True, squeeze=False)
    for ax, country_id in zip(axes.flatten(), country_ids):
        emissions = ds.emissions.sel(country_id=country_id)
        ax.plot(
            ds.year,
            emissions / emissions.sel(year=reference_years[country_id]),
            label="$\mathrm{CO_2}$ emissions",
            linewidth=2.4
        )
        for factor in ["gdp", "carbon-intensity", "energy-intensity", "population"]:
            series = ds.contributions.sel(factor=factor, country_id=country_id).to_series()
            ax.plot(
                ds.year,
                series.cumprod() / series.cumprod().loc[reference_years[country_id]],
                label=NICE_FACTOR_NAMES[factor],
                linestyle="--"
            )
        ax.set_title(pycountry.countries.lookup(country_id).name)

    y_min = min(ax.get_ylim()[0] for ax in axes.flatten())
    y_max = min(ax.get_ylim()[1] for ax in axes.flatten())
    for ax, _ in zip(axes.flatten(), country_ids):
        y_min
        ax.vlines(
            x=crises_years,
            ymin=y_min,
            ymax=y_max,
            linestyles='--',
            linewidths=0.5
        )

    axes[0, 0].legend(framealpha=1.0, ncol=2)
    for i, ax in enumerate(axes.flatten()):
        if i % 2 == 0:
            ax.set_ylabel(f"Change since peak")
    if len(country_ids) < len(axes.flatten()):
        # Last axis is empty.
        axes[-1, 1].set_axis_off()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1)
    sns.despine(fig)
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
        path_to_peak_years=snakemake.input.peak_years,
        country_ids=snakemake.params.country_ids,
        from_year=int(snakemake.params.from_year),
        to_year=int(snakemake.params.to_year),
        crises_years=snakemake.params.crises_years,
        path_to_output=snakemake.output[0]
    )
