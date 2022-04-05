import math
from itertools import chain

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import pycountry

from crisis import Crisis

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
GREY = "#7F7F7F"

NICE_FACTOR_NAMES = {
    "gdp": "GDP",
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity",
    "population": "Population"
}

LINE_STYLES = {
    "gdp": "-",
    "carbon-intensity": "--",
    "energy-intensity": "-.",
    "population": ":"
}


def timeseries(path_to_contributions, path_to_emissions, crises_countries,
               years_before_crisis_start, years_after_crisis_start,
               all_crises, path_to_output):
    country_ids = list(chain(*[crises_countries[crisis] for crisis in crises_countries.keys()]))
    crises = {country: crisis for crisis, countries in crises_countries.items() for country in countries}
    nrows = math.ceil(len(country_ids) / 2)

    ds = (
        read_data(path_to_contributions, path_to_emissions)
        .sel(country_id=country_ids)
    )

    fig = plt.figure(figsize=(8, nrows * 2))
    axes = fig.subplots(nrows, 2, sharex=False, sharey=True, squeeze=False)
    for ax, country_id in zip(axes.flatten(), country_ids):
        crisis = all_crises[crises[country_id]]
        plot_timeseries(ds, country_id, crisis, years_before_crisis_start, years_after_crisis_start, ax)

    axes[0, 0].legend(framealpha=1.0, ncol=2)
    for i, ax in enumerate(axes.flatten()):
        if i % 2 == 0:
            ax.set_ylabel("Change since crisis")
    if len(country_ids) < len(axes.flatten()):
        # Last axis is empty.
        axes[-1, 1].set_axis_off()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1)
    sns.despine(fig, right=False)
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


def plot_timeseries(ds, country_id, crisis, years_before_crisis_start, years_after_crisis_start, ax):
    period = crisis.national_period(country_id)
    ds_crisis = ds.sel(year=range(period.pre_from_year, period.post_to_year + 1))
    reference_year = period.from_year
    emissions = ds_crisis.emissions.sel(country_id=country_id)
    ax.plot(
        ds_crisis.year,
        emissions / emissions.sel(year=reference_year),
        label="$\mathrm{CO_2}$ emissions",
        linewidth=2.4
    )
    for factor in NICE_FACTOR_NAMES.keys():
        series = ds_crisis.contributions.sel(factor=factor, country_id=country_id).to_series()
        ax.plot(
            ds_crisis.year,
            series.cumprod() / series.cumprod().loc[reference_year],
            label=NICE_FACTOR_NAMES[factor],
            linestyle=LINE_STYLES[factor]
        )
    ax.axvspan(
        xmin=period.from_year,
        xmax=period.to_year + 1,
        ymin=0,
        ymax=1,
        linewidth=0.0,
        alpha=0.2,
        color=GREY,
        label="Crisis"
    )
    ax.set_title(pycountry.countries.lookup(country_id).name)
    ax.set_xlim(period.from_year - years_before_crisis_start, period.from_year + years_after_crisis_start)
    ax.get_xaxis().set_major_locator(MultipleLocator(10))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.get_yaxis().set_tick_params(top=True, direction='in')


if __name__ == "__main__":
    timeseries(
        path_to_contributions=snakemake.input.contributions,
        path_to_emissions=snakemake.input.emissions,
        crises_countries=snakemake.params.crises_countries,
        years_before_crisis_start=int(snakemake.params.years_before_crisis_start),
        years_after_crisis_start=int(snakemake.params.years_after_crisis_start),
        all_crises={crisis_slug: Crisis.from_config(crisis_slug, snakemake.params.all_crises[crisis_slug])
                    for crisis_slug in snakemake.params.all_crises},
        path_to_output=snakemake.output[0]
    )
