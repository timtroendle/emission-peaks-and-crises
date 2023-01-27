import gettext

import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import pycountry

from crisis import Crisis


german = gettext.translation('iso3166', pycountry.LOCALES_DIR, languages=['de'])
german.install()

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
    "gdp": "BIP",
    "carbon-intensity": "Kohlenstoffintensität",
    "energy-intensity": "Energieintensität",
    "population": "Bevölkerung"
}

LINE_STYLES = {
    "gdp": "-",
    "carbon-intensity": "--",
    "energy-intensity": "-.",
    "population": ":"
}


def timeseries(path_to_contributions: str, path_to_emissions: str, country_id: str, crisis: Crisis,
               years_before_crisis_start: int, years_after_crisis_start: int, plot_crisis: bool,
               path_to_output: str):

    ds = (
        read_data(path_to_contributions, path_to_emissions)
    )
    fig = plt.figure(figsize=(7.09, 3)) # TODO double-check dimensions
    ax = fig.subplots(1, 1, sharex=False, squeeze=True)

    plot_timeseries(ds, country_id, crisis, years_before_crisis_start, years_after_crisis_start, plot_crisis, ax)

    ax.legend(
        framealpha=1.0,
        ncol=6,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.1),
        frameon=False
    )
    ax.set_ylabel("Veränderung seit\nReferenzjahr")
    sns.despine(fig, right=False)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.2)
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


def plot_timeseries(ds, country_id, crisis, years_before_crisis_start, years_after_crisis_start, plot_crisis, ax):
    period = crisis.national_period(country_id)
    ds_crisis = ds.sel(country_id=country_id, year=range(period.pre_from_year, period.post_to_year + 1)).dropna("year", how="any")
    emissions = ds_crisis.emissions
    ax.plot(
        ds_crisis.year,
        emissions / emissions.isel(year=0),
        label="$\mathrm{CO_2}$-Emissionen",
        linewidth=2.4
    )
    for factor in NICE_FACTOR_NAMES.keys():
        series = ds_crisis.contributions.sel(factor=factor).to_series()
        series.iloc[0] = 1 # to establish as reference year
        ax.plot(
            ds_crisis.year,
            series.cumprod(),
            label=NICE_FACTOR_NAMES[factor],
            linestyle=LINE_STYLES[factor]
        )
    if plot_crisis:
        ax.axvspan(
            xmin=period.from_year,
            xmax=period.to_year + 1,
            ymin=0,
            ymax=1,
            linewidth=0.0,
            alpha=0.2,
            color=GREY,
            label="Krise"
        )
    ax.set_title(_(country_name(country_id)))
    ax.set_xlim(period.from_year - years_before_crisis_start, period.from_year + years_after_crisis_start)
    ax.get_xaxis().set_major_locator(MultipleLocator(5))
    ax.get_xaxis().set_minor_locator(MultipleLocator(1))
    ax.get_yaxis().set_tick_params(top=True, direction='in')


def country_name(country_id):
    name = pycountry.countries.lookup(country_id).name
    if name == "Russian Federation":
        name = "Russia"
    elif name == "Korea, Republic of":
        name = "South Korea"
    return name


if __name__ == "__main__":
    timeseries(
        path_to_contributions=snakemake.input.contributions,
        path_to_emissions=snakemake.input.emissions,
        country_id=snakemake.wildcards.country_id,
        crisis=Crisis.from_config(snakemake.wildcards.crisis, snakemake.params.all_crises[snakemake.wildcards.crisis]),
        years_before_crisis_start=int(snakemake.params.years_before_crisis_start),
        years_after_crisis_start=int(snakemake.params.years_after_crisis_start),
        plot_crisis=bool(snakemake.params.plot_crisis),
        path_to_output=snakemake.output[0]
    )
