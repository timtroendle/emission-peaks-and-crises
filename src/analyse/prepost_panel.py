from itertools import chain

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pycountry
import xarray as xr

from crisis import Crisis


INCREASE_COLOR = "#E64B35"
DECREASE_COLOR = "#4DBBD5"
BAR_WIDTH = 0.8
ZERO_LINE_COLOR = "black"

FACTORS = {
    "gdp-and-population": "GDP",
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity",
    "energy-and-carbon-intensity": "Structural change"
}


def plot_prepost_panel(path_to_prepost_growth, crises_countries, path_to_plot):
    da = read_growth_rates(path_to_prepost_growth)

    n_rows = len(crises_countries)
    fig = plt.figure(figsize=(7.09, n_rows * 0.5 * 0.88625))
    axes = fig.subplots(n_rows, 4, sharex=True, sharey=True)

    for col_id, factor in enumerate(FACTORS.keys()):
        for ax, (crisis, country_id) in zip(axes[:, col_id], crises_countries):
            data = da.sel(factor=factor, country_id=country_id, crisis=crisis.slug)
            plot_growth_rates_as_arrows(ax, data.sel(period="pre").item(), data.sel(period="post").item())

    for ax, (crisis, country_id) in zip(axes[:, 0], crises_countries):
        period = crisis.national_period(country_id)
        country = country_name(country_id)
        ax.set_ylabel(
            f"{country}, {period.from_year}\n({period.pre_from_year}–{period.post_to_year})",
            rotation='horizontal',
            ha='right',
            va='center'
        )

    for ax, factor in zip(axes[0, :], FACTORS.values()):
        ax.set_title(factor)
    for ax, factor in zip(axes[-1, :], FACTORS.values()):
        ax.tick_params(bottom=False, labelbottom=True)
        ax.get_xaxis().set_ticklabels(["pre-crisis", "post-crisis"])

    sns.despine(fig, left=True, right=True, bottom=True, top=True)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0.1)
    fig.savefig(path_to_plot)


def country_name(country_id):
    name = pycountry.countries.lookup(country_id).name
    if name == "Russian Federation":
        name = "Russia"
    elif name == "Korea, Republic of":
        name = "South Korea"
    return name


def read_growth_rates(path_to_data):
    da = (
        xr
        .open_dataset(path_to_data)["growth_rate"]
        .sel(factor=list(FACTORS.keys()))
    )
    da = da * 100 # to percent
    return da


def plot_contribution_factor_as_bars(ax, pre, post):
    pre_color = INCREASE_COLOR if pre >= 0 else DECREASE_COLOR
    post_color = INCREASE_COLOR if post >= 0 else DECREASE_COLOR
    bars = ax.bar(
        [0, 1, 2],
        [pre, post],
        color=[pre_color, post_color],
        width=BAR_WIDTH
    )
    ax.bar_label(bars, labels=[f"{number:.1f}%" for number in [pre, post]])
    ax.hlines(
        y=0,
        xmin=0 - BAR_WIDTH / 2 * 1.5,
        xmax=1 + BAR_WIDTH / 2 * 1.5,
        color=ZERO_LINE_COLOR,
        linewidth=0.75
    )


def plot_growth_rates_as_arrows(ax, pre, post):
    line_color = INCREASE_COLOR if post > pre else DECREASE_COLOR
    pre_color = INCREASE_COLOR if pre > 0 else DECREASE_COLOR
    post_color = INCREASE_COLOR if post > 0 else DECREASE_COLOR
    ax.plot(
        [0, 1],
        [pre, post],
        color=line_color,
        linewidth=1.75
    )
    ax.plot(
        [0],
        [pre],
        marker="o",
        color=pre_color,
        linewidth=1.75
    )
    ax.plot(
        [1],
        [post],
        marker="o",
        color=post_color,
        linewidth=1.75
    )
    ax.hlines(
        y=0,
        xmin=0 - BAR_WIDTH / 2 * 0.5,
        xmax=1 + BAR_WIDTH / 2 * 1.75,
        color=ZERO_LINE_COLOR,
        linewidth=0.25
    )
    ax.annotate(
        text="$\mathrm{\Delta}$: " + f"{post - pre:.1f}%",
        xy=(1 + BAR_WIDTH / 2 * 1.75, 0.2),
        ha='right',
        va='bottom'
    )
    ax.get_yaxis().set_ticks([0])
    ax.get_xaxis().set_ticks([0, 1])
    ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=True)
    for tick in ax.get_yaxis().get_major_ticks():
        tick.set_pad(-8)


if __name__ == "__main__":
    all_crises = {crisis_slug: Crisis.from_config(crisis_slug, snakemake.params.all_crises[crisis_slug])
                  for crisis_slug in snakemake.params.all_crises}
    plot_prepost_panel(
        path_to_prepost_growth=snakemake.input.growth_rates,
        crises_countries=[
            (all_crises[crisis_slug], country_id)
            for crisis_slug in snakemake.params.crises_countries.keys()
            for country_id in snakemake.params.crises_countries[crisis_slug]
        ],
        path_to_plot=snakemake.output[0]
    )
