from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import seaborn as sns

from crisis import Crisis


BLUE = '#4F6DB8'
RED = '#A01914'
GREY = "#7F7F7F"

PEAK_YEAR_THRESHOLD = 2010 # We cannot identify peaks after 2010.


def plot_peak_timeline(path_to_emissions, paths_to_flags, crises, crises_names, high_income, middle_income,
                       path_to_plot, path_to_peak_csv):
    emissions = pd.read_csv(path_to_emissions, index_col=0)
    emissions = emissions[high_income + middle_income]
    rolling_emissions = emissions.rolling(window=5, center=False).mean().dropna(how='all')
    peak_years = (
        rolling_emissions
        .idxmax()
        .rename("peak_year")
        .sort_values()
    )
    peak_years.to_csv(path_to_peak_csv, index=True, header=True)
    n_peaked_all, n_not_peaked_all = cumsum_peaked(rolling_emissions)
    n_peaked_high, n_not_peaked_high = cumsum_peaked(rolling_emissions[high_income])

    fig = plt.figure(figsize=(8, 3.5))
    ax = fig.add_subplot()
    for crisis in crises:
        ax.axvspan(
            xmin=crisis.global_period.from_year,
            xmax=crisis.global_period.to_year + 1,
            ymin=0,
            ymax=1,
            linewidth=0.0,
            alpha=0.2,
            color=GREY
        )
    ax.bar(n_peaked_all.index, n_peaked_all.values, label="Medium-income with emission peak", color=BLUE, alpha=0.4)
    ax.bar(n_peaked_high.index, n_peaked_high.values, label="High-income with emission peak", color=BLUE)
    ax.bar(n_not_peaked_all.index, n_not_peaked_all.values, label="Medium-income without emission peak", color=RED, alpha=0.4)
    ax.bar(n_not_peaked_high.index, n_not_peaked_high.values, label="High-income without emission peak", color=RED)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 30)
    for crisis in crises:
        ax.annotate(
            text=crises_names[crisis],
            xy=(crisis.global_period.from_year + 0.25, ax.get_ylim()[1] - 1),
            va="top"
        )

    flags = {Path(path).stem: plt.imread(path)
             for path in paths_to_flags}
    for year, countries in peak_years.groupby(peak_years.values).groups.items():
        plot_flags(
            year=year,
            flags=[flags[country] for country in countries],
            ax=ax
        )

    ax.set_xlabel("Year")
    ax.set_ylabel("Number of countries")
    ax.legend(loc="lower right")
    sns.despine()
    fig.tight_layout()
    abs_ticklabels = [t for t in ax.get_yaxis().get_ticklabels()]
    for t in abs_ticklabels:
        t.set_text(int(abs(t.get_position()[1])))
    ax.get_yaxis().set_ticklabels(abs_ticklabels)
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(1))
    fig.savefig(path_to_plot, dpi=300)


def cumsum_peaked(emissions):
    peaks = emissions.idxmax()
    n_peaked = (
        peaks
        .value_counts()
        .sort_index()
        .cumsum()
        .reindex(emissions.index)
        .fillna(method="ffill")
        .fillna(0)
        [:-1] # all countries "peak" in the last year
    )
    n_peaked.rename("Number of countries with emission peaks", inplace=True)
    n_peaked.index.rename("Year", inplace=True)
    n_not_peaked = -(emissions.shape[1] - n_peaked)
    return n_peaked.loc[:PEAK_YEAR_THRESHOLD], n_not_peaked.loc[:PEAK_YEAR_THRESHOLD]


def plot_flags(year, flags, ax):
    if year < PEAK_YEAR_THRESHOLD:
        y = 60
        for flag in flags:
            plot_flag((year, y), flag, ax)
            y = y - 4


def plot_flag(coords, flag, ax):
    im = OffsetImage(flag, zoom=0.06)
    im.image.axes = ax

    ab = AnnotationBbox(im, coords,  xybox=(0., -16.), frameon=False,
                        xycoords='data',  boxcoords="offset points", pad=0)

    ax.add_artist(ab)


if __name__ == "__main__":
    crises = [Crisis.from_config(crisis_slug, snakemake.params.all_crises[crisis_slug])
              for crisis_slug in snakemake.params.crises_slugs]
    plot_peak_timeline(
        path_to_emissions=snakemake.input.emissions,
        paths_to_flags=snakemake.input.flags,
        path_to_plot=snakemake.output.plot,
        crises=crises,
        crises_names={crisis: name for crisis, name in zip(crises, snakemake.params.crises_names)},
        high_income=snakemake.params.high_income,
        middle_income=snakemake.params.middle_income,
        path_to_peak_csv=snakemake.output.csv
    )
