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
    "gdp": "GDP",
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity",
    "energy-and-carbon-intensity": "Technological change"
}


def plot_prepost_panel(path_to_prepost_contributions, country_crisis_tuples, path_to_plot):
    da = read_contribution_factors(path_to_prepost_contributions)

    n_rows = len(country_crisis_tuples)
    fig = plt.figure(figsize=(8, n_rows * 0.5))
    axes = fig.subplots(n_rows, 4, sharex=True, sharey=True)

    for col_id, factor in enumerate(FACTORS.keys()):
        for ax, (country, crisis) in zip(axes[:, col_id], country_crisis_tuples):
            data = da.sel(factor=factor, country_id=country, crisis=crisis.slug)
            plot_contribution_factor_as_arrows(ax, data.sel(period="pre").item(), data.sel(period="post").item())

    for ax, (country_id, crisis) in zip(axes[:, 0], country_crisis_tuples):
        ax.set_ylabel(
            pycountry.countries.lookup(country_id).name + "\n" + f"({crisis.pre_from_year}–{crisis.post_to_year})",
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


def read_contribution_factors(path_to_data):
    da = (
        xr
        .open_dataset(path_to_data)["contributions"]
        .sel(factor=list(FACTORS.keys()))
    )
    da = (da - 1) * 100
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


def plot_contribution_factor_as_arrows(ax, pre, post):
    color = INCREASE_COLOR if post > pre else DECREASE_COLOR
    ax.plot(
        [0, 1],
        [pre, post],
        marker='o',
        color=color,
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
        path_to_prepost_contributions=snakemake.input.contributions,
        country_crisis_tuples=[(country_id, all_crises[crisis_slug])
                               for (country_id, crisis_slug) in snakemake.params.country_crisis_tuples],
        path_to_plot=snakemake.output[0]
    )
