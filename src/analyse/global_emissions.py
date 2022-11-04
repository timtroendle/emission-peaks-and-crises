import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from crisis import Crisis

LINE_COLOR = "#E64B35"
GREY = "#7F7F7F"

INDEX_WORLD = "WLD"


def plot_global_emissions(path_to_emissions, crises, crises_names, path_to_plot):
    global_emissions = pd.read_csv(path_to_emissions, index_col=0).loc[:, INDEX_WORLD].div(1000)

    fig = plt.figure(figsize=(8, 2.5))
    ax = fig.add_subplot(111)

    ax.plot(global_emissions, color=LINE_COLOR)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 5)
    ax.set_xlim(1965, 2021)
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(1))
    ax.set_ylabel("Global CO₂ emissions (Gt)")

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
    for crisis in crises:
        ax.annotate(
            text=crises_names[crisis],
            xy=(crisis.global_period.from_year + 0.25, ax.get_ylim()[1] - 1),
            va="top"
        )

    sns.despine()
    fig.tight_layout()
    fig.savefig(path_to_plot, dpi=300)


if __name__ == "__main__":
    crises = [Crisis.from_config(crisis_slug, snakemake.params.all_crises[crisis_slug])
              for crisis_slug in snakemake.params.crises_slugs]
    plot_global_emissions(
        path_to_emissions=snakemake.input.emissions,
        crises=crises,
        crises_names={crisis: name for crisis, name in zip(crises, snakemake.params.crises_names)},
        path_to_plot=snakemake.output[0]
    )
