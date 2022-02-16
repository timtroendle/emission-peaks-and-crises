import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


NATURE_PALETTE = [
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
sns.set_palette(NATURE_PALETTE)

INDEX_WORLD = "WLD"

def plot_global_emissions(path_to_emissions, path_to_plot):
    global_emissions = pd.read_csv(path_to_emissions, index_col=0).loc[:, INDEX_WORLD].div(1000)

    fig = plt.figure(figsize=(8, 3.5))
    ax = fig.add_subplot(111)

    ax.plot(global_emissions)

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(1))
    ax.set_ylabel("Global CO₂ emissions (Gt)")
    ax.annotate(
        text="Financial crisis",
        xy=(2008, global_emissions.loc[2008]),
        xytext=(2008, 35),
        textcoords="data",
        horizontalalignment="right",
        arrowprops={"arrowstyle": "->"}
    )
    ax.annotate(
        text="Fall of the Soviet Union",
        xy=(1991, global_emissions.loc[1990]),
        xytext=(1991, 15),
        horizontalalignment="center",
        arrowprops={"arrowstyle": "->"}
    )
    ax.annotate(
        text="                   Oil crises",
        xy=(1973, global_emissions.loc[1973]),
        xytext=(1973, 25),
        horizontalalignment="center",
        arrowprops={"arrowstyle": "->"}
    )
    ax.annotate(
        text="      ",
        xy=(1979, global_emissions.loc[1979]),
        xytext=(1979, 25),
        horizontalalignment="center",
        arrowprops={"arrowstyle": "->"}
    )

    sns.despine()
    fig.savefig(path_to_plot, dpi=300)


if __name__ == "__main__":
    plot_global_emissions(
        path_to_emissions=snakemake.input.emissions,
        path_to_plot=snakemake.output[0]
    )
