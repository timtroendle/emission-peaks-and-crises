import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import country_converter as coco


BLUE = '#4F6DB8'
RED = '#A01914'


def plot_peak_timeline(path_to_emissions, path_to_plot, crises, path_to_peak_csv):
    cc = coco.CountryConverter()
    emissions = pd.read_csv(path_to_emissions, index_col=0)
    emissions.drop(
        columns=list(filter(lambda x: x not in cc.data.ISO3.values, emissions.columns)),
        inplace=True
    )
    peak_years = (
        emissions
        .idxmax()
        .rename("peak_year")
        .sort_values()
    )
    peak_years.to_csv(path_to_peak_csv, index=True, header=True)
    rolling_emissions = emissions.rolling(window=5, center=False).mean().dropna(how='all')
    n_peaked_all, n_not_peaked_all = cumsum_peaked(rolling_emissions)

    fig = plt.figure(figsize=(8, 3.5))
    ax = fig.add_subplot()
    ax.bar(n_peaked_all.index, n_peaked_all.values, label="Countries with emission peaks", color=BLUE, alpha=1)
    ax.bar(n_not_peaked_all.index, n_not_peaked_all.values, label="Countries without emission peaks", color=RED, alpha=1)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 10)
    ax.vlines(
        x=list(crises.values()),
        ymin=ax.get_ylim()[0],
        ymax=ax.get_ylim()[1],
        linestyles='dotted',
        linewidths=0.75
    )
    for name, year in crises.items():
        ax.annotate(
            s=name,
            xy=(year + 0.5, ax.get_ylim()[1]),
            va="top"
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
    n_not_peaked = -(emissions.notna().sum(axis=1) - n_peaked)
    return n_peaked, n_not_peaked


if __name__ == "__main__":
    plot_peak_timeline(
        path_to_emissions=snakemake.input.emissions,
        path_to_plot=snakemake.output.plot,
        crises=dict(zip(snakemake.params.crises_names, snakemake.params.crises_years)),
        path_to_peak_csv=snakemake.output.csv
    )
