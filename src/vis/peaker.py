import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import country_converter as coco


BLUE = '#4F6DB8'
RED = '#A01914'


def plot_peak_timeline(path_to_emissions, path_to_plot, path_to_peak_csv):
    cc = coco.CountryConverter()
    emissions = pd.read_csv(path_to_emissions, index_col=0)
    emissions.drop(
        columns=list(filter(lambda x: x not in cc.data.ISO3.values, emissions.columns)),
        inplace=True
    )
    (
        emissions
        .idxmax()
        .rename("peak_year")
        .sort_values()
        .to_csv(path_to_peak_csv, index=True, header=True)
    )
    rolling_emissions = emissions.rolling(window=5, center=False).mean().dropna(how='all')
    n_peaked_all, n_not_peaked_all = cumsum_peaked(rolling_emissions)
    n_peaked_eu15, n_not_peaked_eu15 = cumsum_peaked(rolling_emissions[cc.data.ISO3[cc.data.EU15.notna()].values])

    fig = plt.figure(figsize=(8, 3.5))
    ax = fig.add_subplot()
    ax.bar(n_peaked_all.index, n_peaked_all.values, label="All with emission peak", color=BLUE, alpha=0.4)
    ax.bar(n_peaked_eu15.index, n_peaked_eu15.values, label="EU15 with emission peak", color=BLUE)
    ax.bar(n_not_peaked_all.index, n_not_peaked_all.values, label="All without emission peak", color=RED, alpha=0.4)
    ax.bar(n_not_peaked_eu15.index, n_not_peaked_eu15.values, label="EU15 without emission peak", color=RED)
    ax.set_xlabel("Year")
    ax.set_ylabel("Number of countries")
    ax.legend()
    sns.despine()
    fig.tight_layout()
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
        path_to_peak_csv=snakemake.output.csv
    )
