import seaborn as sns
import xarray as xr
import matplotlib.pyplot as plt

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


def contribution(path_to_contributions, from_year, to_year, country_id, path_to_output):
    assert from_year <= to_year
    contributions = (
        xr
        .open_dataset(path_to_contributions)["contributions"]
        .sel(country_id=country_id, year=[from_year, to_year])
        .sum("year")
        .to_series()
    )
    contributions["emissions"] = contributions.sum()

    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    sns.barplot(
        data=(
            contributions
            .reindex(["emissions", "gdp", "population", "carbon-intensity", "energy-intensity"])
            .reset_index()
        ),
        x="factor",
        y="contributions",
        palette=[COLOR_PALETTE[0]] + [COLOR_PALETTE[1]] * 4
    )
    sns.despine(bottom=True)
    ax.hlines(0, xmin=-1, xmax=4.5)
    ax.set_xlim(-0.5, 4.5)
    ax.set_xlabel("")
    ax.set_ylabel("Change in carbon emissions (Mt)")
    fig.savefig(path_to_output)


if __name__ == "__main__":
    contribution(
        path_to_contributions=snakemake.input.contributions,
        from_year=int(snakemake.wildcards.from_year),
        to_year=int(snakemake.wildcards.to_year),
        country_id=snakemake.wildcards.country_id,
        path_to_output=snakemake.output[0]
    )
