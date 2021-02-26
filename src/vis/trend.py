import xarray as xr
import matplotlib.pyplot as plt


def plot_trend_connecting_dot_plot(path_to_trend_data, crisis_name, path_to_plot):
    df = (
        xr.open_dataset(path_to_trend_data)
        ["trend"]
        .sel(variable="emissions")
        .to_dataframe()
        .mul(100)
        .drop(columns=["variable"])
        .unstack()
        .loc[:, "trend"]
        .drop(columns=["crisis"])
        .sort_values(by='post_crisis')
    )
    country_range = range(1, len(df.index) + 1)

    fig = plt.figure(figsize=(8, len(df.index) / 4))
    ax = fig.subplots()

    ax.hlines(y=country_range, xmin=df['pre_crisis'], xmax=df['post_crisis'], color='grey', alpha=0.4)
    ax.vlines(x=0, ymin=country_range[0], ymax=country_range[-1], color='grey', linestyle=":", alpha=0.4)
    ax.scatter(
        df['post_crisis'].where(df["post_crisis"] <= df["pre_crisis"]),
        country_range,
        color='green',
        alpha=0.4,
        label='decrease'
    )
    ax.scatter(
        df['post_crisis'].where(df["post_crisis"] > df["pre_crisis"]),
        country_range,
        color='red',
        alpha=0.4,
        label="increase"
    )

    ax.set_yticks(country_range)
    ax.set_yticklabels(df.index)
    ax.set_title(
        f"Emissions growth before and after {crisis_name}",
        ha='left',
        x=0,
        y=1.025
    )
    ax.text(
        s="Circles show post-crisis growth, the beginnings of the lines pre-crisis growth.",
        ha='left',
        va="bottom",
        x=0,
        y=1.01,
        transform=ax.transAxes
    )
    ax.set_xlabel('Emissions growth (%)')
    ax.set_ylabel('Country')
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    plot_trend_connecting_dot_plot(
        path_to_trend_data=snakemake.input.trend,
        crisis_name=snakemake.params.crisis_name,
        path_to_plot=snakemake.output[0]
    )
