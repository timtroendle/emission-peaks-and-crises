import xarray as xr
import matplotlib.pyplot as plt


def decoupling_index(path_to_trend, crisis_name, path_to_plot):
    da = xr.open_dataset(path_to_trend)["trend"]
    index = da.sel(variable="emissions") / da.sel(variable="gdp")
    df = (
        index
        .sel(period=["pre_crisis", "post_crisis"])
        .to_series()
        .unstack()
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
        f"Decoupling index before and after {crisis_name}",
        ha='left',
        x=0,
        y=1.025
    )
    ax.text(
        s="Circles show post-crisis indices, the beginnings of the lines pre-crisis indices.",
        ha='left',
        va="bottom",
        x=0,
        y=1.01,
        transform=ax.transAxes
    )
    ax.set_xlabel(f'Decoupling index')
    ax.set_ylabel('Country')
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    decoupling_index(
        path_to_trend=snakemake.input.trend,
        crisis_name=snakemake.params.crisis_name,
        path_to_plot=snakemake.output[0]
    )
