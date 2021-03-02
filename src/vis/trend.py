import xarray as xr
import matplotlib.pyplot as plt

VARIABLE_NAMES = {
    "emissions": "Emissions",
    "gdp": "GDP",
    "population": "Population",
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity"
}


def plot_trend_connecting_dot_plot(path_to_trend_data, crisis_name, variable, r2_threshold, path_to_plot):
    ds = xr.open_dataset(path_to_trend_data)
    trend = (
        ds
        ["trend"]
        .sel(variable=variable)
        .to_dataframe()
        .mul(100)
        .drop(columns=["variable"])
        .unstack()
        .loc[:, "trend"]
        .drop(columns=["crisis"])
        .sort_values(by='post_crisis')
    )
    r2 = (
        ds
        ["r_squared"]
        .sel(variable=variable)
        .to_series()
    )
    country_range = range(1, len(trend.index) + 1)

    fig = plt.figure(figsize=(8, len(trend.index) / 4))
    ax = fig.subplots()

    ax.hlines(y=country_range, xmin=trend['pre_crisis'], xmax=trend['post_crisis'], color='grey', alpha=0.4)
    ax.vlines(x=0, ymin=country_range[0], ymax=country_range[-1], color='grey', linestyle=":", alpha=0.4)
    ax.scatter(
        trend['post_crisis'].where(trend["post_crisis"] <= trend["pre_crisis"]),
        country_range,
        color='green',
        alpha=0.4,
        label='decrease'
    )
    ax.scatter(
        trend['post_crisis'].where(trend["post_crisis"] > trend["pre_crisis"]),
        country_range,
        color='red',
        alpha=0.4,
        label="increase"
    )
    ax.scatter(
        trend['post_crisis'].where(r2 < r2_threshold),
        country_range,
        linewidth=1,
        color="k",
        marker="x"
    )

    ax.set_yticks(country_range)
    ax.set_yticklabels(trend.index)
    ax.set_title(
        f"{VARIABLE_NAMES[variable]} growth before and after {crisis_name}",
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
    ax.text(
        s=f"x: $r^2$ < {r2_threshold}",
        transform=ax.transAxes,
        ha="right",
        x=0.975,
        y=0.025
    )
    ax.set_xlabel(f'{VARIABLE_NAMES[variable]} growth (%)')
    ax.set_ylabel('Country')
    fig.savefig(path_to_plot)


if __name__ == "__main__":
    plot_trend_connecting_dot_plot(
        path_to_trend_data=snakemake.input.trend,
        crisis_name=snakemake.params.crisis_name,
        variable=snakemake.wildcards.variable,
        r2_threshold=snakemake.params.r2_threshold,
        path_to_plot=snakemake.output[0]
    )
