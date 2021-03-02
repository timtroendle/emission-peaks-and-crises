import xarray as xr
import matplotlib.pyplot as plt

VARIABLE_NAMES = {
    "emissions": "Emissions",
    "gdp": "GDP",
    "population": "Population",
    "carbon-intensity": "Carbon intensity",
    "energy-intensity": "Energy intensity"
}
ALPHA_NOT_SIGNIFICANT = 0.2
ALPHA_SIGNIFICANT = 0.4


def plot_trend_connecting_dot_plot(path_to_trend_data, crisis_name, variable,
                                   r2_threshold, p_value_threshold, path_to_plot):
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
        .reindex_like(trend)
    )
    p_value = (
        ds
        ["p_value"]
        .sel(variable="emissions", period=["pre_crisis", "post_crisis"])
        .to_dataframe()
        .drop(columns=["variable"])
        .unstack()
        .loc[:, "p_value"]
        .reindex_like(trend)
    )
    pre_significant = p_value["pre_crisis"] <= p_value_threshold
    post_significant = p_value["post_crisis"] <= p_value_threshold

    fig = plt.figure(figsize=(8, len(trend.index) / 4))
    ax = fig.subplots()

    ax.scatter(
        trend['post_crisis'].where(trend["post_crisis"] <= trend["pre_crisis"]),
        trend.index,
        color='green',
        alpha=ALPHA_NOT_SIGNIFICANT,
        label='decrease'
    )
    ax.scatter(
        trend['post_crisis'].where((trend["post_crisis"] <= trend["pre_crisis"]) & post_significant),
        trend.index,
        color='green',
        alpha=ALPHA_SIGNIFICANT,
        label='decrease'
    )
    ax.scatter(
        trend['post_crisis'].where(trend["post_crisis"] > trend["pre_crisis"]),
        trend.index,
        color='red',
        alpha=ALPHA_NOT_SIGNIFICANT,
        label="increase"
    )
    ax.scatter(
        trend['post_crisis'].where((trend["post_crisis"] > trend["pre_crisis"]) & post_significant),
        trend.index,
        color='red',
        alpha=ALPHA_NOT_SIGNIFICANT,
        label="increase"
    )
    ax.scatter(
        trend['post_crisis'].where(r2 < r2_threshold),
        trend.index,
        linewidth=1,
        color="k",
        marker="x"
    )
    ax.scatter(
        trend['pre_crisis'].where(r2 < r2_threshold),
        trend.index,
        linewidth=1,
        color="k",
        marker="x"
    )

    ax.hlines(
        y=trend.index,
        xmin=trend['pre_crisis'],
        xmax=trend['post_crisis'],
        color='grey',
        zorder=0,
        alpha=ALPHA_NOT_SIGNIFICANT
    )
    ax.hlines(
        y=trend.index.where(pre_significant).dropna(),
        xmin=trend['pre_crisis'].where(pre_significant).dropna(),
        xmax=trend['post_crisis'].where(pre_significant).dropna(),
        color='grey',
        zorder=0,
        alpha=ALPHA_SIGNIFICANT
    )
    ax.vlines(
        x=0,
        ymin=trend.index[0],
        ymax=trend.index[-1],
        color='grey',
        linestyle=":",
        alpha=ALPHA_SIGNIFICANT
    )

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
        s=f"x: $r^2$ < {r2_threshold}\nSaturated values: p <= {p_value_threshold}",
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
        p_value_threshold=snakemake.params.p_value_threshold,
        path_to_plot=snakemake.output[0]
    )
