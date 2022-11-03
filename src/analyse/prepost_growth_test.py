import typing

import pandas as pd
import xarray as xr
import matplotlib as mpl
from matplotlib import figure
import seaborn.objects as so
import scipy.stats


def test_prepost_growth(path_to_prepost_growth: str, paths_to_output: typing.Any):
    all_data_points, datas = read_data(path_to_prepost_growth)
    plot_qq(datas, paths_to_output.qq)
    plot_hist(all_data_points, paths_to_output.hist)
    perform_test(datas, paths_to_output.test)


def plot_qq(datas: dict[str, pd.Series], path_to_plot: str):
    fig = figure.Figure(figsize=(8, 4))
    axes = fig.subplots(2, 3).flatten()
    for ax, name, data in zip(axes, datas.keys(), datas.values()):
        scipy.stats.probplot(data, dist="norm", plot=ax)
        ax.set_title(name)
    fig.tight_layout()
    fig.savefig(path_to_plot)


def plot_hist(all_data_points: pd.DataFrame, path_to_plot: str):
    fig = figure.Figure(figsize=(8, 5))
    p = (
        so
        .Plot(all_data_points.replace({"gdp-and-population": "GDP", "energy-and-carbon-intensity": "Economic structure"}), x="diffs", color="factor")
        .layout(engine="tight")
        .facet(all_data_points.groupby(["factor", "peaked"]).ngroup(), wrap=3, order=[3, 4, 5, 0, 1, 2]) # necessary because of facet titles
        .theme(mpl.rc_params())
        .add(so.Bar(), so.Hist("count"))
        .on(fig)
        .label(
            x="Change in growth during crisis",
            y="Country-crisis (#)",
            color=str.capitalize,
            title="",
        )
    )
    p.show() # necessary to populate fig
    fig.get_axes()[0].set_title("No peak-and-decline")
    fig.get_axes()[1].set_title("Peak-and-decline")
    fig.get_axes()[2].set_title("Peaked before")
    fig.savefig(path_to_plot)


def perform_test(datas: dict[str, pd.Series], path_to_csv: str):
    names = datas.keys()
    values = datas.values()
    tests = [scipy.stats.wilcoxon(data, alternative="less") for data in values]
    test_result = pd.DataFrame(
        index=names,
        data={
            "n": [len(data) for data in values],
            "median": [data.median() for data in values],
            "q25": [data.quantile(0.25) for data in values],
            "q75": [data.quantile(0.75) for data in values],
            "statistic_V": [test.statistic for test in tests],
            "p_value": [test.pvalue for test in tests],
            "method": ["Wilcoxon signed rank test"],
            "alternative": ["less"]
        }
    )
    test_result.to_csv(path_to_csv, index=True, header=True)


def read_data(path_to_prepost_growth: str) -> tuple[pd.DataFrame, dict[str, pd.Series]]:
    ds = xr.open_dataset(path_to_prepost_growth)
    peak_and_decline = ds["peaked"] == "peak and decline"
    no_peak_and_decline = ds["peaked"] == "no peak and decline"
    all_types = (ds["peaked"] == "peak and decline") | (ds["peaked"] == "no peak and decline") | (ds["peaked"] == "peaked before")

    gdp_diff_peak = ds.sel(factor="gdp-and-population").where(peak_and_decline).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])
    gdp_diff_nopeak = ds.sel(factor="gdp-and-population").where(no_peak_and_decline).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])
    gdp_diff_all = ds.sel(factor="gdp-and-population").where(all_types).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])

    sc_diff_peak = ds.sel(factor="energy-and-carbon-intensity").where(peak_and_decline).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])
    sc_diff_nopeak = ds.sel(factor="energy-and-carbon-intensity").where(no_peak_and_decline).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])
    sc_diff_all = ds.sel(factor="energy-and-carbon-intensity").where(all_types).sel(period="pre", drop=True).to_dataframe().dropna(subset=["diffs"])

    all_data_points = pd.concat([
        gdp_diff_all,
        sc_diff_all]
    )
    datas = [gdp_diff_nopeak, gdp_diff_peak, gdp_diff_all, sc_diff_nopeak, sc_diff_peak, sc_diff_all]
    datas = [data["diffs"] for data in datas]
    names = ["GDP no peak", "GDP peak", "GDP all", "Economic structure no peak", "Economic structure peak", "Economic structure all"]
    return all_data_points, dict(zip(names, datas))


if __name__ == "__main__":
    test_prepost_growth(
        path_to_prepost_growth=snakemake.input.prepost,
        paths_to_output=snakemake.output
    )
