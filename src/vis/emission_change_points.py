from itertools import cycle

import pandas as pd
import matplotlib.pylab as plt
import ruptures as rpt
from ruptures.utils import pairwise

COLOR_CYCLE = ["#4286f4", "#f44174"]
ALPHA = 0.2 # transparency of the colored background


def emissions_change_point_plots(path_to_emissions, country_code, penalty, path_to_output):
    emissions = pd.read_csv(path_to_emissions, index_col=0, parse_dates=True)
    signal = (emissions.diff() / emissions).mul(100).loc[:, country_code].dropna()
    algo = rpt.Pelt(model="l1").fit(signal.values)
    change_points = algo.predict(pen=penalty)
    change_points = signal.iloc[[0] + change_points[:-1] + [change_points[-1] - 1]].index
    fig = plt.figure(figsize=(8, 3))
    ax = fig.add_subplot(111)
    ax.set_ylabel("Carbon emissions (Mt)")
    emissions.loc[:, country_code].reindex_like(signal).plot(legend=True)
    for (start, end), col in zip(pairwise(change_points), cycle(COLOR_CYCLE)):
        ax.axvspan(max(signal.index[0], start), end, facecolor=col, alpha=ALPHA)
    fig.savefig(path_to_output)


if __name__ == "__main__":
    emissions_change_point_plots(
        path_to_emissions=snakemake.input.emissions,
        country_code=snakemake.wildcards.country,
        penalty=snakemake.params.penalty,
        path_to_output=snakemake.output[0]
    )
