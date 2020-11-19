"""Calculate the contribution of kaya factors to emissions.

This is based on the logarithmic mean division method as defined in

Wang, Q., Zhao, M., Li, R., &#38; Su, M. (2018). Decomposition and decoupling analysis of
carbon emissions from economic growth: A comparative study of China and the United States.
<i>Journal of Cleaner Production</i>, <i>197</i>, 178–184.
https://doi.org/10.1016/j.jclepro.2018.05.285
"""

import numpy as np
import pandas as pd


def contribution(path_to_emissions, path_to_kaya_factor, path_to_output):
    emissions = pd.read_csv(path_to_emissions, index_col=0)
    kaya_factor = pd.read_csv(path_to_kaya_factor, index_col=0)

    nominator = emissions - emissions.shift(1)
    denominator = emissions.applymap(np.log10) - emissions.shift(1).applymap(np.log10)
    factor = (kaya_factor / kaya_factor.shift(1)).applymap(np.log10)
    contribution = (nominator / denominator) * factor
    contribution.to_csv(path_to_output, index=True, header=True)


if __name__ == "__main__":
    contribution(
        path_to_emissions=snakemake.input.emissions,
        path_to_kaya_factor=snakemake.input.kaya_factor,
        path_to_output=snakemake.output[0]
    )
