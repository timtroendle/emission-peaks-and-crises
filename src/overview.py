import pandas as pd

from crisis import Crisis


def crisis_analysis(path_to_carbon_emissions, path_to_gdp, crisis, path_to_output):
    emissions = pd.read_csv(path_to_carbon_emissions, index_col=0).loc[:, crisis.country_ids]
    gdp = pd.read_csv(path_to_gdp, index_col=0).loc[:, crisis.country_ids]
    periods = [
        (
            slice(crisis.pre_from_year - 1, crisis.pre_to_year),
            f"1-pre-crisis-{crisis.pre_from_year}-{crisis.pre_to_year}"
        ),
        (
            slice(crisis.from_year - 1, crisis.to_year),
            f"2-crisis-{crisis.from_year}-{crisis.to_year}"
        ),
        (
            slice(crisis.post_from_year - 1, crisis.post_to_year),
            f"3-post-crisis-{crisis.post_from_year}-{crisis.post_to_year}"
        )
    ]
    period_data = [
        analyse_period(emissions.loc[period[0], :], gdp.loc[period[0], :], period[1])
        for period in periods
    ]
    (
        pd
        .concat(period_data)
        .pivot(columns="period")
        .to_csv(path_to_output, index=True, header=True)
    )


def analyse_period(emissions, gdp, name):
    emissions_0 = emissions.iloc[0]
    emissions_1 = emissions.iloc[-1]
    gdp_0 = gdp.iloc[0]
    gdp_1 = gdp.iloc[-1]
    delta_emissions = (emissions_1 - emissions_0) / emissions_0
    delta_gdp = (gdp_1 - gdp_0) / gdp_0
    decoupling_index = delta_emissions / delta_gdp
    return pd.DataFrame({
        "carbon_emissions": delta_emissions,
        "gdp": delta_gdp,
        "decoupling_index": decoupling_index,
        "coupling_type": coupling_type(delta_emissions, delta_gdp, decoupling_index),
        "period": name
    })


def coupling_type(delta_emissions, delta_gdp, decoupling_index):
    coupling_type = pd.Series(dtype="object", index=decoupling_index.index)
    rising_gdp = delta_gdp > 0
    falling_gdp = delta_gdp < 0
    rising_emissions = delta_emissions > 0
    falling_emissions = delta_emissions < 0
    coupled = (decoupling_index > 0.8) & (decoupling_index < 1.2)
    decoupled_small = decoupling_index < 0.8
    decoupled_large = decoupling_index > 1.2

    coupling_type[rising_gdp & falling_emissions & decoupled_small] = "strong_decoupling"
    coupling_type[rising_gdp & rising_emissions & decoupled_small] = "weak_decoupling"
    coupling_type[falling_gdp & falling_emissions & decoupled_large] = "recessive_decoupling"
    coupling_type[rising_gdp & rising_emissions & decoupled_large] = "expansive_negative_decoupling"
    coupling_type[falling_gdp & rising_emissions & decoupled_small] = "strong_negative_decoupling"
    coupling_type[falling_gdp & falling_emissions & decoupled_small] = "weak_negative_decoupling"
    coupling_type[rising_gdp & rising_emissions & coupled] = "expansive_coupling"
    coupling_type[falling_gdp & falling_emissions & coupled] = "recessive_coupling"
    return coupling_type


if __name__ == "__main__":
    crisis_analysis(
        path_to_carbon_emissions=snakemake.input.emissions,
        path_to_gdp=snakemake.input.gdp,
        crisis=Crisis.from_config(snakemake.params.crisis),
        path_to_output=snakemake.output[0]
    )
