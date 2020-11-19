from dataclasses import dataclass
import pycountry


@dataclass
class Crisis:
    pre_from_year: int
    pre_to_year: int
    from_year: int
    to_year: int
    post_from_year: int
    post_to_year: int
    ids_countries: list

    @classmethod
    def from_config(cls, config):
        return Crisis(
            pre_from_year=config["pre-from-year"],
            pre_to_year=config["from-year"] - 1,
            from_year=config["from-year"],
            to_year=config["to-year"],
            post_from_year=config["to-year"] + 1,
            post_to_year=config["post-to-year"],
            ids_countries=[pycountry.countries.lookup(country_name).alpha_3
                           for country_name in config["countries"]]
        )


CRISES = {
    crisis_name: Crisis.from_config(config["crises"][crisis_name])
    for crisis_name in config["crises"]
}


def contribution_factor_barplots(wildcards):
    crisis_name = wildcards["crisis"]
    crisis = CRISES[crisis_name]
    pre_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.pre_from_year}-{crisis.pre_to_year}-{country_id}.png"
        for country_id in crisis.ids_countries
    ]
    crisis_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.from_year}-{crisis.to_year}-{country_id}.png"
        for country_id in crisis.ids_countries
    ]
    post_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.post_from_year}-{crisis.post_to_year}-{country_id}.png"
        for country_id in crisis.ids_countries

    ]
    return pre_plots + crisis_plots + post_plots


rule crisis_analysis:
    message: "Analyse {wildcards.crisis}."
    input:
        contribution_factor_barplots
    output: touch("build/crises/{crisis}/analysis.done")


rule decoupling_index:
    message: "Calculate the decoupling index."
    input:
        src = "src/decoupling.py",
        emissions = rules.emissions.output[0],
        gdp = rules.gdp.output[0]
    output: "build/decoupling-index.csv"
    conda: "../envs/default.yaml"
    script: "../src/decoupling.py"


rule contribution:
    message: "Calculate the contribution of {wildcards.kaya_factor} to emissions."
    input:
        src = "src/contribution.py",
        emissions = rules.emissions.output[0],
        kaya_factor = "build/{kaya_factor}.csv"
    wildcard_constraints:
        kaya_factor = "((population)|(energy-intensity-in-ej-per-usd)|(gdp-in-usd-per-capita)|(carbon-intensity-in-mt-per-ej))",
    output: "build/contribution-{kaya_factor}.csv"
    conda: "../envs/default.yaml"
    script: "../src/contribution.py"


rule contributions:
    message: "Merge all contribution factors."
    input:
        src = "src/contributions.py",
        population = "build/contribution-population.csv",
        energy_intensity = "build/contribution-energy-intensity-in-ej-per-usd.csv",
        gdp = "build/contribution-gdp-in-usd-per-capita.csv",
        carbon_intensity = "build/contribution-carbon-intensity-in-mt-per-ej.csv"
    output: "build/contributions.nc"
    conda: "../envs/default.yaml"
    script: "../src/contributions.py"


rule plot_contribution_bar_chart:
    message: "Plot a bar chart visualising contribution factors for country {wildcards.country_id}."
    input:
        src = "src/vis/contribution.py",
        contributions = rules.contributions.output[0],
    output: "build/crises/{crisis}/contribution/{from_year}-{to_year}-{country_id}.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution.py"
