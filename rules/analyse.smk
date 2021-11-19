from dataclasses import dataclass

from src.crisis import Crisis

CRISES = {
    crisis_name: Crisis.from_config(config["crises"][crisis_name], config["countries"])
    for crisis_name in config["crises"]
}


def contribution_factor_barplots(wildcards):
    crisis_name = wildcards["crisis"]
    crisis = CRISES[crisis_name]
    pre_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.pre_from_year}-{crisis.pre_to_year}-{country_id}.png"
        for country_id in crisis.country_ids
    ]
    crisis_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.from_year}-{crisis.to_year}-{country_id}.png"
        for country_id in crisis.country_ids
    ]
    post_plots = [
        f"build/crises/{{crisis}}/contribution/{crisis.post_from_year}-{crisis.post_to_year}-{country_id}.png"
        for country_id in crisis.country_ids

    ]
    return pre_plots + crisis_plots + post_plots


def contribution_factor_timeseries(wildcards):
    crisis_name = wildcards["crisis"]
    crisis = CRISES[crisis_name]
    plots = [
        f"build/crises/{{crisis}}/contribution-timeseries/{crisis.pre_from_year}-{crisis.post_to_year}-{country_id}.png"
        for country_id in crisis.country_ids
    ]
    return plots


def trend_analysis(wildcards):
    crisis_name = wildcards["crisis"]
    crisis = CRISES[crisis_name]
    plots = [
        f"build/crises/{{crisis}}/trend-emissions.png",
        f"build/crises/{{crisis}}/trend-population.png",
        f"build/crises/{{crisis}}/trend-gdp.png",
        f"build/crises/{{crisis}}/trend-carbon-intensity.png",
        f"build/crises/{{crisis}}/trend-energy-intensity.png",
        f"build/crises/{{crisis}}/decoupling-index.png",
        f"build/crises/{{crisis}}/trend-change-in-percentage-points.csv",
        f"build/crises/{{crisis}}/r-squared.csv"
    ]
    return plots


rule crisis_analysis:
    message: "Analyse {wildcards.crisis}."
    input:
        "build/crises/{crisis}/overview.csv",
        trend_analysis
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


rule multiplicative_contributions:
    message: "Calculate multiplicative contributions to emissions."
    input:
        src = "src/mul_contributions.py",
        population = rules.population.output[0],
        energy_intensity = rules.energy_intensity.output[0],
        gdp = rules.gdp_per_capita.output[0],
        carbon_intensity = rules.carbon_intensity.output[0]
    output:
        nc = "build/multiplicative-contributions.nc",
        csv = "build/multiplicative-contributions.csv"
    conda: "../envs/default.yaml"
    script: "../src/mul_contributions.py"


rule plot_contribution_bar_chart:
    message:
        "Plot a bar chart visualising contribution factors for country {wildcards.country_id} " +
        "during {wildcards.crisis}."
    input:
        src = "src/vis/contribution.py",
        contributions = rules.contributions.output[0],
    output: "build/crises/{crisis}/contribution/{from_year}-{to_year}-{country_id}.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution.py"


rule plot_emission_change_points:
    message: "Plot change points of emission time series for country {wildcards.country}."
    input:
        src = "src/vis/emission_change_points.py",
        emissions = rules.emissions.output[0]
    params: penalty = 3 # penalty value used in the PELT method
    output: "build/change-points/emissions-{country}.png"
    conda: "../envs/changepoint.yaml"
    script: "../src/vis/emission_change_points.py"


rule plot_contribution_timeseries:
    message: "Plot timeseries of contributions for country {wildcards.country_id} during {wildcards.crisis}."
    input:
        src = "src/vis/contribution_timeseries.py",
        emissions = rules.emissions.output[0],
        contributions = rules.multiplicative_contributions.output.nc
    output: "build/crises/{crisis}/contribution-timeseries/{from_year}-{to_year}-{country_id}.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution_timeseries.py"


rule plot_peaker:
    message: "Plot timeline of peaks."
    input:
        src = "src/vis/peaker.py",
        emissions = rules.emissions.output[0]
    params:
        crises_years = [1974, 1979, 1990, 2008],
        crises_names = ["First and", "Second oil crisis", "Soviet Union collapse", "Financial crisis"]
    output:
        plot = "build/peaker.png",
        csv = "build/peaker.csv"
    conda: "../envs/default.yaml"
    script: "../src/vis/peaker.py"


rule plot_contribution_timeseries_panelA:
    message: "Plot timeseries of contributions, panel A."
    input:
        src = "src/vis/contribution_timeseries.py",
        emissions = rules.emissions.output[0],
        contributions = rules.multiplicative_contributions.output.nc,
        peak_years = rules.plot_peaker.output.csv
    params:
        country_ids = ["BEL", "DEU", "GBR", "LUX", "IRL", "ITA", "JPN", "NLD", "PRT", "USA"],
        from_year = 1970,
        to_year = 2019,
        crises_years = [1974, 1979, 1990, 2008]
    output: "build/contribution-timeseries/panelA.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelB:
    message: "Plot timeseries of contributions, panel B."
    input:
        src = "src/vis/contribution_timeseries.py",
        emissions = rules.emissions.output[0],
        contributions = rules.multiplicative_contributions.output.nc,
        peak_years = rules.plot_peaker.output.csv
    params:
        country_ids = ["FRA", "SWE", "ESP", "ARG", "BRA"],
        from_year = 1970,
        to_year = 2019,
        crises_years = [1974, 1979, 1990, 2008]
    output: "build/contribution-timeseries/panelB.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelC:
    message: "Plot timeseries of contributions, panel C."
    input:
        src = "src/vis/contribution_timeseries.py",
        emissions = rules.emissions.output[0],
        contributions = rules.multiplicative_contributions.output.nc,
        peak_years = rules.plot_peaker.output.csv
    params:
        country_ids = ["GRC", "SVN"],
        from_year = 1970,
        to_year = 2019,
        crises_years = [1974, 1979, 1990, 2008]
    output: "build/contribution-timeseries/panelC.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution_timeseries_panel.py"

rule overview:
    message: "Create overview table for {wildcards.crisis}."
    input:
        src = "src/overview.py",
        emissions = rules.emissions.output[0],
        gdp = rules.gdp_per_capita.output[0]
    params:
        crisis = lambda wildcards: config["crises"][wildcards["crisis"]],
        countries = config["countries"]
    output: "build/crises/{crisis}/overview.csv"
    conda: "../envs/default.yaml"
    script: "../src/overview.py"


rule trend:
    message: "Determine trends in time series for {wildcards.crisis}."
    input:
        src = "src/trend.py",
        emissions = rules.emissions.output[0],
        gdp = rules.gdp.output[0],
        population = rules.population.output[0],
        carbon_intensity = rules.carbon_intensity.output[0],
        energy_intensity = rules.energy_intensity.output[0]
    params:
        crisis = lambda wildcards: config["crises"][wildcards["crisis"]],
        countries = config["countries"]
    output: "build/crises/{crisis}/trend.nc"
    conda: "../envs/default.yaml"
    script: "../src/trend.py"


rule results_as_csv:
    message: "Transform results into csv format for {wildcards.crisis}."
    input:
        src = "src/csv_results.py",
        trend = rules.trend.output[0]
    output:
        trend = "build/crises/{crisis}/trend.csv",
        trend_change = "build/crises/{crisis}/trend-change-in-percentage-points.csv",
        r_squared = "build/crises/{crisis}/r-squared.csv",
        p_value = "build/crises/{crisis}/worst-p-value.csv"
    conda: "../envs/default.yaml"
    script: "../src/csv_results.py"


rule plot_trend:
    message: "Plot pre- and post-crisis trends for {wildcards.variable} during {wildcards.crisis}."
    input:
        src = "src/vis/trend.py",
        trend = rules.trend.output[0]
    params:
        crisis_name = lambda wildcards: CRISES[wildcards["crisis"]].name,
        r2_threshold = config["report"]["r2_threshold"],
        p_value_threshold = config["report"]["p_value_threshold"]
    output: "build/crises/{crisis}/trend-{variable}.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/trend.py"


rule plot_decoupling_index:
    message: "Plot pre- and post-crisis decoupling_index during {wildcards.crisis}."
    input:
        src = "src/vis/decoupling_index.py",
        trend = rules.trend.output[0]
    params: crisis_name = lambda wildcards: CRISES[wildcards["crisis"]].name
    output: "build/crises/{crisis}/decoupling-index.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/decoupling_index.py"


rule method_comparison:
    message: "Apply different timeseries analysis methods for comparison -- {wildcards.crisis}."
    input:
        src = "src/methods.py",
        emissions = rules.emissions.output[0],
        gdp = rules.gdp.output[0],
        population = rules.population.output[0],
        carbon_intensity = rules.carbon_intensity.output[0],
        energy_intensity = rules.energy_intensity.output[0]
    params:
        crisis = lambda wildcards: config["crises"][wildcards["crisis"]],
        countries = config["countries"]
    output: "build/crises/{crisis}/methods.nc"
    conda: "../envs/default.yaml"
    script: "../src/methods.py"
