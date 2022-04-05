from dataclasses import dataclass


rule decoupling_index:
    message: "Calculate the decoupling index."
    input:
        src = "src/analyse/decoupling.py",
        emissions = "build/emissions-in-mt-bp.csv",
        gdp = rules.gdp.output[0]
    output: "build/decoupling-index.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/decoupling.py"


rule plot_global_emissions:
    message: "Plot global emissions."
    input:
        src = "src/analyse/global_emissions.py",
        emissions = "build/emissions-in-mt-bp.csv"
    params:
        crises_slugs = ["first-oil-crisis", "second-oil-crisis", "soviet-union-collapse", "financial-crisis"],
        all_crises = config["crises"],
        crises_names = ["First and", "Second oil crisis", "Soviet Union collapse", "Financial crisis"],
    output: "build/global-emissions.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/global_emissions.py"


rule multiplicative_contributions:
    message: "Calculate multiplicative contributions to emissions."
    input:
        src = "src/analyse/mul_contributions.py",
        population = rules.population.output[0],
        energy_intensity = rules.energy_intensity.output[0],
        gdp = rules.gdp_per_capita.output[0],
        carbon_intensity = rules.carbon_intensity.output[0],
        energy_and_carbon_intensity = rules.energy_and_carbon_intensity.output[0],
        gdp_and_population = rules.gdp.output[0]
    output:
        nc = "build/multiplicative-contributions.nc",
        csv = "build/multiplicative-contributions.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/mul_contributions.py"


rule prepost_growth_rates:
    message: "Calculate growth rates pre and post all crises."
    input:
        src = "src/analyse/prepost_growth.py",
        population = rules.population.output[0],
        energy_intensity = rules.energy_intensity.output[0],
        gdp = rules.gdp_per_capita.output[0],
        carbon_intensity = rules.carbon_intensity.output[0],
        energy_and_carbon_intensity = rules.energy_and_carbon_intensity.output[0],
        gdp_and_population = rules.gdp.output[0]
    params:
        crises = config["crises"]
    output:
        nc = "build/prepost-growth-rates.nc",
        csv = "build/prepost-growth-rates.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_growth.py"


rule plot_prepost_growth_peak_and_decline:
    message: "Plot growth rates pre and post crises in peak-and-decline countries."
    input:
        src = "src/analyse/prepost_panel.py",
        growth_rates = rules.prepost_growth_rates.output[0]
    params:
        crises_countries = config["highlights"]["peak-and-decline"],
        all_crises = config["crises"],
    output: "build/prepost-growth-peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_panel.py"


rule plot_prepost_growth_no_peak_and_decline:
    message: "Plot growth rates pre and post crises in no peak-and-decline countries."
    input:
        src = "src/analyse/prepost_panel.py",
        growth_rates = rules.prepost_growth_rates.output[0]
    params:
        crises_countries = config["highlights"]["no-peak-and-decline"],
        all_crises = config["crises"],
    output: "build/prepost-growth-no-peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_panel.py"


rule plot_peaker:
    message: "Plot timeline of peaks."
    input:
        src = "src/analyse/peaker.py",
        emissions = "build/emissions-in-mt-{source}.csv"
    params:
        crises_slugs = ["first-oil-crisis", "second-oil-crisis", "soviet-union-collapse", "financial-crisis"],
        all_crises = config["crises"],
        crises_names = ["First and", "Second oil crisis", "Soviet Union collapse", "Financial crisis"],
        high_income = ["AUS", "AUT", "BEL", "CAN", "CHL", "CZE", "DEU", "DNK", "ESP", "EST", "FIN",
                       "FRA", "GBR", "GRC", "HUN", "IRL", "ISL", "ISR", "ITA", "JPN", "KOR", "LTU",
                       "LUX", "LVA", "NLD", "NOR", "NZL", "POL", "PRT", "SAU", "SVK", "SVN", "SWE",
                       "CHE", "USA"],
        middle_income = ["ARG", "BRA", "CHN", "COL", "IDN", "IND", "MEX", "RUS", "TUR", "ZAF"]
    wildcard_constraints:
        source = "((bp)|(iea))"
    output:
        plot = "build/peaker-{source}.png",
        csv = "build/peaker-{source}.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/peaker.py"


rule plot_contribution_timeseries_peak_and_decline:
    message: "Plot timeseries of contributions, peak and decline."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc
    params:
        crises_countries = config["highlights"]["peak-and-decline"],
        years = 20,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_no_peak_and_decline:
    message: "Plot timeseries of contributions, no peak and decline."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        crises_countries = config["highlights"]["no-peak-and-decline"],
        years = 20,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/no-peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_all:
    message: "Plot timeseries of contributions, all."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        crises_countries = config["highlights"]["all"],
        years = 50,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/all.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"
