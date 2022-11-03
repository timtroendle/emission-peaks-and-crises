from dataclasses import dataclass


rule decoupling_index:
    message: "Calculate the decoupling index."
    input:
        emissions = "build/emissions-in-mt-bp.csv",
        gdp = rules.gdp.output[0]
    output: "build/decoupling-index.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/decoupling.py"


rule plot_global_emissions:
    message: "Plot global emissions."
    input:
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


rule prepost_growth_rates_peaks_marked:
    message: "Mark prepost growth rates with emission peaks."
    input:
        prepost = rules.prepost_growth_rates.output.nc
    params:
        peak_and_decline = config["highlights"]["peaks"]["peak-and-decline"],
        no_peak_and_decline = config["highlights"]["peaks"]["no-peak-and-decline"],
        peaked_before = config["highlights"]["peaks"]["peaked-before"]
    output:
        nc = "build/prepost-growth-rates-marked.nc",
        csv = "build/prepost-growth-rates-marked.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_growth_marked.py"


rule test_prepost_growth:
    message: "Run tests on the growth rates during crises."
    input:
        prepost = rules.prepost_growth_rates_peaks_marked.output[0]
    output:
        qq = "build/prepost-growth-qq-plots.png",
        hist = "build/prepost-growth-histograms.png",
        test = "build/prepost-growth-tests.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_growth_test.py"


rule plot_prepost_growth_peak_and_decline:
    message: "Plot growth rates pre and post crises in peak-and-decline countries."
    input:
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
        growth_rates = rules.prepost_growth_rates.output[0]
    params:
        crises_countries = config["highlights"]["no-peak-and-decline-prepost"],
        all_crises = config["crises"],
    output: "build/prepost-growth-no-peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_panel.py"


rule plot_peaker:
    message: "Plot timeline of peaks."
    input:
        emissions = "build/emissions-in-mt-{source}.csv",
        flags = expand(
            "build/flag/{country}.png",
            country=config["groups"]["high-income"] + config["groups"]["middle-income"]
        )
    params:
        crises_slugs = ["first-oil-crisis", "second-oil-crisis", "soviet-union-collapse", "financial-crisis"],
        all_crises = config["crises"],
        crises_names = ["First and", "Second oil crisis", "Soviet Union collapse", "Financial crisis"],
        high_income = config["groups"]["high-income"],
        middle_income = config["groups"]["middle-income"]
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
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc
    params:
        crises_countries = config["highlights"]["peak-and-decline"],
        years_before_crisis_start = 10,
        years_after_crisis_start = 18,
        all_crises = config["crises"],
        plot_crisis = True,
        share_y_axis = False
    output: "build/contribution-timeseries/peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_no_peak_and_decline:
    message: "Plot timeseries of contributions, no peak and decline."
    input:
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        crises_countries = config["highlights"]["no-peak-and-decline-timeline"],
        years_before_crisis_start = 10,
        years_after_crisis_start = 9,
        all_crises = config["crises"],
        plot_crisis = True,
        share_y_axis = False
    output: "build/contribution-timeseries/no-peak-and-decline.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_all:
    message: "Plot timeseries of contributions, all."
    input:
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        crises_countries = config["highlights"]["all"],
        years_before_crisis_start = 32,
        years_after_crisis_start = 18,
        all_crises = config["crises"],
        plot_crisis = False,
        share_y_axis = False
    output: "build/contribution-timeseries/all.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_non_crisis_peaker:
    message: "Plot timeseries of contributions, non crisis peaker."
    input:
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc
    params:
        crises_countries = config["highlights"]["non-crisis-peaker"],
        years_before_crisis_start = 10,
        years_after_crisis_start = 10,
        all_crises = config["crises"],
        plot_crisis = False,
        share_y_axis = False
    output: "build/contribution-timeseries/non-crisis-peaker.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"
