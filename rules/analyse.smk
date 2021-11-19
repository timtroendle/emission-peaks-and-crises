from dataclasses import dataclass

from src.analyse.crisis import Crisis

CRISES = {
    crisis_name: Crisis.from_config(config["crises"][crisis_name], config["countries"])
    for crisis_name in config["crises"]
}


rule decoupling_index:
    message: "Calculate the decoupling index."
    input:
        src = "src/analyse/decoupling.py",
        emissions = rules.emissions.output[0],
        gdp = rules.gdp.output[0]
    output: "build/decoupling-index.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/decoupling.py"


rule multiplicative_contributions:
    message: "Calculate multiplicative contributions to emissions."
    input:
        src = "src/analyse/mul_contributions.py",
        population = rules.population.output[0],
        energy_intensity = rules.energy_intensity.output[0],
        gdp = rules.gdp_per_capita.output[0],
        carbon_intensity = rules.carbon_intensity.output[0]
    output:
        nc = "build/multiplicative-contributions.nc",
        csv = "build/multiplicative-contributions.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/mul_contributions.py"


rule plot_peaker:
    message: "Plot timeline of peaks."
    input:
        src = "src/analyse/peaker.py",
        emissions = rules.emissions.output[0]
    params:
        crises_years = [1974, 1979, 1990, 2008],
        crises_names = ["First and", "Second oil crisis", "Soviet Union collapse", "Financial crisis"]
    output:
        plot = "build/peaker.png",
        csv = "build/peaker.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/peaker.py"


rule plot_contribution_timeseries_panelA:
    message: "Plot timeseries of contributions, panel A."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
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
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelB:
    message: "Plot timeseries of contributions, panel B."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
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
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelC:
    message: "Plot timeseries of contributions, panel C."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
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
    script: "../src/analyse/contribution_timeseries_panel.py"
