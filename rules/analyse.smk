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


rule multiplicative_contributions:
    message: "Calculate multiplicative contributions to emissions."
    input:
        src = "src/analyse/mul_contributions.py",
        population = rules.population.output[0],
        energy_intensity = rules.energy_intensity.output[0],
        gdp = rules.gdp_per_capita.output[0],
        carbon_intensity = rules.carbon_intensity.output[0],
        energy_and_carbon_intensity = rules.energy_and_carbon_intensity.output[0]
    output:
        nc = "build/multiplicative-contributions.nc",
        csv = "build/multiplicative-contributions.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/mul_contributions.py"


rule prepost_multiplicative_contributions:
    message: "Calculate average contributions pre and post all crises."
    input:
        src = "src/analyse/prepost_contributions.py",
        contributions = rules.multiplicative_contributions.output[0]
    params:
        crises = config["crises"]
    output:
        nc = "build/prepost-contributions.nc",
        csv = "build/prepost-contributions.csv"
    conda: "../envs/default.yaml"
    script: "../src/analyse/prepost_contributions.py"


rule plot_prepost_contributions:
    message: "Plot contributions pre and post crises."
    input:
        src = "src/analyse/prepost_panel.py",
        contributions = rules.prepost_multiplicative_contributions.output[0]
    params:
        countries_in_crises = {
            "first-oil-crisis": ["BEL", "FRA", "DEU", "LUX", "GBR"],
            "financial-crisis": ["ESP", "GRC", "IRL", "ITA", "JPN", "NLD", "PRT", "SVN", "USA"],
            "second-oil-crisis": ["SWE"],
            "brazilian-crisis": ["BRA"],
            "argentinian-crisis": ["ARG"]
        }
    output: "build/prepost-contributions.png"
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


rule plot_contribution_timeseries_panelA:
    message: "Plot timeseries of contributions, panel A."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc
    params:
        country_ids_to_crises = {
            "BEL": "first-oil-crisis",
            "DEU": "first-oil-crisis",
            "GBR": "first-oil-crisis",
            "LUX": "first-oil-crisis",
            "IRL": "financial-crisis",
            "ITA": "financial-crisis",
            "JPN": "financial-crisis",
            "NLD": "financial-crisis",
            "PRT": "financial-crisis",
            "USA": "financial-crisis"
        },
        years = 20,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/panelA.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelB:
    message: "Plot timeseries of contributions, panel B."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        country_ids_to_crises = {
            "FRA": "first-oil-crisis",
            "SWE": "second-oil-crisis",
            "ESP": "financial-crisis",
            "ARG": "argentinian-crisis",
            "BRA": "brazilian-crisis"
        },
        years = 20,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/panelB.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"


rule plot_contribution_timeseries_panelC:
    message: "Plot timeseries of contributions, panel C."
    input:
        src = "src/analyse/contribution_timeseries_panel.py",
        emissions = "build/emissions-in-mt-bp.csv",
        contributions = rules.multiplicative_contributions.output.nc,
    params:
        country_ids_to_crises = {
            "GRC": "financial-crisis",
            "SVN": "financial-crisis",
        },
        years = 20,
        all_crises = config["crises"]
    output: "build/contribution-timeseries/panelC.png"
    conda: "../envs/default.yaml"
    script: "../src/analyse/contribution_timeseries_panel.py"
