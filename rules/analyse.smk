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
    output: "build/contribution/{from_year}-{to_year}-{country_id}.png"
    conda: "../envs/default.yaml"
    script: "../src/vis/contribution.py"
