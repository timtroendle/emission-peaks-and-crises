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
