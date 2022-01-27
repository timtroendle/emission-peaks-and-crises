import sys

import pytest
import pandas as pd
import xarray as xr


def run_test(path_to_test_dir: str, path_to_emissions: str, path_to_contribution_factors: str,
             path_to_output: str, countries: list):
    exit_code = pytest.main(
        [
            path_to_test_dir,
            f"--html={path_to_output}",
            "--self-contained-html",
            "--verbose"
        ],
        plugins=[
            _create_config_plugin(
                path_to_emissions=path_to_emissions,
                path_to_contribution_factors=path_to_contribution_factors,
                countries=countries
            )
        ]
    )
    sys.exit(exit_code)


def _create_config_plugin(path_to_emissions: str, path_to_contribution_factors: str, countries: list):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture()
        def emissions(self):
            return pd.read_csv(path_to_emissions, index_col=0)

        @pytest.fixture(scope="session")
        def contribution_factors(self):
            return xr.open_dataset(path_to_contribution_factors)["contributions"]

        @pytest.fixture(params=countries)
        def country_id(self, request):
            return request.param

        @pytest.fixture(scope="session", params=xr.open_dataset(path_to_contribution_factors).factor)
        def contribution_factor(self, request):
            return request.param

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(
        path_to_test_dir=snakemake.input.test_dir,
        path_to_emissions=snakemake.input.emissions,
        path_to_contribution_factors=snakemake.input.contribution_factors,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
