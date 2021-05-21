from dataclasses import dataclass

import pycountry


@dataclass
class Crisis:
    name: str
    pre_from_year: int
    pre_to_year: int
    from_year: int
    to_year: int
    post_from_year: int
    post_to_year: int
    country_ids: list

    @classmethod
    def from_config(cls, config, country_names):
        return Crisis(
            name=config["name"],
            pre_from_year=config["pre-from-year"],
            pre_to_year=config["from-year"] - 1,
            from_year=config["from-year"],
            to_year=config["to-year"],
            post_from_year=config["to-year"] + 1,
            post_to_year=config["post-to-year"],
            country_ids=[country_to_country_code(country_name)
                         for country_name in country_names]
        )


def country_to_country_code(country_name):
    if country_name == "World":
        return "WLD"
    elif country_name == "OECD":
        return "OED"
    elif country_name == "Non-OECD":
        return "NOE"
    else:
        return pycountry.countries.lookup(country_name).alpha_3
