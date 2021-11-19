from dataclasses import dataclass


@dataclass
class Crisis:
    name: str
    slug: str
    pre_from_year: int
    pre_to_year: int
    post_from_year: int
    post_to_year: int

    @classmethod
    def from_config(cls, slug, config):
        return Crisis(
            name=config["name"],
            slug=slug,
            pre_from_year=config["pre-from-year"],
            pre_to_year=config["first-year"] - 1,
            post_from_year=config["first-year"],
            post_to_year=config["post-to-year"]
        )
