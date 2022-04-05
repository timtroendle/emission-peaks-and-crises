from dataclasses import dataclass, field
from typing import Mapping
import warnings


@dataclass(eq=True, frozen=True)
class CrisisPeriod:
    pre_from_year: int
    pre_to_year: int
    from_year: int
    to_year: int
    post_from_year: int
    post_to_year: int

    @classmethod
    def from_config(cls, config):
        return CrisisPeriod(
            pre_from_year=config["pre-from-year"],
            pre_to_year=config["first-year"] - 1,
            from_year=config["first-year"],
            to_year=config["final-year"],
            post_from_year=config["first-year"],
            post_to_year=config["post-to-year"]
        )


class hashabledict(dict):
  def __key(self):
    return tuple((k,self[k]) for k in sorted(self))
  def __hash__(self):
    return hash(self.__key())
  def __eq__(self, other):
    return self.__key() == other.__key()


@dataclass(eq=True, frozen=True)
class Crisis:
    name: str
    slug: str = field(repr=False)
    global_period: CrisisPeriod = field(repr=False)
    national_periods: Mapping[str, CrisisPeriod] = field(repr=False)

    @classmethod
    def from_config(cls, slug, config):
        if "national" in config:
            national_periods = hashabledict({
                country_id: CrisisPeriod.from_config(period)
                for country_id, period in config["national"].items()
            })
        else:
            national_periods = hashabledict({})
        return Crisis(
            name=config["name"],
            slug=slug,
            global_period=CrisisPeriod.from_config(config["global"]),
            national_periods=national_periods
        )


    def national_period(self, country_id):
        try:
            return self.national_periods[country_id]
        except KeyError: # take global period by default
            return self.global_period
