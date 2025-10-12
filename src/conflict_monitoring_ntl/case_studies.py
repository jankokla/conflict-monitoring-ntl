import datetime
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Literal

import geopandas as gpd
import yaml
from dateutil.relativedelta import relativedelta


@dataclass(frozen=True)
class CaseStudy:
    country: str
    city: str
    date: datetime.date
    region_type: Literal["rural", "peri-urban", "urban"]
    shapefile_path: Path
    timezone: str

    @property
    def gdf(self):
        shapes_path = Path(__file__).parents[2] / "data" / "shapes"
        return gpd.read_file(shapes_path / self.region_type / self.shapefile_path)

    @staticmethod
    def all(namespace: SimpleNamespace) -> list["CaseStudy"]:
        """Return all CaseStudy instances in a SimpleNamespace."""
        return [v for v in vars(namespace).values() if isinstance(v, CaseStudy)]

    def __str__(self):
        """Return a formatted string summary of the case study."""
        lines = [
            f"Case Study: {self.city} ({self.region_type})",
            f"  Country     : {self.country}",
            f"  City        : {self.city}",
            f"  Date        : {self.date.strftime('%Y-%m-%d')}",
        ]
        return "\n".join(lines)

    def get_date_range(
        self,
        days: int = 0,
        months: int = 0,
        step_days=1,
    ) -> list[datetime.date]:

        start_date = self.date - relativedelta(days=days, months=months)
        end_date = self.date + relativedelta(days=days, months=months)

        current = start_date
        out = []
        while current <= end_date:
            out.append(current)
            current += datetime.timedelta(days=step_days)
        return out


def _load_case_studies(region_type):
    yaml_path = Path(__file__).parents[2] / "config" / "case-studies.yaml"
    with open(yaml_path) as f:
        data = yaml.safe_load(f)
    studies = {
        entry["id"]: CaseStudy(
            entry["country"],
            entry["city"],
            datetime.date.fromisoformat(str(entry["date"])),
            region_type,
            Path(entry["shapefile"]),
            entry["timezone"],
        )
        for entry in data[region_type]
    }
    return SimpleNamespace(**studies)


CaseStudiesRural: SimpleNamespace = _load_case_studies("rural")
CaseStudiesPeriUrban: SimpleNamespace = _load_case_studies("peri-urban")
CaseStudiesUrban: SimpleNamespace = _load_case_studies("urban")

# to get all: CaseStudy.all(CaseStudyRural)
