"""this module will define the sample object grouping experiments and results"""

import json

from .constants import SAMPLE_DIR, STANDARD_SITE_DENSITY, STANDARD_SPECIFIC_CAPACITANCE
from .measurement import Measurement, all_measurements

SAMPLE_TYPES = {
    "Ru": {
        "hydrous": [
            "Taiwan1G",
        ],
        "metallic": ["Melih", "Bernie"],
        "foam": "Evans",
        "rutile": ["Reshma4", "Maundy", "Stoff", "Sofie", "Mette", "John"],
        "amorphous": ["Reshma1", "Nancy", "Easter", "Taiwan"],
    },
    "Ir": {
        "hydrous": ["Legend1C", "Decade1G"],
        "metallic": "Jazz",
        "rutile": ["Folk", "Champ", "Legend"],
        "amorphous": ["Goof", "Decade"],
    },
    "Pt": {"metallic": "Trimi"},
}

SAMPLE_ISOTOPES = {
    "16": ["Reshma", "Folk"],
    "18": [
        "Maundy",
        "Stoff",
        "Sofie",
        "Mette",
        "John",
        "Easter",
        "Taiwan",
        "Champ",
        "Legend",
        "Goof",
        "Decade",
    ],
    "(check!)": ["Melih", "Bernie", "Evans", "Nancy", "Jazz", "Trimi"],
}


def get_element_and_type(name, get="both"):
    for element, oxide_types in SAMPLE_TYPES.items():
        for oxide_type, sample_names in oxide_types.items():
            for sample_name in (
                [sample_names] if isinstance(sample_names, str) else sample_names
            ):
                if name.startswith(sample_name):
                    if get == "element":
                        return element
                    elif get in ["type", "crystallinity", "oxide_type"]:
                        return oxide_type
                    return element, oxide_type
    return "unknown"


def get_isotope(name):
    for isotope, sample_names in SAMPLE_ISOTOPES.keys():
        for sample in [sample_names] if isinstance(sample_names, str) else sample_names:
            if name.startswith(sample):
                return isotope
    return "unknown"


class Sample:
    def __init__(self, name, synthesis_date=None, history=None):
        self.name = name
        self.synthesis_date = synthesis_date
        self.history = history or {}

    def as_dict(self):
        return dict(
            name=self.name, synthesis_date=self.synthesis_date, history=self.history
        )

    def save(self):
        file = SAMPLE_DIR / (self.name + ".json")
        self_as_dict = self.as_dict()
        with open(file, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, path_to_file):
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        return cls(**self_as_dict)

    @classmethod
    def open(cls, name):
        path_to_file = SAMPLE_DIR / (name + ".json")
        return cls.load(path_to_file)

    def __repr__(self):
        return f"{self.__class__}({self.name})"

    @property
    def isotope(self):
        return get_isotope(self.name)

    @property
    def element(self):
        return get_element_and_type(self.name, get="element")

    @property
    def site_density(self):
        """Density of sites in mol/cm^2. TODO: module dictionary with elements"""
        return STANDARD_SITE_DENSITY  # 1 site/nm^2 in mol/cm^2

    @property
    def specific_capacitance(self):
        """Specific capacitance in Farad/cm^2. TODO: module dictionary with elements"""
        return STANDARD_SPECIFIC_CAPACITANCE

    @property
    def oxide_type(self):
        return get_element_and_type(self.name, get="oxide_type")

    def describe(self):
        return f"{self.name} is {self.element} {self.oxide_type} with oxygen {self.isotope}"

    @property
    def description(self):
        return self.describe()

    @property
    def measurement_ids(self):
        m_ids = []
        for m in all_measurements():
            if m.sample_name == self.name:
                m_ids += [m.id]
        return m_ids

    @property
    def measurements(self):
        m_ids = self.measurement_ids
        ms = [Measurement.open(m_id) for m_id in m_ids]
        # sort by tstamp:
        ms = sorted(ms)
        return ms

    @property
    def tofs(self):
        from .tof import all_tofs

        return [tof for tof in all_tofs() if tof.sample_name == self.name]

    @property
    def tof_sets(self):
        from .tof import all_tof_sets

        return [t_set for t_set in all_tof_sets() if t_set.sample_name == self.name]

    @classmethod
    def all_samples(self):
        all_samples = [sample.stem for sample in SAMPLE_DIR.rglob('*.json')]
        all_samples.sort()
        return all_samples
