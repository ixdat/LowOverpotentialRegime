from pathlib import Path

# -------------- natural constants -------------- #
STANDARD_ALPHA = 0.99804  # the natural ratio ^{16}O/(^{16}O + ^{18}O) of oxygen atoms
# determined from control standard experiments in generate_standard_experiments.py


# -------------- project-specific constants -------------- #
PROJECT_START_TIMESTAMP = 1533502800  # August 25, 2018 in unix time

EXPERIMENT_TAGS = {
    "y": "yes, purely systematic",  # 30 minutes at one current density, taken out at
    # OCP just after, ICPMS samples at ~2, 10, 20 minutes.
    "s": "starts systematic",  # starts with one current density and ICPMS samples at
    # 2 and 10 minutes, but play around after
    "k": "shortened systematic (<30 minutes)",  # systematic but for less than 30
    # minutes
    "c": "composite systematic",  # short systematic measurements, typically 10
    # minutes, at different current densities, with ICPMS samples taken interspersed
    "p": "some kind of constant-potential systematic",
    "f": "failed systematic",  # something is wrong.
    "n": "not systematic at all.",
    "d": "duplicate",
    "b": "broken",  # something is wrong with the file
    "a": "activity",
    "m": "missing (start, CV's) for activity",
    "t": "short (truncated) activity",
    "o": "oxygen-saturated activity",
    "l": "long-term (overnight) activity",
    "16": "in un-labeled (H2(16)O) electrolyte",
    "18": "in labeled (H2(18)O) electrolyte",
    "q": "[Quit and save progress]",
}
STANDARD_EXPERIMENT_TAGS = ["y", "k", "s", "c"]

# -------------- table stuff (directories and counter files) -------------- #
PROJECT_DIR = Path(__file__).absolute().parent.parent.parent
EXPORT_DATA_DIR = PROJECT_DIR / "data_export"
DEFAULT_DATA_DIR = PROJECT_DIR / "data_raw"
# Default location to store raw data. If you store the raw data somewhere else,
#   please create a settings.py file and define `DATA_DIR` in that file.

ELOG_DIR = PROJECT_DIR / "tables/elog"

SAMPLE_DIR = PROJECT_DIR / "tables/samples"

MEASUREMENT_DIR = PROJECT_DIR / "tables/measurements"
MEASUREMENT_ID_FILE = MEASUREMENT_DIR / "LAST_MEASUREMENT_ID.pyoer20"

CALIBRATION_DIR = PROJECT_DIR / "tables/calibrations"
CALIBRATION_ID_FILE = CALIBRATION_DIR / "LAST_CALIBRATION_ID.pyoer20"

ICPMS_DIR = PROJECT_DIR / "tables/icpms"
ICPMS_ID_FILE = ICPMS_DIR / "LAST_ICPMS_ID.pyoer20"
ICPMS_CALIBRATION_ID_FILE = ICPMS_DIR / "LAST_ICPMS_CALIBRATION_ID.pyoer20"

EXPERIMENT_DIR = PROJECT_DIR / "tables/experiments"
EXPERIMENT_ID_FILE = EXPERIMENT_DIR / "LAST_EXPERIMENT_ID.pyOER20"

TOF_DIR = PROJECT_DIR / "tables/tofs"
TOF_ID_FILE = TOF_DIR / "LAST_TOF_ID.pyOER20"

AVOGADROS_CONSTANT = 6.02217e23  # [1/mol]
FARADAY_CONSTANT = 96485  # [C/mol]
GAS_CONSTANT = 8.31446  # [J/(mol*K)]
STANDARD_TEMPERATURE = 298.15  # [K]

STANDARD_SITE_DENSITY = 5 / (1e-7**2) / AVOGADROS_CONSTANT
# ^ 5 sites per nm^2 in [mol/cm^2]
# five sites per sq nm is based on Ru(110) cus site density, from Rao2017.

STANDARD_SPECIFIC_CAPACITANCE = 115e-6  # specific capacitance in [Far/cm^2]
# 115 uF/cm^2 is based on average of single crystal capacitances.

STANDARD_ELECTRODE_AREA = 0.196  # electrode area in [cm^2]

STANDARD_SYMMETRY_FACTOR = 0.5  # default exponential weight in Tafel eq., unitless.
