"""This module links Measuremnt to sensitivity factor and electrolyte isotope ratio"""

from pathlib import Path
import json
import numpy as np
from matplotlib import pyplot as plt
from EC_MS.utils.extraction_class import Extraction
from EC_MS import Chem
from .measurement import Measurement
from .constants import CALIBRATION_DIR, CALIBRATION_ID_FILE, PROJECT_START_TIMESTAMP
from .tools import (
    singleton_decorator,
    CounterWithFile,
    fit_exponential,
)


if not CALIBRATION_DIR.exists():
    # This block of code makes it so I can delete calibrations/ and start over quickly.
    Path.mkdir(CALIBRATION_DIR)
    with open(CALIBRATION_ID_FILE, "w") as counter_file:
        counter_file.write("0")


def all_calibrations(calibration_dir=CALIBRATION_DIR):
    """returns an iterator that yields measurements in order of their id"""
    N_calibrations = CalibrationCounter().last()
    for n in range(1, N_calibrations + 1):
        try:
            calibration = Calibration.open(n, calibration_dir=calibration_dir)
        except FileNotFoundError as e:
            print(f"all_calibrations() is skipping {n} due to error = \n{e}")
        else:
            yield calibration


@singleton_decorator  # ... this is not really necessary as the file is read each time.
class CalibrationCounter(CounterWithFile):
    """Counts calibrations. 'id' increments the counter. 'last()' retrieves last id"""

    _file = CALIBRATION_ID_FILE


calibration_counter = CalibrationCounter()


class Calibration:
    """Class for referencing calibration raw data and storing a calibration result"""

    # -------- all of this should be a template like a baseclass ---------------- #

    def __init__(
        self,
        c_id=None,
        m_id=None,
        tspan=None,
        cal_tspans=None,
        t_bg=None,
        F=None,
        alpha=None,
        category=None,
        isotope=None,
        **kwargs,
    ):
        """Initiate the calibration

        Args:
            c_id (int): calibration number (self.id)
            m_id (int): measurement number
            tspan (tspan): timespan of the measurement used for the calibration
            F (dict): sensitivity factors
            alpha (float): isotope ratio in electrolyte
            category (string): typically either "good" or "bad"
            isotope (string): typically "16" or "18"
        """
        if c_id is None:
            c_id = calibration_counter.id
        self.id = c_id
        self.m_id = m_id
        self.tspan = tspan
        self.cal_tspans = cal_tspans
        self.t_bg = t_bg
        if F is None:
            F = {}
        self.F = F
        self.alpha = alpha
        self.category = category
        self.isotope = isotope
        self._measurement = None  # measurement is a managed property
        self._extraction = None  # extraction is a managed property
        self.extra_stuff = kwargs
        self.name = self.make_name()

    def as_dict(self):
        """Return the dictionary representation of the calibration"""
        self_as_dict = {
            "c_id": self.id,
            "m_id": self.m_id,
            "tspan": self.tspan,
            "cal_tspans": self.cal_tspans,
            "t_bg": self.t_bg,
            "F": self.F,
            "alpha": self.alpha,
            "category": self.category,
            "isotope": self.isotope,
        }
        return self_as_dict

    @classmethod
    def load(cls, file_name, calibration_dir=CALIBRATION_DIR):
        path_to_file = Path(calibration_dir) / file_name
        with open(path_to_file) as f:
            self_as_dict = json.load(f)
        self_as_dict.update(file_loaded_from=path_to_file)
        if "id" in self_as_dict:
            self_as_dict["c_id"] = self_as_dict.pop("id")
        return cls(**self_as_dict)

    @classmethod
    def open(cls, c_id, calibration_dir=CALIBRATION_DIR):
        """Opens the measurement given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(calibration_dir).iterdir()
                if path.stem.startswith(f"c{c_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no calibration with id = m{c_id}")
        return cls.load(path_to_file)

    def save(self, file_name=None, calibration_dir=CALIBRATION_DIR):
        self_as_dict = self.as_dict()
        if not file_name:
            file_name = self.make_name() + ".json"
        with open(Path(calibration_dir) / file_name, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    def save_with_rename(self, file_name=None):
        if "file_loaded_from" in self.extra_stuff:
            Path(str(self.extra_stuff["file_loaded_from"])).unlink()
        self.save(file_name=file_name)

    # -------- now starts the calibration-specific stuff ---------------- #
    @property
    def measurement(self):
        if not self._measurement:
            self._measurement = Measurement.open(self.m_id)
        return self._measurement

    @property
    def meas(self):
        return self.measurement.meas

    def make_name(self):
        c_id = self.id
        date = self.measurement.date
        sample = self.measurement.sample_name
        category = self.category
        name = f"c{c_id} is a {category} cal with {sample} on {date}"
        return name

    def calibration_curve(self, ax=None, **kwargs):
        """Calibrate O2 using EC_MS.Dataset.calibration_curve() and self.tspans"""
        F_O2 = 0
        for mass in ["M32", "M34", "M36"]:
            O2, ax = self.meas.ecms_calibration_curve(
                mol="O2",
                mass=mass,
                n_el=4,
                tspans=self.cal_tspans,
                t_bg=self.t_bg,
                out=["Molecule", "ax"],
                ax=ax,
                **kwargs,
            )
            F_O2 += O2.F_cal
        self.F["O2"] = F_O2
        return F_O2

    def cal_F_O2(self):
        """calibrate the O2 signal based on assumption of OER during self.tspan"""
        Y_cum = 0
        for mass in ["M32", "M34", "M36"]:
            x, y = self.meas.grab_signal(
                mass=mass,
                tspan=self.tspan,
                t_bg=self.t_bg,
                unit="A",
            )
            Y = np.trapz(y, x)
            Y_cum += Y
        t, I = self.meas.grab("current", tspan=self.tspan)
        I *= 1e-3  # from [mA] to [A]
        Q = np.trapz(I, t)
        n = Q / (4 * Chem.Far)
        F_O2 = Y_cum / n

        self.F["O2"] = F_O2

        return F_O2

    def cal_alpha(self):
        """calibrate the isotope ratio based on OER during self.tspan"""

        if "16" in self.isotope:
            x_34, y_34 = self.meas.grab_signal("M34", tspan=self.tspan)
            x_32, y_32 = self.meas.grab_signal("M32", tspan=self.tspan)
            beta = np.mean(y_34) / np.mean(y_32)
            alpha = 2 / (2 + beta)
        elif "18" in self.isotope:
            x_34, y_34 = self.meas.grab_signal("M34", tspan=self.tspan)
            x_36, y_36 = self.meas.grab_signal("M36", tspan=self.tspan)
            gamma = np.mean(y_34) / np.mean(y_36)
            alpha = gamma / (2 + gamma)
        else:
            alpha = None

        self.alpha = alpha

        return alpha


class CalibrationSeries:
    """A class to describe the trend of calibrations and predict F based on tstamp"""

    def __init__(self, c_id_list=None, tau=None, y0=None, y1=None):
        """Initiate with a list of calibration ids and, optionally, fit params

        Args:
            c_id_list (list of int): The list of id's for the calibrations
                Defaults to all calibrations categorized as "good"
            tau (float): time constant for sensitivity factor convergence
            y0 (float): convergence value of sensitivity factor
            y1 (float): value of sensitivity factor at t=0
        """
        if not c_id_list:
            c_id_list = [c.id for c in all_calibrations() if c.category == "good"]
        self.c_id_list = c_id_list
        self.tau = tau
        self.y0 = y0
        self.y1 = y1
        # function: returns sensitivity factor given timestamp
        self._F_of_tstamp = None

    def as_dict(self):
        self_as_dict = {
            "c_id_list": self.c_id_list,
            "tau": self.tau,
            "y0": self.y0,
            "y1": self.y1,
        }
        return self_as_dict

    def save(self):
        self_as_dict = self.as_dict()
        with open(CALIBRATION_DIR / "TREND.json", "w") as f:
            json.dump(self_as_dict, f)

    @classmethod
    def load(cls):
        with open(CALIBRATION_DIR / "TREND.json", "r") as f:
            self_as_dict = json.load(f)
        return cls(**self_as_dict)

    def calibrations(self):
        """Generator yielding the calibrations in the CalibrationSeries"""
        for c_id in self.c_id_list:
            yield Calibration.open(c_id)

    def sensitivity_trend(self, ax="new"):
        """return and plot the sensitivity factors vs time"""

        time_vec = np.array([])
        F_vec = np.array([])

        if ax == "new":
            fig, ax = plt.subplots()
            ax.set_xlabel("project time / [days]")
            ax.set_ylabel("O2 sensitivity / [C/mol]")

        t_project_start = PROJECT_START_TIMESTAMP
        for calibration in self.calibrations():

            t = calibration.meas.tstamp - t_project_start
            F = calibration.F["O2"]
            time_vec = np.append(time_vec, t)
            F_vec = np.append(F_vec, F)

            if ax:
                if calibration.isotope == 16:
                    color = "k"
                elif calibration.isotope == 18:
                    color = "g"
                else:
                    color = "r"  # something's wrong

                sample = calibration.measurement.sample_name
                if sample == "Trimi1":
                    marker = "s"  # Platinum as squares
                elif "Jazz" in sample or "Folk" in sample or "Emil" in sample:
                    marker = "o"  # Iridium and IrO2 as circles
                else:
                    marker = "*"  # something's wrong

                ax.plot(t / (24 * 60 * 60), F, color=color, marker=marker)

        return time_vec, F_vec

    def make_F_of_tstamp(self):
        """Make a function for returning F given a tstamp"""
        tau, y0, y1 = self.tau, self.y0, self.y1
        if tau is None or y0 is None or y1 is None:
            F_of_tstamp = None
        else:

            def F_of_tstamp(tstamp):
                t = tstamp - PROJECT_START_TIMESTAMP
                return y0 + (y1 - y0) * np.exp(-t / tau)

        self._F_of_tstamp = F_of_tstamp
        return F_of_tstamp

    @property
    def F_of_tstamp(self):
        if not self._F_of_tstamp:
            return self.make_F_of_tstamp()
        return self._F_of_tstamp

    def fit_exponential(self, ax="new"):
        time_vec, F_vec = self.sensitivity_trend(ax=ax)

        tau, y0, y1 = fit_exponential(time_vec, F_vec)

        if ax:
            ax = plt.gca()
            t_fit = np.arange(0, 2 * 365 * 24 * 60 * 60, 1000)  # 2 years from start
            F_fit = y0 + (y1 - y0) * np.exp(-t_fit / tau)
            ax.plot(t_fit / (24 * 60 * 60), F_fit, "k--")
        self.tau = tau
        self.y0 = y0
        self.y1 = y1
        self.make_F_of_tstamp()

        return tau, y0, y1

    def calc_F_from_trend(self, measurement):
        """Return the sensitivity factor for measurement based on calibration trend

        Args:
            measurement (Measurement): the measurement
        """
        tstamp = measurement.meas.tstamp
        return self.F_of_tstamp(tstamp)
