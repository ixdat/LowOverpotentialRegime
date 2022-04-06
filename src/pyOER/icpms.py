# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:56:04 2020

@author: scott
"""

from pathlib import Path
import json
import numpy as np
from matplotlib import pyplot as plt
from .constants import ICPMS_DIR, ICPMS_ID_FILE, ICPMS_CALIBRATION_ID_FILE
from .tools import singleton_decorator, CounterWithFile
from EC_MS import Chem


Measurement = None  # .measurement.Measurement imported first call avoid circular import


@singleton_decorator
class ICPMSCounter(CounterWithFile):
    """Counts icpms samples. 'id' increments the counter. 'last()' retrieves last id"""

    _file = ICPMS_ID_FILE


@singleton_decorator
class ICPMSCalCounter(CounterWithFile):
    """Counts icpms cals. 'id' increments the counter. 'last()' retrieves last id"""

    _file = ICPMS_CALIBRATION_ID_FILE


def all_icpms_points(icpms_dir=ICPMS_DIR):
    """returns an iterator that yields measurements in order of their id"""
    N_measurements = ICPMSCounter().last()
    for n in range(1, N_measurements):
        try:
            icpms_point = ICPMSPoint.open(n, icpms_dir=icpms_dir)
        except FileNotFoundError as e:
            # print(f"skipping {n} due to error = \n{e}")
            continue
        else:
            yield icpms_point


class ICPMSPoint:
    def __init__(
        self,
        i_id,
        ic_id,
        m_id,
        element,
        mass,
        signal,
        sampling_time,
        initial_volume,
        dilution,
        description,
    ):
        """Initiate an ICP-MS measurement

        Args:
            i_id (int or None): The ICPMS sample id. Defaults to ICPMSCounter.id
            ic_id (int): The id of the corresponding ICPMS calibration.
            m_id (int): The corresponding EC measurement id.
            element (str): The metal element measured
            mass (str): The mass at which the ICPMS was measuring for element.
            signal (float): The (averaged) ICPMS signal
            dilution (float): The ratio of concentration in the initial volume to
                concentration in the ICPMS sample as measured
            initial_volume (float): The volume of the containing electrolyte before any
                dilution, in [m^3]
            sampling_time (float): The time relative to Measurement.open(m_id).tstamp
                that the ICPMS sample was taken, in [s]
            description (str): useful info on how the sample was taken.
        """
        if not i_id:
            i_id = ICPMSCounter().id
        self.id = i_id
        self.ic_id = ic_id
        self.m_id = m_id
        self.element = element
        self.mass = mass
        self.signal = signal
        self.dilution = dilution
        self.initial_volume = initial_volume
        self.sampling_time = sampling_time
        self.description = description
        if ic_id is not None and ic_id > 0:  # I use -1 sometimes.
            self.calibration = ICPMSCalibration.open(ic_id)
        self._measurement = None  # opened when needed
        self.icpms_dir = ICPMS_DIR

    def as_dict(self):
        """Dictionary representation of the ICPMS spec"""
        self_as_dict = dict(
            id=self.id,
            ic_id=self.ic_id,
            m_id=self.m_id,
            element=self.element,
            mass=self.mass,
            signal=self.signal,
            dilution=self.dilution,
            initial_volume=self.initial_volume,
            sampling_time=self.sampling_time,
            description=self.description,
        )
        return self_as_dict

    def save(self, file_name=None):
        """Save the ICPMS measurement as .json in self.icpms_dir / file_name"""
        self_as_dict = self.as_dict()
        if not file_name:
            file_name = f"{self}.json"
        # print(f"saving measurement '{file_name}'")
        path_to_measurement = Path(self.icpms_dir) / file_name
        with open(path_to_measurement, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, file_name, icpms_dir=ICPMS_DIR):
        """Loads the measurement given its file path"""
        path_to_file = Path(icpms_dir) / file_name
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        if "id" in self_as_dict:
            self_as_dict["i_id"] = self_as_dict.pop("id")
        return cls(**self_as_dict)

    @classmethod
    def open(cls, i_id, icpms_dir=ICPMS_DIR):
        """Opens the measurement given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(icpms_dir).iterdir()
                if path.stem.startswith(f"i{i_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no icpms sample with id = {i_id}")
        return cls.load(path_to_file)

    def __repr__(self):
        s = (
            f"i{self.id} is {self.element} from "
            f"{self.sample} {self.description} on {self.date}"
        )

        s = s.replace("?", " hhmmmm ")  # can't save otherwise
        s = s.replace("->", "to")  # can't save otherwise
        return s

    @property
    def measurement(self):
        if not self._measurement and self.m_id:
            global Measurement
            if not Measurement:
                from .measurement import Measurement

            self._measurement = Measurement.open(self.m_id)
        return self._measurement

    @property
    def date(self):
        if self.m_id:
            return self.measurement.date

    @property
    def sample(self):
        if self.m_id:
            return self.measurement.sample_name

    @property
    def concentration(self):
        """Return the concentration of the element in original sample in [mol/m^3]"""
        c1 = self.calibration.signal_to_concentration(self.signal)
        c0 = c1 * self.dilution
        return c0

    @property
    def amount(self):
        """Return the amount element in initial volume in [mol]"""
        return self.concentration * self.initial_volume


class ICPMSCalibration:
    def __init__(
        self,
        ic_id,
        date,
        element,
        mass,
        ppbs,
        signals,
        wash_signals,
    ):
        """Initiate an ICP-MS calibration

        Args:
            ic_id (int or None): the id of the ICPMS calibration. Defaults to counter.
            date (str): the scott-formatted date the calibration was measured on
            element (str): the element
            mass (str): the mass
            ppbs (np.array): concentrations in [ppb]
            signals (np.array): signals in [counts]
        """
        if not ic_id:
            ic_id = ICPMSCalCounter().id
        self.id = ic_id
        self.date = date
        self.element = element
        self.mass = mass
        self.ppbs = ppbs
        self.signals = signals
        self.wash_signals = wash_signals
        self.icpms_dir = ICPMS_DIR
        self._calibration_curve = None

    def as_dict(self):
        """Dictionary representation of the ICPMS calibration"""
        self_as_dict = dict(
            id=self.id,
            date=self.date,
            element=self.element,
            mass=self.mass,
            ppbs=list(self.ppbs),
            signals=list(self.signals),
            wash_signals=list(self.wash_signals),
        )
        return self_as_dict

    def save(self, file_name=None):
        """Save the ICPMS calibration as .json in self.icpms_dir / file_name"""
        self_as_dict = self.as_dict()
        if not file_name:
            file_name = (
                f"ic{self.id} is icpms calibration for {self.element} on {self.date} "
            )
        # print(f"saving measurement '{file_name}'")
        path_to_measurement = Path(self.icpms_dir) / file_name
        with open(path_to_measurement, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, file_name, icpms_dir=ICPMS_DIR):
        """Loads the measurement given its file path"""
        path_to_file = Path(icpms_dir) / file_name
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        if "id" in self_as_dict:
            self_as_dict["ic_id"] = self_as_dict.pop("id")
        self_as_dict["ppbs"] = np.array(self_as_dict["ppbs"])
        self_as_dict["signals"] = np.array(self_as_dict["signals"])
        self_as_dict["wash_signals"] = np.array(self_as_dict["wash_signals"])
        return cls(**self_as_dict)

    def make_calibration_curve(self):
        """Make self._calibration_curve best fit of ln(self.ppbs) to ln(self.signals)"""
        ln_ppbs = np.log(self.ppbs)
        ln_signals = np.log(self.signals - self.bg)

        p = np.polyfit(ln_signals, ln_ppbs, deg=1)
        # print(p)  # debugging

        def calibration_curve(counts):
            """Return the concentration in [ppb] given signal in [counts]"""
            ln_counts = np.log(counts)
            ln_ppb = p[0] * ln_counts + p[1]
            ppb = np.exp(ln_ppb)
            return ppb

        self._calibration_curve = calibration_curve

    @classmethod
    def open(cls, ic_id, icpms_dir=ICPMS_DIR):
        """Opens the measurement given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(icpms_dir).iterdir()
                if path.stem.startswith(f"ic{ic_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no icpms calibration with id = {ic_id}")
        return cls.load(path_to_file)

    @property
    def bg(self):
        """Background signal / [counts]"""
        return np.mean(self.wash_signals)

    @property
    def calibration_curve(self):
        """The calibration curve converting counts to parts per billion"""
        if not self._calibration_curve:
            self.make_calibration_curve()
        return self._calibration_curve

    @property
    def dl_signal(self):
        """Detection limit signal / [counts]"""
        return self.bg + 3 * np.std(self.wash_signals)

    @property
    def dl_concentration(self):
        """Detection limit concentration / [mol/m^3]"""
        return self.calibration_curve(self.dl_signal)

    def signal_to_concentration(self, signal):
        """Return concentration in [mol/m^3] of ICPMS sample given its signal"""
        ppb_concentration = self.calibration_curve(signal - self.bg)

        if ppb_concentration < self.dl_concentration:
            print(
                f"WARNING! ICPMS implied {self.element} concentration in ICPMS sample "
                + f"is {ppb_concentration} ppb, which is below "
                + f"the detection limit of {self.dl_concentration} ppb"
            )
        kg_per_m3 = ppb_concentration * 1e-6
        kg_per_mol = Chem.get_mass(self.element) * 1e-3
        concentration = kg_per_m3 / kg_per_mol

        return concentration

    def plot_calibration(self, ax=None):
        """Plot the ICPMS calibration (as fig A.4 of Scott's PhD thesis)

        Args:
            ax (plt.Axis):
        """
        ppbs = self.ppbs
        signals = self.signals
        calibration_curve = self.calibration_curve
        if not ax:
            fig, ax = plt.subplots()
            ax.set_xlabel(f"amount of {self.element} / [ppb]")
            ax.set_ylabel(f"counts at {self.mass}")
            ax.set_xscale("log")
            ax.set_yscale("log")

        ax.plot(ppbs, signals, "ks", markersize=7)
        x_fit, y_fit = calibration_curve(signals), signals

        if self.wash_signals is not None:
            wash = self.wash_signals
            mean = np.mean(wash)
            std = np.std(wash)

            y_fit = np.append(mean, y_fit)
            x_fit = np.append(calibration_curve(mean), x_fit)
            ax.plot(x_fit, y_fit, "r--")

            xlim = ax.get_xlim()
            ax.plot(xlim, [mean, mean], "k--")
            ax.plot(xlim, [mean + 3 * std, mean + 3 * std], "k:")
            ax.set_xlim(xlim)

        else:
            ax.plot(x_fit, y_fit, "r--")
