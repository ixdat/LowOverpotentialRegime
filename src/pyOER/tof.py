"""this module implements methods and classes around point results"""
from pathlib import Path
import json
import numpy as np

from .calc import (
    calc_OER_rate,
    calc_dissolution_rate,
    calc_exchange_rate,
    calc_potential,
    calc_current,
)
from .tools import singleton_decorator, CounterWithFile
from .constants import TOF_DIR, TOF_ID_FILE, FARADAY_CONSTANT
from .experiment import open_experiment


@singleton_decorator
class TOFCounter(CounterWithFile):
    """Counts measurements. 'id' increments the counter. 'last()' retrieves last id"""

    _file = TOF_ID_FILE


def all_tofs(tof_dir=TOF_DIR):
    """returns an iterator that yields measurements in order of their id"""
    N_tofs = TOFCounter().last()
    for n in range(1, N_tofs + 1):
        try:
            tof = TurnOverFrequency.open(n, tof_dir=tof_dir)
        except FileNotFoundError as e:
            continue
        else:
            yield tof


def all_tof_sets(tof_dir=TOF_DIR):
    tof_sets = {}
    for tof in all_tofs(tof_dir=tof_dir):
        e_id = tof.e_id
        tspan = tuple(int(t) for t in tof.tspan)
        if not (e_id, tspan) in tof_sets:
            tof_sets[(e_id, tspan)] = TurnOverSet()
        tof_sets[(e_id, tspan)].add_tof(tof)
    yield from tof_sets.values()


class TurnOverSet:
    def __init__(
        self,
        t_ids=None,
    ):
        """Initiate a set of turn-over-frequencies taken from one point in an experiment

        Args:
            t_ids (dict): {tof_type: t_id} where tof_type is the action the tof
                describes ("activity", "exchange", "dissolution") and t_id is the id
        """
        self.t_ids = t_ids or {}
        self._tofs = {}
        self.experiment = None
        self.tspan = None

    def __contains__(self, item):
        return item in self.t_ids

    def __repr__(self):
        return f"TurnOverSet({self.t_ids})"

    def add_tof(self, tof):
        tof_type = tof.tof_type
        self.t_ids[tof_type] = tof.id
        self._tofs[tof_type] = tof
        if not self.experiment:
            self.experiment = tof.experiment
        elif not (tof.experiment.id == self.experiment.id):
            raise TypeError(f"can't add {tof} to {self} as experiment is not the same")
        if not self.tspan:
            self.tspan = tof.tspan
        elif not (tof.tspan[0] == self.tspan[0]):
            raise TypeError(f"can't add {tof} to {self} as tspans are not the same")

    def get_tof(self, item):
        if item in self._tofs:
            return self._tofs[item]
        elif item in self.t_ids:
            if self.experiment:
                self._tofs[item] = TurnOverFrequency.open(
                    self.t_ids[item], experiment=self.experiment
                )
            else:
                self._tofs[item] = TurnOverFrequency.open(self.t_ids[item])
                self.experiment = self._tofs[item].experiment
            return self._tofs[item]
        raise KeyError(f"{self} does not have tof for {item}")

    def __getitem__(self, item):
        return self.get_tof(item)

    def __iter__(self):
        yield from self._tofs.values()

    def __getattr__(self, item):
        try:
            return self.get_tof(item)
        except KeyError as e:
            raise AttributeError(e)

    @property
    def sample(self):
        return self.experiment.sample

    @property
    def sample_name(self):
        return self.experiment.sample_name


class TOFCollection:
    """A group of turnoverfrequencies"""

    def __init__(self, tof_list=None):
        self.tof_list = tof_list

    def __iter__(self):
        yield from self.tof_list

    def mean_rate(self):
        return np.mean(np.array([tof.rate for tof in self]))

    def std(self):
        return np.mean(np.array([tof.rate for tof in self]))


class TurnOverFrequency:
    def __init__(
        self,
        tof_type=None,
        rate=None,
        tof=None,
        current=None,
        potential=None,
        e_id=None,
        experiment=None,
        sample_name=None,
        tspan=None,
        r_id=None,
        rate_calc_kwargs=None,
        description=None,
        t_id=None,
        amount=None,
    ):
        """Iinitiate a TurnOverFrequency

        Args:
            tof_type (str): The type of TOF. Options are 'activity', 'exchange', and
                'dissolution'.
            rate (float): The un-normalized rate, if known, in [mol/s]
            e_id (int): The id of the associated experiment
            experiment (Experiment): optionally, the Experiment itself can be given to
                save time.
            sample_name (str): Sample name. Only needed if no experiment is given.
            tspan (timespan): The time interval over which to integrate/average
            r_id (int): The id of the associated roughness measurement
            rate_calc_kwargs (dict): Extra kwargs for the relevant rate calc. function.
            description (str): free-form description of the TOF point
            t_id (int): The principle key. Defaults to incrementing the counter
        """
        self.tof_type = tof_type
        self.e_id = e_id
        self.tspan = tspan
        self.r_id = r_id
        self._experiment = experiment
        self._rate = rate
        self._tof = tof
        self._potential = potential
        self._current = current
        self._sample_name = sample_name
        self.description = description
        self.rate_calc_kwargs = rate_calc_kwargs or {}
        self.id = t_id or TOFCounter().id
        self._rate = rate
        self._amount = amount

    def as_dict(self):
        """The dictionary represnetation of the TOF's metadata"""
        return dict(
            tof_type=self.tof_type,
            rate=self._rate,  # result!
            tof=self._tof,  # result!
            potential=self._potential,  # result!
            current=self._current,  # result!
            e_id=self.e_id,
            tspan=self.tspan,
            r_id=self.r_id,
            rate_calc_kwargs=self.rate_calc_kwargs,
            description=self.description,
            t_id=self.id,
            amount=self._amount,
            sample_name=self._sample_name,
        )

    def save(self):
        """Save the TOF's metadata to a .json file"""
        self_as_dict = self.as_dict()
        path_to_file = TOF_DIR / f"{self}.json"
        with open(path_to_file, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    @classmethod
    def load(cls, path_to_file, **kwargs):
        """Load a TOF from the metadata stored in a file"""
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    @classmethod
    def open(cls, t_id, tof_dir=TOF_DIR, **kwargs):
        """Opens the measurement given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(tof_dir).iterdir()
                if path.stem.startswith(f"t{t_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no TurnOverFrequency with id = t{t_id}")
        return cls.load(path_to_file, **kwargs)

    def __repr__(self):
        if not self.experiment:
            return f"t{self.id} is {self.tof_type} without experiment in pyOER"
        return f"t{self.id} is {self.tof_type} on {self.sample_name} on {self.date}"

    # -------- table-joining properties ------------ #
    @property
    def experiment(self):
        if not self._experiment:
            if not self.e_id:
                print(f"TOF with id={self.id} has no attached experiment.")
                return None
            self._experiment = open_experiment(self.e_id)
        return self._experiment

    @property
    def measurement(self):
        if not self.experiment:
            return None
        return self.experiment.measurement

    @property
    def date(self):
        if not self.measurement:
            return None
        return self.measurement.date

    @property
    def sample(self):
        return self.measurement.sample

    @property
    def sample_name(self):
        if self._sample_name:
            return self._sample_name
        if not self.measurement:
            return None
        return self.measurement.sample_name

    @property
    def element(self):
        return self.sample.element

    @property
    def t_interval(self):
        """float: The length of electrolysis time covered by the TOF"""
        return self.tspan[-1] - self.tspan[0]

    @property
    def rate_calculating_function(self):
        """The function that this TOF uses to calculate its rate"""
        if self.tof_type == "activity":
            return calc_OER_rate
        if self.tof_type == "exchange":
            return calc_exchange_rate
        if self.tof_type == "dissolution":
            return calc_dissolution_rate
        raise TypeError(f"no associated rate function for tof_type={self.tof_type}")

    def calc_rate(self, **kwargs):
        """Calculate and return the relevant rate in [mol/s]"""
        if not self.experiment:
            return
        rate_calc_kwargs = self.rate_calc_kwargs
        rate_calc_kwargs.update(kwargs)
        rate = self.rate_calculating_function(
            experiment=self.experiment, tspan=self.tspan, **rate_calc_kwargs
        )
        self._rate = rate
        return rate

    def calc_amount(self, **kwargs):
        """Calculate and return the relevant amout in [mol]"""
        if not self.experiment:
            return
        amount = self.calc_rate(**kwargs) * self.t_interval  # noqa
        self._amount = amount
        return amount

    @property
    def rate(self):
        """The rate (activity, dissolution, or exchange) in [mol/s]"""
        if not self._rate:
            self.calc_rate()
        return self._rate

    def calc_tof(self):
        if not self.experiment:
            return
        self._tof = self.rate / self.experiment.n_sites
        return self._tof

    @property
    def tof(self):
        """The estimated turn-over frequency in [s^-1]

        TOF is estimated (via self.rate / self.experiment.n_sites) as:
        tof = partial_capacitance_normalized_current / (4 * FARADAYS_CONSTANT)   \
            * STANDARD_SPECIFIC_CAPACITANCE / STANDARD_SITE_DENSITY
        giving units [A/F] / [C/mol] * [F/cm^2] / [mol/cm^2] =  [A/C] = [s^-1]
        """
        if not self._tof:
            self.calc_tof()
        return self._tof

    def calc_potential(self):
        if not self.experiment:
            return
        self._potential = calc_potential(self.experiment, self.tspan)
        return self._potential

    @property
    def amount(self):
        """The rate (activity, dissolution, or exchange) in [mol]"""
        if not self._amount:
            self.calc_amount()
        return self._amount

    @property
    def potential(self):
        """The potential vs RHE in [V]"""
        if not self._potential:
            self.calc_potential()
        return self._potential

    def calc_current(self):
        self._current = calc_current(self.experiment, self.tspan)
        return self._current

    @property
    def current(self):
        """The potential vs RHE in [V]"""
        if not self._current:
            self.calc_current()
        return self._current

    def get_tof_triplet(self):
        act_tof = None
        diss_tof = None
        exc_tof = None
        for tof in self.experiment.get_tofs():
            if not tof.tspan == self.tspan:
                continue
            if tof.tof_type == "activity":
                act_tof = tof
            elif tof.tof_type == "dissolution":
                diss_tof = tof
            elif tof.tof_type == "exchange":
                exc_tof = tof
        return act_tof, diss_tof, exc_tof

    def calc_faradaic_efficiency(self, n_el=4, mol=None):
        rate = self.calc_rate(mol=mol) if mol else self.rate
        FE = rate * n_el * FARADAY_CONSTANT / self.current
        return FE
