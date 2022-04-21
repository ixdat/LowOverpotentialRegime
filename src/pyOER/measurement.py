"""Define the Measurement class containing metadata and pointers to raw data.

"""
from pathlib import Path, PureWindowsPath, PurePosixPath
import json
import re
import time
import datetime

from ixdat import Measurement as Meas
from .constants import MEASUREMENT_DIR, MEASUREMENT_ID_FILE, STANDARD_ELECTRODE_AREA
from .tools import singleton_decorator, CounterWithFile, FLOAT_MATCH
from .settings import DATA_DIR


@singleton_decorator
class MeasurementCounter(CounterWithFile):
    """Counts measurements. 'id' increments the counter. 'last()' retrieves last id"""

    _file = MEASUREMENT_ID_FILE


def all_measurements(measurement_dir=MEASUREMENT_DIR):
    """returns an iterator that yields measurements in order of their id"""
    N_measurements = MeasurementCounter().last()
    for n in range(1, N_measurements + 1):
        try:
            measurement = Measurement.open(n, measurement_dir=measurement_dir)
        except FileNotFoundError as e:
            continue
            # print(f"itermeasurement skipping {n} due to error = \n{e}")
        else:
            yield measurement


class Measurement:
    """Wrapper around ixdat measurement w metadata for its context in this project."""

    def __init__(
        self,
        m_id=None,
        measurement_dir=MEASUREMENT_DIR,
        copied_at=None,
        name=None,
        sample=None,
        technique=None,
        isotope=None,
        date=None,
        analysis_date=None,
        raw_data_path=None,
        exported_data_path=None,
        meas=None,
        linked_measurements=None,
        elog_number=None,
        elog=None,
        EC_tag=None,
        category=None,
        **kwargs,
    ):
        """Initiate Measurement object

        Intended use is to load the ixdat measurement (meas) separately using
        Measurement.load_data()

        Args:
            m_id (int): the unique id of the measurement
            name (str): the name of the measurement
            measurement_dir (Path-like): where to SAVE the measurement metadata
            copied_at (float): the time at which the meas was read
            raw_data_path (Path-like): path to file to load the raw data from pkl
            exported_data_path (Path-like): path to file to SAVE the raw data as pkl
            meas (ixdat.techniques.ECMSMeasurement): the ixdat object with the data
            linked_measurements (dict): measurements to link to this one
            kwargs (dict): gets added, not used. Here so that I can add extra stuff when
                saving just to improve readability of the json
        """
        if not m_id:
            m_id = MeasurementCounter().id
        self.id = m_id
        self.name = name
        self.sample_name = sample
        self.technique = technique
        self.isotope = isotope
        self.date = date
        self.analysis_date = analysis_date
        self.measurement_dir = measurement_dir
        self.copied_at = copied_at  # will be replaced by time.time()
        self.raw_data_path = raw_data_path
        self.exported_data_path = exported_data_path
        self._meas = meas  # meas is a managed property
        self.linked_measurements = linked_measurements
        self.extra_stuff = kwargs
        self.elog_number = elog_number
        self._elog = elog  # elog is a managed property
        self.EC_tag = EC_tag
        self.category = category

    def as_dict(self):
        self_as_dict = dict(
            id=self.id,
            name=self.name,
            sample=self.sample_name,
            technique=self.technique,
            isotope=self.isotope,
            date=self.date,
            analysis_date=self.analysis_date,
            copied_at=self.copied_at,
            raw_data_path=str(self.raw_data_path),
            exported_data_path=str(self.exported_data_path),
            linked_measurements=self.linked_measurements,
            elog_number=self.elog_number,
            EC_tag=self.EC_tag,
            category=self.category,
            # do not put the meas into self_as_dict!
        )
        if self.copied_at:  # just for human-readability of measurement .json
            self_as_dict["copied_on"] = datetime.datetime.fromtimestamp(
                self.copied_at
            ).strftime("%Y-%m-%d %H:%M:%S")
        return self_as_dict

    def save(self, file_name=None, save_dataset=False):
        self_as_dict = self.as_dict()
        if not file_name:
            if not self.name:
                self.make_name()
            file_name = self.name + ".json"
        print(f"saving measurement '{file_name}'")
        path_to_measurement = Path(self.measurement_dir) / file_name
        with open(path_to_measurement, "w") as f:
            json.dump(self_as_dict, f, indent=4)
        if save_dataset:
            self.export_data()

    def save_with_rename(self, file_name=None):
        if "file_loaded_from" in self.extra_stuff:
            Path(str(self.extra_stuff["file_loaded_from"])).unlink()
        self.save(file_name=file_name)

    @classmethod
    def load(cls, file_name, measurement_dir=MEASUREMENT_DIR):
        """Loads the measurement given its file path"""
        path_to_file = Path(measurement_dir) / file_name
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        self_as_dict.update(file_loaded_from=path_to_file)
        if "id" in self_as_dict:
            self_as_dict["m_id"] = self_as_dict.pop("id")

        # Change hardcoded path to relative path
        # Should work both Windows -> Linux/Mac and Linux/Mac -> Windows
        # as long as the path specified in measurement file is absolute
        # Only tested Windows -> Linux and Windows -> Windows
        for key in ["raw_data_path", "exported_data_path"]:
            try:
                self_as_dict[key]
            except KeyError:
                continue
            if self_as_dict[key] != "None":
                if "\\" in self_as_dict[key]:
                    path_type = PureWindowsPath
                elif "/" in self_as_dict[key]:
                    path_type = PurePosixPath
                else:
                    print(type(self_as_dict[key]), repr(self_as_dict[key]))
                    raise TypeError(
                        f"Could not detect whether {self_as_dict[key]} \
                        was Windows or Posix path"
                    )
                path_parts_list = path_type(self_as_dict[key]).parts
                for i, part in enumerate(path_parts_list):
                    if part == DATA_DIR.name:
                        self_as_dict[key] = str(
                            DATA_DIR.joinpath(*path_parts_list[i + 1 :])
                        )
                        # pycharm complains about using "i" after the loop.
                        break
                else:
                    print("Could not convert PureWindowsPath to local path!")
        # Probably not used, but might as well save the correct path here too
        self_as_dict["measurement_dir"] = str(MEASUREMENT_DIR)
        return cls(**self_as_dict)

    @classmethod
    def open(cls, m_id, measurement_dir=MEASUREMENT_DIR):
        """Opens the measurement given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(measurement_dir).iterdir()
                if path.stem.startswith(f"m{m_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no measurement with id = m{m_id}")
        return cls.load(path_to_file)

    def __repr__(self):
        return self.make_name()

    def make_name(self):
        """make a name for self from its id, sample, category, date, and technique"""
        category_string = "u"  # u for uncategorized
        if isinstance(self.category, (list, tuple)):
            category_string = self.category[0]
            for cat in self.category[1:]:
                category_string += f" and {cat}"
        elif self.category:
            category_string = self.category
        self.name = (
            f"m{self.id} is {self.sample_name} {category_string} on {self.date}"
            f" by {self.technique}"
        )
        return self.name

    @property
    def sample(self):
        if not self.sample_name:
            return
        from .sample import Sample

        try:
            return Sample.open(self.sample_name)
        except FileNotFoundError:
            return Sample(name=self.sample_name)

    @property
    def meas(self):
        """The ixdat ECMSMeasurement associated with the measurement"""
        if not self._meas:
            self.load_data()
        return self._meas

    @property
    def elog(self):
        if not self._elog:
            self.open_elog()
        return self._elog

    @property
    def tstamp(self):
        return self.meas.tstamp

    def __gt__(self, other):
        return self.tstamp > other.tstamp

    def __ge__(self, other):
        return self.tstamp >= other.tstamp

    def load_data(self):
        """load the ixdat meas from the EC_MS pkl file"""
        data_path = str(Path(self.raw_data_path).expanduser())
        self._meas = Meas.read(data_path, reader="EC_MS")
        if not self._meas.series_list:
            raise IOError(f"Dataset in {self.raw_data_path} loaded empty.")
        return self._meas

    def export_data(self):
        """SAVE the meas in the new data directory"""
        name = self.name if self.name else self.make_name()
        path_to_pkl = self.exported_data_path / (name + ".csv")
        self.meas.export(file_name=path_to_pkl)
        self.copied_at = time.time()

    def plot(self, *args, **kwargs):
        """shortcut to self.meas.plot"""
        kwargs.update(legend=False)
        return self.meas.plot_measurement(*args, **kwargs)

    def cut_meas(self, *args, **kwargs):
        self._meas = self.meas.cut(*args, **kwargs)

    def open_elog(self):
        from .elog import ElogEntry

        try:
            self._elog = ElogEntry.open(self.elog_number)
        except FileNotFoundError as e:
            print(f"'{self.name}' has no elog due to: {e}!")
            return

    def print_notes(self):
        if not self.elog:
            print(f"'{self}' has no elog.")
        else:
            print(f"\n### printing from '{self.elog}' ###\n")
            print(f"elog.field_data = {self.elog.field_data}")
            print(f"\n######## start of elog notes for '{self}' ###########\n")
            notes = self.elog.notes
            if self.EC_tag:
                try:
                    EC_tag_match = re.search(fr"\n\s+{self.EC_tag}", notes)
                except TypeError:
                    print(f"[problem searching for '{self.EC_tag}' in:\n{notes}]\n")
                    return
                # ^ note, EC_tag has the "..." already in it.
                if EC_tag_match:
                    notes = (
                        notes[0 : EC_tag_match.start()]
                        + "\n# ======================================== #\n"
                        + "# ===   MEASUREMENT NOTES START HERE   === #\n"
                        + "# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv #"
                        + notes[EC_tag_match.start() :]
                    )
            print(notes)
            print(f"\n\n######## end of elog notes for '{self}' ###########\n")

    @property
    def RE_vs_RHE(self):
        try:
            RE_vs_RHE_str = self.elog.field_data["RE_vs_RHE"]
            if RE_vs_RHE_str == "":
                raise KeyError
        except AttributeError:
            print(r"WARNING!!! Measurement '{self}' has no elog :(")
            RE_vs_RHE = None
        except (KeyError, TypeError):
            print(f"WARNING!!! No RE_vs_RHE in ({self.elog}), the elog for '{self}'")
            RE_vs_RHE = None
        else:
            RE_vs_RHE = float(re.search(FLOAT_MATCH, RE_vs_RHE_str).group())
        return RE_vs_RHE

    @property
    def A_el(self):
        return STANDARD_ELECTRODE_AREA

    @property
    def R_Ohm(self):
        try:
            R_ohm_str = self.elog.field_data["Resistor"]
            if R_ohm_str == "":
                raise KeyError
        except AttributeError:
            print(f"WARNING!!! Measurement '{self}' has no elog :(")
            R_ohm = None
        except (KeyError, TypeError):
            print(f"WARNING!!! No R_ohm in ({self.elog}), the elog for '{self}'")
            R_ohm = None
        else:
            R_ohm = float(re.search(FLOAT_MATCH, R_ohm_str).group())
        return R_ohm

    def get_icpms_points(self):
        """Return a list of ICPMSPoints from the measurement"""
        from .icpms import all_icpms_points

        ips = []  # icpms points
        ts = []  # sampling times, for sorting
        for ip in all_icpms_points():
            if ip.m_id == self.id and "duplicate" not in ip.description:
                ips += [ip]
                ts += [ip.sampling_time]

        if len(ips) > 1:
            # sort them by sampling time:
            ts, indeces = zip(*sorted((t_i, i) for i, t_i in enumerate(ts)))
            ips = [ips[index] for index in indeces]
        return ips

    def get_standard_experiment(self):
        from .experiment import all_standard_experiments

        for se in all_standard_experiments():
            if se.m_id == self.id:
                return se
        print(f"'{self}' is not a standard experiment.")
        return None
