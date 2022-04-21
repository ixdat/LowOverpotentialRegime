"""ISS handler module for pyOER.

Simple usage [deprecated]:
You have ISS of samples, "Reshma1" and "Reshma2". You can load all
these samples by loading "Reshma" without a postfix. The following
piece of code will load ISS experiments for both sample series,
create a plot of the isotopic oxygen ratios for every spectrum, and
opening a plot verifying how accurate the peak decomposition is.

---- Code begin ----
import pyOER

# Load samples
experiment_chain = pyOER.ISS('Reshma')

# Plot isotopic oxygen ratios
experiment_chain.plot_fit_ratios(True)

# Spectrum number 7 appears to be the only outlier, so compare
# spectrum 6 and spectrum 7:
experiment_chain.plot_fit(6)
experiment_chain.plot_fit(7)

# The 5% O-18 could be explained by improper background subtraction
# and is therefore seemingly within the fitting error.
---- Code end ----
"""
import json
import pickle5 as pickle
import pathlib
import datetime

import numpy as np
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from .tools import weighted_smooth as smooth

# from .tools import smooth
from .tools import get_range, dict_from_json
from .settings import DATA_DIR

from collections.abc import Mapping


def deep_update(d1, d2):
    """https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth"""
    if all((isinstance(d, Mapping) for d in (d1, d2))):
        for k, v in d2.items():
            d1[k] = deep_update(d1.get(k), v)
        return d1
    return d2


class ISSIterator:
    """Iterator class for ISS"""

    def __init__(self, iss_handler):
        # Reference to main class
        self._handle = iss_handler
        self._initial_index = self._handle.active
        self._index = 0
        # Loop through sets sorted by (earliest) date
        self._items = self._handle.keys

    def __next__(self):
        """Return next dataset."""
        if self._index < len(self._handle):
            self._handle.active = self._items[self._index]
            self._index += 1
            return self._handle.data
        # Restore active dataset
        self._handle.active = self._initial_index
        raise StopIteration


class ISS:
    """ISS handler"""

    def __init__(self, sample=None, fit=None, verbose=False):
        """Main interface for the ISS data"""
        self.verbose = verbose
        self.json_path = (
            pathlib.Path(__file__).absolute().parent.parent.parent / "tables" / "leis"
        )
        self.data_path = DATA_DIR / "Data" / "ISS" / "organized_pickles"
        self.extras_path = DATA_DIR / "Data" / "ISS" / "pickled_pickles"
        self._active = None
        self._relative = None
        self._set_active = None
        self._old_active = None
        self.plot_list = []
        self.sample = sample
        self._meta = None
        # All reference data
        with open(self.data_path / "iss_reference.pickle", "rb") as f:
            self._ref = pickle.load(f)
        self._init_extras()  # _background, _shifted
        self.fit_ratios = {}
        self.fit_coeffs = {}
        if sample is not None:
            self.get_sample(sample)
        if fit is not None:
            self.fit_with_reference(peaks=fit, plot=False)

    def all_leis(self):
        """Return a list of all leis samples."""
        all_samples = [sample.stem for sample in self.json_path.rglob("*.json")]
        all_samples.sort()
        return all_samples

    def __iter__(self):
        """Loop through the datasets sorted from earliest to last."""
        return ISSIterator(self)

    def __len__(self):
        """Return the number of datasets currently loaded for sample."""
        return len(self.keys)

    @property
    def data(self):
        """Return the active dataset."""
        return self._active[self.active]

    @property
    def fit_coeff(self):
        """Return the fit coefficients for active dataset.

        Returns a dict with the available keys."""
        return self.fit_coeffs[self.active]

    @property
    def fit_ratio(self):
        """Return the fitted (isotope) ratios for active dataset.

        Returns a dict with the available keys."""
        return self.fit_ratios[self.active]

    @property
    def keys(self):
        if self._active is None:
            print("Use get_sample(name) to select datasets")
            return
        selection = self.relative_to(datetime.datetime(9999, 1, 1))["before"]
        return [key for (key, date) in selection]

    @property
    def labels(self, print_=True):
        if self._active is None:
            print("Use get_sample(name) to select datasets")
            return
        if print_ is True:
            for key in self.keys:
                print(f"Key: {key}")
                print(f"\tSample: {self._active[key].sample}")
                print(f"\tComment (filename): {self._active[key].comment}")
                print(f"\tComment (from avg): {self._active[key].note[0]}")
                print(f"\tRecorded: {self._active[key].date}")
                print(f"\tPath: {self._active[key].filename}")

    def get_metadata(self, sample=None):
        """Fetch the JSON metadata corresponding to ´sample´.
        If ´None´ (default), fetch JSON corresponding to current sample."""
        if sample is None:
            if self.sample is None:
                raise ValueError("You have to choose a sample")
            sample = self.sample
        # TODO: will give an error if not matching an existing file
        with open(self.json_path / (str(sample) + ".json"), "r") as f:
            metadata = dict_from_json(json.load(f))
        if sample == self.sample:
            self._meta = metadata

        return metadata

    def save_json(self, metadata=None):
        """Save a metadata dictionary as JSON."""
        if metadata is None:
            if self._meta is None:
                raise ValueError("You have to choose a sample")
            metadata = self._meta
        if not isinstance(metadata, dict):
            # Not checking the contents/structure
            raise TypeError("´metadata´ in ´self.save_json´ must be of type dict")
        with open(self.json_path / metadata["file"], "w") as f:
            json.dump(metadata, f, indent=4)

    def show_meta_keys(self):
        """Print and return the available metadata keys."""
        meta_keys = []

        print("Available metadata keys:\n")
        print(" file")
        meta_keys.append("file")
        print(" data")
        meta_keys.append("data")
        for i in self._meta["data"][self.active].keys():
            print(f"  - {i}")
            meta_keys.append(i)
        print(" custom")
        meta_keys.append("custom")
        print("  - results")
        meta_keys.append("results")
        for i in self._meta["custom"]["results"][self.active].keys():
            print(f"     - {i}")
            meta_keys.append(i)
        print("  - measurements")
        meta_keys.append("measurements")
        for i in self._meta["custom"]["measurements"].keys():
            print(f"     - {i}")
            meta_keys.append(i)
        print()
        return meta_keys

    def meta(self, key=None):
        """Serve the contents from metadata dict"""
        if key is None:
            return self.show_meta_keys()

        if key == "file":
            return self._meta["file"]

        if key == "data":
            return self._meta["data"][self.active]

        if key in self._meta["data"][self.active].keys():
            return self._meta["data"][self.active][key]

        if key == "custom":
            return self._meta["custom"]

        if key == "results":
            return self._meta["custom"]["results"][self.active]

        if key in self._meta["custom"]["results"][self.active].keys():
            return self._meta["custom"]["results"][self.active][key][
                self.data.default_scan
            ]

        if key == "measurements":
            return self._meta["custom"]["measurements"]

        if key in self._meta["custom"]["measurements"].keys():
            return self._meta["custom"]["measurements"][key]

    def update_meta(self, key, value):
        """Update a field in the metadata dict."""
        # TODO: make more robust
        if (
            key == "file"
            or key == "data"
            or key in self._meta["data"][self.active].keys()
        ):
            raise KeyError(
                f"The data for {key} is generated from raw data and shouldn't be"
                'changed. Use the "custom" data fields instead!'
            )
        if key is None:
            key = "None"

        if key == "custom":
            if not isinstance(value, dict):
                raise TypeError("´value´ must be of type ´dict´ when used here")
            dictionary = self._meta["custom"]
            dictionary = deep_update(dictionary, value)

        elif key == "results":
            if not isinstance(value, dict):
                raise TypeError("´value´ must be of type ´dict´ when used here")
            dictionary = self._meta["custom"]["results"]
            dictionary = deep_update(dictionary, value)

        elif key == "measurements":
            if not isinstance(value, dict):
                raise TypeError("´value´ must be of type ´dict´ when used here")
            dictionary = self._meta["custom"]["measurements"]
            dictionary = deep_update(dictionary, value)

        elif key.startswith("m") and key[1:].isnumeric():
            self._meta["custom"]["measurements"][key] = value

        elif key in self._meta["custom"]["results"][self.active].keys():
            self._meta["custom"]["results"][self.active][key][
                self.data.default_scan
            ] = value

        else:
            self.show_meta_keys()
            raise KeyError(f'Key "{key}" does not match the structure of the metadata')

    @property
    def active(self):
        """Return active selection in a pretty way"""
        if self._set_active is None:
            raise ValueError("Use ´get_sample´ to select a dataset")
        else:  # return selected dataset
            return self._set_active

    @active.setter
    def active(self, value):
        """Set active selection to a key or list of keys.
        Set to None for certain effect."""
        ### TODO
        # do some check of value here
        ###
        self._set_active = value

    def relative_to(self, timestamp):
        """Take active list and return sorted relative to timestamp
        'timestamp' must be either datetime object or string of type 20A31"""
        a_to_1 = {
            "a": 1,
            "b": 2,
            "c": 3,
            "d": 4,
            "e": 5,
            "f": 6,
            "g": 7,
            "h": 8,
            "i": 9,
            "j": 10,
            "k": 11,
            "l": 12,
        }
        if isinstance(timestamp, str):
            if len(timestamp) == 5:
                year = 2000 + int(timestamp[:2])
                month = a_to_1[timestamp[2].lower()]
                day = int(timestamp[-2:])
                timestamp = datetime.datetime(year, month, day)
        elif isinstance(timestamp, type(datetime.datetime.now())):
            pass
        else:
            print("Timestamp type not understood")
            return
        list_ = [(key, self._active[key].date) for key in self._active.keys()]
        list_.sort(key=lambda x: x[1])
        match = False
        for i, (key, date) in enumerate(list_):
            if timestamp < date:
                match = True
                break
        if not match:
            i += 1
        return {
            "match": timestamp,
            "before": list_[:i],
            "after": list_[i:],
        }

    def get_sample(self, sample):
        """Get all ISS data involving sample_name"""
        self.sample = sample
        try:
            self.get_metadata()
        except ValueError:
            return
        keys = [key for key in self._meta["data"]]
        filenames = [self._meta["data"][key]["pickle_name"] for key in keys]
        self._active = self._load_set(filenames, keys)
        self._init_extras()
        self._load_extras()
        self.active = 0

    def _load_set(self, filenames, keys=None):
        """Take list of filenames and load it into dictionary"""
        iss_dict = dict()
        if keys is None:
            iterator = enumerate(filenames)
        else:
            iterator = list(zip(keys, filenames))
        for i, filename in iterator:
            with open(self.data_path / filename, "rb") as f:
                iss_dict[i] = pickle.load(f)
            if self.verbose:
                print(i, filename)
        return iss_dict

    def _load_extras(self):
        """Load the extra information which is calculated from the raw data."""
        path_to_file = self.extras_path / (self.sample + ".pickle")
        try:
            with open(path_to_file, "rb") as f:
                data = pickle.load(f)
        except IOError:
            print("File not found error:", path_to_file)
            data = None
            self._init_extras()
            return

        # Background data
        for i, region in enumerate(data["background"]):
            self._background[region] = {}
            for key in data["background"][region]:
                self._background[region][key] = data["background"][region][key]

        # Alignment data
        for j, region in enumerate(data["shifted"]):
            self._shifted[region] = {}
            for key in data["shifted"][region]:
                self._shifted[region][key] = data["shifted"][region][key]

        # Update region if unambiguous
        if i == 0 and j == 0:
            self._region = region

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, value):
        if value in self._shifted or value in self._background:
            self._region = value
        else:
            raise ValueError(f"Region ´{value}´ does not exist.")

    def aligned(self, key=None):
        """Alias for ´shifted´."""
        return self.shifted(key)

    def shifted(self, key=None):
        """Return the shifted spectrum of the active dataset according to region."""
        if not self.region in self._shifted:
            print("\nAvailable regions for shifted (aligned) data:")
            for region in self._shifted:
                print(f" - {region}")
            print("Use ´self.region´ to select/activate a region.")
            return

        selector = (self.active, self.data.default_scan)
        if key in self._shifted[self.region][selector]:
            return self._shifted[self.region][selector][key]
        print("\nAvailable keys for shifted (aligned) data:")
        for key in self._shifted[self.region][selector]:
            print(f" - {key}")

    def background(self, key=None):
        """Return the background spectrum of the active dataset according to region."""
        if not self.region in self._background:
            print("\nAvailable regions for background subtracted data:")
            for region in self._background:
                print(f" - {region}")
            print("Use ´self.region´ to select/activate a region.")
            return

        selector = (self.active, self.data.default_scan)
        if key in self._background[self.region][selector]:
            return self._background[self.region][selector][key]
        print("\nAvailable keys for background subtracted data:")
        for key in self._background[self.region][selector]:
            print(f" - {key}")

    def save_extras(self):
        """Save the current version of extra data to file."""
        data = {
            "background": self._background,
            "shifted": self._shifted,
        }
        destination = self.extras_path / (self.sample + ".pickle")
        old_data = None
        if destination.exists():
            with open(destination, "rb") as f:
                old_data = pickle.load(f)
        change = False
        if self.verbose:
            print("\n", change)
        try:
            regions = np.unique(
                list(
                    zip(
                        data["background"].keys(),
                        old_data["background"].keys(),
                    )
                )
            )
            for region in regions:
                index = []
                for i in list(data["background"][region].keys()) + list(
                    old_data["background"][region].keys()
                ):
                    if not i in index:
                        index.append(i)
                for i, j in index:
                    if self.verbose:
                        print("\n", region, i, j)
                    check_1 = np.all(
                        data["background"][region][i, j]["x"]
                        == old_data["background"][region][i, j]["x"]
                    )
                    change = change or not check_1
                    if self.verbose:
                        print(f"\n  1 {check_1}")
                        print("  back x: ", change)

                    new = data["background"][region][i, j]["y"]
                    old = old_data["background"][region][i, j]["y"]
                    check_2 = np.all(new[np.isfinite(new)] == old[np.isfinite(old)])
                    change = change or not check_2
                    if self.verbose:
                        print(f"\n  2 {check_2}")
                        print("  back y: ", change)

                    check_3 = (
                        data["background"][region][i, j]["limits"]
                        == old_data["background"][region][i, j]["limits"]
                    )
                    change = change or not check_3
                    if self.verbose:
                        print(f"\n  3 {check_3}")
                        print("  back limits: ", change)

                    new = data["shifted"][region][i, j]["xy"]
                    old = old_data["shifted"][region][i, j]["xy"]
                    check_4 = np.all(new[np.isfinite(new)] == old[np.isfinite(old)])
                    change = change or not check_4
                    if self.verbose:
                        print(f"\n  4 {check_4}")
                        print("  shift xy: ", change)

                    check_5 = [
                        data["shifted"][region][i, j][key]
                        == old_data["shifted"][region][i, j][key]
                        for key in ["region", "masses", "limits", "good"]
                    ]
                    change = change or not all(check_5)
                    if self.verbose:
                        print(f"\n  5 {check_5}")
                        print("  shift keys: ", change)
                        print()

        except (KeyError, TypeError):
            change = True
        if change:
            print(f"\n{destination.stem}: ", end="")
            print("OVERWRITING: Extra data has changed from file.\n")
            with open(destination, "wb") as f:
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
        else:
            if self.verbose:
                print(f"\n{destination.stem}: ", end="")
                print("Extra data equivalent to file - not overwriting...\n")

    def _init_extras(self):
        self._region = None
        self._background = {}
        self._shifted = {}

    def plot(self, selection=[], mass_lines=[], show=True):
        """Plot selected spectra.
        'selection' must be a list of keys matching 'iss_dict' returned by
        self.load_set.
        """
        self.backup_active()
        if len(selection) == 0:
            # Plot all (sorted)
            selection = self.keys
        plt.figure(f"Raw leis data for: {self.sample}")

        if len(mass_lines) > 0:
            self.data.add_mass_lines(
                mass_lines,
                color="gray",
                labels=False,
            )
        if isinstance(selection, int):
            selection = [selection]
        for key in selection:
            self.active = key
            if self.verbose:
                print(f"\nPlotting {self.sample}: {key}")
                print(f"{self.data.date}    ({self.data.setup})")
                print(f"Comment: {self.data.comment}")
                print(f"Note(s):")
                if isinstance(self.data.note, list):
                    for note in self.data.note:
                        print(f" - {note}")
                elif isinstance(self.data.note, dict):
                    for index, note in self.data.note.items():
                        print(f" - {index}: {note}")
                else:
                    print(f" - {self.data.note}")
            for i in self.data:
                plt.plot(
                    self.data.x,
                    self.data.y,
                    label=f"{key} - {self.sample} ({i})",
                )

        plt.legend()
        plt.xlabel("Energy (eV)")
        plt.ylabel("Counts per second")
        self.restore_active()
        if show is True:
            plt.show()

    def show(self):
        """Show figures"""
        plt.show()

    def backup_active(self):
        self._old_active = self.active

    def restore_active(self):
        self.active = self._old_active

    # Maybe a @classmethod ?
    def fit_with_reference(
        self,
        selection=None,
        peaks=[],
        plot_result=False,
        region_names=["oxygen"],
        align=True,
        recalculate=False,
    ):
        """Fit reference data to selected datasets.

        Input  :
            selection (list): where each element is an integer representing a
                key/index in self.keys. Each of these ISS objects will be fit
                to the  reference data. Defaults to all datasets.
            peaks (list): where each element is an integer representing the
                mass of the elements to fit. If elements should be compared
                group-wise, as in O-16 or O-18 to total oxygen signal, nest
                them in a list of their own, e.g. peaks=[[16, 18], 40]. No
                sensitivity factors are applied.
            plot_result (bool): whether or not to autoplot results. Default False.
            align (bool): has to be True (default) for now.

        Output :
            ratios (dict): which will contain the ratio of peak-to-peak(s)
                for every combination of peak in ´peaks´.

        Example:
        >>> ratio = self.fit_with_reference(peaks=[[16, 18]])
        >>> print(ratio['16/18'])
        2.917
        >>> print(ratio['16'])
        0.7447
            __________________________________________
        Notes:

        REFERENCE DATA:
            self._ref (nested dict)
            self._ref[setup][peak].keys() =

                'xy': xy-data of original spectrum
                'background': peak-subtracted original spectrum
                'peak': background-subtracted original spectrum
                'area': integrated signal of 'peak'
                'region': the [low, high] delimiter used for the peak
                     identification/background subtraction
                'file': full path to original spectrum on host computer
                'iss': original ISS object

        FITTING METHOD:
            for each spectrum in ´selection´:
                for each peak in ´peaks´:
                    subtract background from spectrum using same region as ref;
                    for each nested peak (if any):
                        add scaled nested peaks to background subtracted data for best fit;
                    save fit result to ´results´ dictionary;

        RETURN METHOD:
            for each result in ´results´:
                for each other result in ´results´:
                    if result != other result:
                        save to ´ratios´: (result)/(other result)
            for each nested peak:
                save to ´ratios´:
                    'peak1' = peak1 / (peak1 + peak2)
                    'peak2' = peak2 / (peak1 + peak2)
        """

        # Main loop
        coeffs = {}
        ratios = {}
        for data in self:
            if selection:
                if self.active not in selection:
                    continue
            if plot_result:
                plt.figure(f"Fitting: {data.sample} - {self.active}")

            # Define region limits by comparing with reference file
            ref = self._ref[data.setup]
            for peak, region_name in list(zip(peaks, region_names)):
                # Peak is single
                if isinstance(peak, int):
                    region = ref[peak]["region"]
                    N = 1
                # Peak is a group
                elif isinstance(peak, list):
                    region = ref[peak[0]]["region"]
                    N = len(peak)
                    for _ in peak:
                        if ref[_]["region"] != region:
                            raise ValueError(
                                (
                                    f'Grouped peaks "{peak}" are not defined by the '
                                    "same region"
                                )
                            )
                else:
                    raise TypeError(
                        (
                            f"Item in kwarg peaks is not understood: {peak}\n"
                            "Must be an integer or list of integers."
                        )
                    )

                if self.verbose:
                    print("Selected: ", self.active)
                    print("Region name: ", region_name)
                    print("Region: ", region)

                if not region_name in self._shifted:
                    self._shifted[region_name] = {}
                    if region_name != "oxygen":
                        raise NotImplementedError
                if not self.active in self._shifted[region_name] or recalculate:
                    try:
                        aligned_data = align_spectra(
                            [self.data],
                            limits=[350, 520],
                            masses=[16, 18],
                            key=region_name,
                            plot_result=plot_result,
                            verbose=self.verbose,
                        )
                    except ValueError:
                        aligned_data = {"good": False}
                    for scan, value in aligned_data.items():
                        self._shifted[region_name][(self.active, scan)] = value
                else:
                    print("Not running ´align_spectra´...")
                self.region = region_name

                # Subtract background and make accessible afterwards
                if self.active == 0:
                    if self.verbose:
                        print(self.data.sample)
                if not region_name in self._background:
                    self._background[region_name] = {}
                for scan in data:
                    if not self.shifted("good"):
                        if self.verbose:
                            print(
                                f"Skipping bad data.. {data.sample} ({self.active}, {scan})"
                            )
                        results = {
                            self.active: {
                                "O16": {scan: None},
                                "O18": {scan: None},
                                "c_16": {scan: 0},
                                "c_18": {scan: 0},
                            }
                        }
                        self.update_meta("results", results)
                        self._background[region_name][(self.active, scan)] = {
                            "x": [],
                            "y": [],
                            "limits": [region],
                        }
                        continue  # skip bad data set
                    coeffs[(self.active, scan)] = {}
                    ratios[(self.active, scan)] = {}
                    if self.verbose:
                        print("Good: ", (self.active, scan), self.shifted("good"))
                    if (
                        not (self.active, scan) in self._background[region_name]
                        or recalculate
                    ):
                        if self.shifted("good"):
                            background = subtract_single_background(
                                self.shifted("xy"),
                                ranges=[region],
                            )
                            self._background[region_name][(self.active, scan)] = {
                                "x": self.shifted("x"),
                                "y": background,
                                "limits": [region],
                            }
                    else:
                        if self.verbose:
                            print(
                                "Background subtraction already performed. Skipping..."
                            )

                    isolated_peak = data.y - self.background("y")
                    isolated_peak[np.isnan(isolated_peak)] = 0.1
                    if plot_result:
                        plt.plot(
                            data.x,
                            self.background("y"),
                            "k:",
                            # label='Background',
                        )
                        plt.plot(
                            data.x,
                            data.y,
                            "k-",
                            label="Raw data",
                        )
                        plt.plot(
                            self.background("x"),
                            self.background("y"),
                            "b:",
                            # label='Background',
                        )
                        plt.plot(
                            self.shifted("x"),
                            data.y,
                            "b-",
                            label="Aligned data",
                        )

                    # Create a common x-axis for comparisons
                    pseudo_x = np.linspace(
                        region[0],
                        region[1],
                        (region[1] - region[0]) * 10 + 1,
                    )
                    interp_dat = interp1d(
                        self.shifted("x"),
                        isolated_peak,
                        kind="linear",
                    )
                    interp_back = interp1d(
                        self.background("x"),
                        self.background("y"),
                        kind="linear",
                    )
                    interp_ref = {}
                    interp_ref[16] = interp1d(
                        ref[16]["x"],
                        ref[16]["peak"],
                        kind="linear",
                    )
                    interp_ref[18] = interp1d(
                        ref[18]["x"],
                        ref[18]["peak"],
                        kind="linear",
                    )
                    mask = get_range(pseudo_x, *region)
                    dat_x = pseudo_x[mask]
                    dat_y = interp_dat(dat_x)
                    if plot_result:
                        plt.plot(
                            ref[16]["x"],
                            ref[16]["peak"],
                            "r:",
                            label="O16 ref",
                        )
                        plt.plot(
                            ref[18]["x"],
                            ref[18]["peak"],
                            "g:",
                            label="O18 ref",
                        )

                    def func(x, *args):
                        """Fitting function"""
                        signal = x * 0
                        for arg, i in list(zip(args, peak)):
                            signal += arg * interp_ref[i](x)
                        return signal

                    # Fit reference to data
                    fit, _ = curve_fit(
                        func,
                        dat_x,
                        dat_y,
                        p0=[2.0] * N,
                        bounds=(0, 3),
                    )
                    fitted_signal = interp_back(dat_x)
                    for i in range(len(peak)):
                        coeffs[(self.active, scan)][peak[i]] = fit[i]
                        fitted_signal += interp_ref[peak[i]](dat_x) * fit[i]

                    if plot_result:
                        plt.plot(
                            dat_x,
                            fitted_signal,
                            "y-",
                            label="Best fit",
                        )
                    # Calculate output ratios
                    total = 0
                    all_peaks = []
                    for peak in peaks:
                        if isinstance(peak, list):
                            for peak_ in peak:
                                all_peaks.append(peak_)
                        else:
                            all_peaks.append(peak)

                    """
                    for peak1 in all_peaks:
                        for peak2 in all_peaks:
                            if peak1 == peak2:
                                continue
                            if self.shifted('good'):
                                ratios[(self.active, scan)][f'{peak1}/{peak2}'] = (
                                    coeffs[(self.active, scan)][peak1]
                                    / coeffs[(self.active, scan)][peak2]
                                    * ref[peak1]['area']
                                    / ref[peak2]['area']
                                    )
                            else:
                                ratios[(self.active, scan)][f'{peak1}/{peak2}'] = None
                    """

                    # Group ratios
                    for peak in peaks:
                        if not isinstance(peak, list):
                            continue
                        total = 0
                        for peak_ in peak:
                            if self.shifted("good"):
                                ratios[(self.active, scan)][f"{peak_}"] = (
                                    ref[peak_]["area"]
                                    * coeffs[(self.active, scan)][peak_]
                                )
                                total += ratios[(self.active, scan)][f"{peak_}"]
                            else:
                                ratios[(self.active, scan)][f"{peak_}"] = None

                        if self.shifted("good"):
                            for peak_ in peak:
                                ratios[(self.active, scan)][f"{peak_}"] /= total

                    if plot_result:
                        data.add_mass_lines(all_peaks)
                        plt.legend()
                    # Save in object
                    self.fit_ratios[(self.active, scan)] = ratios[(self.active, scan)]
                    self.fit_coeffs[(self.active, scan)] = coeffs[(self.active, scan)]
                    if self.shifted("good"):
                        results = {
                            self.active: {
                                "O16": {scan: ratios[(self.active, scan)]["16"] * 100},
                                "O18": {scan: ratios[(self.active, scan)]["18"] * 100},
                                "c_16": {scan: coeffs[(self.active, scan)][16]},
                                "c_18": {scan: coeffs[(self.active, scan)][18]},
                            }
                        }
                        self.update_meta("results", results)

        return ratios, coeffs

    def plot_fit_ratios(self, show_plot=False):
        """Make a plot of O16/18 ratios for instance"""
        # Make sure references have been fitted
        if len(self.fit_ratios.keys()) == 0:
            if self.verbose:
                print('Calling method "fit_with_reference(peaks=[[16, 18]])')
            self.fit_with_reference(peaks=[[16, 18]])

        # Prepare plot
        fig = plt.figure("Fit ratios plot title")
        ax = fig.add_axes([0.05, 0.15, 0.9, 0.6])
        colors = ["k", "r", "g", "b", "m"] * 10

        # Plot all O-16 ratios
        plot_data = []
        counter = 0

        for _ in self:
            i = self.active
            # Skip bad data
            if not self.shifted("good"):
                counter += 1
                continue
            # Plot good data
            plt.plot(
                counter,
                self.fit_ratios[i]["16"] * 100,
                "o",
                color=colors[0],
            )
            plot_data.append(
                [
                    self.data.sample,
                    self.active,
                    self.data.sample,
                    self.data.date,
                    counter,
                    self.fit_ratios[i]["16"],
                    self.fit_ratios[i]["18"],
                ]
            )
            counter += 1

        # Plot formatting
        xticks = [i for (gen_name, data_object, name, date, i, r1, r2) in plot_data]
        dates = [
            date_formatter(date)
            for (gen_name, data_object, name, date, i, r1, r2) in plot_data
        ]
        xlabels = [
            f"{gen_name} {name.lstrip(gen_name)} - {active}"
            for (gen_name, active, name, date, i, r1, r2) in plot_data
        ]

        # Some of the following secondary axis methods requires matplotlib > 3.1.x
        secaxx = ax.secondary_xaxis("top")
        secaxy = ax.secondary_yaxis("right")

        # Update canvas
        fig.canvas.draw()

        secaxy.set_ylabel("O-18 ratio (%)")
        yticks = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)
        secaxy.set_ticks(yticks)
        yticks.reverse()
        yticks = [str(i) for i in yticks]
        secaxy.set_yticklabels(yticks)
        secaxx.set_xticks(xticks)
        secaxx.set_xticklabels(dates, rotation=90, fontsize=12)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=90, fontsize=12)
        ax.set_ylabel("O-16 ratio (%)")
        plt.grid(True)
        if show_plot is True:
            plt.show()

    def plot_fit(self, index=None, labels=True, show=True):
        """Visually verify the automatic fit to reference data"""

        # Temporarily change the active dataset
        self.backup_active()
        if index is not None:
            self.active = index

        # Make sure references have been fitted
        if len(self.meta("results").keys()) == 0:
            if self.verbose:
                print('Calling method "fit_with_reference(peaks=[[16, 18]])')
            self.fit_with_reference(peaks=[[16, 18]])

        # Initialize figure
        plt.figure(
            f"Peak Deconvolution _ {self.sample} - {self.active, self.data.default_scan}"
        )

        # Compared x arrays are shifted with respect to each other.
        x_common = np.linspace(0, 1000, num=1001)
        setup = self.data.setup
        ref1 = self._ref[setup][16]
        ref2 = self._ref[setup][18]

        # Raw + background
        plt.plot(
            self.data.x,
            self.data.y,
            "m-",
            label="Raw unaligned",
        )
        plt.plot(
            self.shifted("x"),
            self.data.y,
            "k-",
            label="Raw aligned",
        )
        plt.plot(
            x_common,
            get_common_y(
                x_common,
                self.shifted("x"),
                self.shifted("y"),
            ),
            "b-",
            label="Aligned+smoothed",
        )
        background = get_common_y(
            x_common,
            self.background("x"),
            self.background("y"),
        )
        plt.plot(
            x_common,
            background,
            "b:",
            label="Background",
        )

        y_ref1 = get_common_y(
            x_common,
            ref1["x"],
            ref1["peak"] * self.meta("c_16"),
        )
        y_ref2 = get_common_y(
            x_common,
            ref2["x"],
            ref2["peak"] * self.meta("c_18"),
        )

        # Total fit
        plt.plot(
            x_common,
            background + y_ref1 + y_ref2,
            "y-",
            label="Sum of components",
        )
        plt.plot(
            x_common,
            y_ref1,
            "r-",
            label="O-16 component",
        )
        plt.plot(
            x_common,
            y_ref2,
            "g-",
            label="O-18 component",
        )
        self.data.add_mass_lines([16, 18, 101], labels=labels)

        # Show
        plt.title(f"{self.sample}: {self.active}")
        plt.xlabel("Energy (eV)")
        plt.ylabel("Counts per second")
        plt.xlim(300, 800)
        mask = get_range(self.shifted("x"), 300, 800)
        plt.ylim(0, max(self.shifted("y")[mask]))
        plt.legend()
        if show:
            plt.show()

        # Change back the active dataset
        self.restore_active()


def get_common_y(x_common, x, y, num=4):
    y_common = np.interp(x_common, x, smooth(y, num), left=0, right=0)
    return y_common


def date_formatter(date, latex=True):
    """Take datetime object and return string of YYADD"""
    YY = date.year - 2000
    M = date.month
    DD = date.day
    hh = date.hour
    mm = date.minute
    ss = date.second
    translate = {
        1: "a",
        2: "b",
        3: "c",
        4: "d",
        5: "e",
        6: "f",
        7: "g",
        8: "h",
        9: "i",
        10: "j",
        11: "k",
        12: "l",
    }
    string = f"{YY}{translate[M].upper()}{DD} {hh}:{mm}:{ss}"
    string = (
        (r"$\bf{" if latex else "")
        + f"{str(YY).zfill(2)}{translate[M].upper()}{str(DD).zfill(2)}"
        + (r"}$" if latex else "")
        + f"   {str(hh).zfill(2)}:{str(mm).zfill(2)}:{str(ss).zfill(2)}"
    )
    return string


def subtract_single_background(xy, ranges=[], avg=3, verbose=False):
    """Subtract the background from a single spectrum"""
    x = xy[:, 0]
    y = xy[:, 1]
    background = np.copy(y)
    for limit in ranges:
        indice = get_range(x, *limit)
        # if first index is chosen
        # OR
        # if last ten indice are included
        if indice[0] == 0 or indice[-1] > len(x) - 10:
            if verbose:
                print("Uhh", indice[0], indice[-1], limit)
                print(f"Searching for indice within limits: {limit}")
                print(
                    f"First and last index: {indice[0]} and {indice[-1]} out of total {len(x) - 1}"
                )
                print(f"This is x = [{x[indice[0]]} and {x[indice[-1]]}]")
            background[indice] = 0
        elif len(indice) == 0:
            if verbose:
                print("Did not find data within limit: {}".format(limit))
        else:
            y1 = np.average(y[indice[0] - avg : indice[0] + avg])
            y2 = np.average(y[indice[-1] - avg : indice[-1] + avg])
            a_coeff = (y2 - y1) / (limit[1] - limit[0])
            b_coeff = y1 - a_coeff * limit[0]
            background[indice] = x[indice] * a_coeff + b_coeff
    return background


def align_spectra(
    iss_data,
    limits=[350, 520],
    masses=[16, 18],
    key="oxygen",
    plot_result=False,
    verbose=False,
    func_type="skewed",
):
    """Shift the iss data within 'limits' region to snap maximum signal
    unto nearest mass in list 'masses'.

        function (str): One of 'parabola', 'gauss' or 'skewed' (default). Determines the
        type of function used to align the spectra."""

    from scipy.optimize import curve_fit
    from scipy.special import erf

    if plot_result:
        import matplotlib.pyplot as plt

    def parabola(x, a, b, c):
        """2nd degree polynomial"""
        return a * x**2 + b * x + c

    def gauss(x, A, x0, sigma):
        """Gauss function or normal distribution"""
        return A * np.exp(-((x - x0) ** 2) / 2 / sigma**2)

    def skewed(x, A, x0, sigma, alpha):
        """Skewed gauss function"""
        return 2 / sigma * gauss(x, A, x0, sigma) * erf(alpha * (x - x0) / sigma)

    if verbose:
        print('Entering function "align_spectra"')
    if plot_result:
        old_ax = plt.gca()

    return_data = []
    for data in iss_data:
        # Initialize attributes
        shifted = {}
        smoothed = {}

        info = {}
        for scan in data:
            # Get index of region of interest
            index = get_range(data.x, *limits)
            # Find maximum in region
            ys = smooth(data.y, num=4)
            smoothed[key] = ys
            maximum = max(ys[index])
            if not np.isfinite(maximum):
                good = False
                # TODO: Add more information about data set
                print("´align_spectra´: no data found within set limits. Skipping...")
                info.update(
                    {
                        scan: {
                            "x": None,
                            "y": None,
                            "xy": None,
                            "region": key,
                            "masses": masses,
                            "limits": limits,
                            "good": good,
                        }
                    }
                )
                continue
            good = True
            i_max = np.where(ys == maximum)[0]
            x_max = data.x[i_max][0]

            # Estimate fitting parameters
            width = 20  # Estimate of peak width
            if func_type == "skewed":
                p0 = [maximum, x_max + 10, width, -1]
                function = skewed
            elif func_type == "parabola":
                a = -0.5 * maximum / width**2
                b = -2 * a * x_max
                c = maximum + b**2 / 4 / a
                p0 = [a, b, c]
                function = parabola
            elif func_type == "gauss":
                p0 = [maximum, x_max, width]
                function = gauss
            else:
                raise ValueError(f"func_type {func_type} not a valid option")
            new_index = get_range(data.x, x_max - 15, x_max + 15)
            fit, _ = curve_fit(
                function,
                data.x[new_index],
                data.y[new_index],
                p0=p0,
                maxfev=100000,
            )
            if verbose:
                print("Result of fit: ", fit)
            if plot_result:
                plt.figure(
                    (
                        f"Aligning {data.sample} {date_formatter(data.date, latex=False)}"
                        " - mass {masses}"
                    )
                )
                plt.plot(
                    data.x[index],
                    data.y[index],
                    "k-",
                    label="Raw data",
                )
                plt.plot(
                    data.x[new_index],
                    function(data.x[new_index], *p0),
                    "g-",
                    label="Estimated max",
                )
                plt.plot(
                    data.x[new_index],
                    function(data.x[new_index], *fit),
                    "r-",
                    label="Best fit max",
                )
            if function == parabola:
                new_x_max = -fit[1] / 2 / fit[0]
                if verbose:
                    print(f'Raw "maximum" x: {x_max}\nFitted x: {new_x_max}')
            elif function == gauss:
                new_x_max = fit[1]
                if verbose:
                    print(f'Raw "maximum" x: {x_max}\nFitted x: {new_x_max}')
            elif function == skewed:
                fit_x = np.linspace(min(data.x), max(data.x), num=16000)
                fit_y = skewed(fit_x, *fit)
                fit_i = np.where(fit_y == max(fit_y))[0]
                new_x_max = fit_x[fit_i][0]
                if verbose:
                    print(f'Raw "maximum" x: {x_max}\nFitted x: {new_x_max}')
            x_max = new_x_max

            # Find difference from reference
            energies = data.convert_energy(np.array(masses))
            distances = x_max - energies
            distance = distances[np.where(abs(distances) == min(abs(distances)))[0][0]]
            # If distance is too big, something is wrong with the algorithm
            if verbose:
                print(f"Distance between fitted maximum and expected: {distance} eV")
            max_distance = 30
            if abs(distance) > max_distance:
                msg = (
                    f"align_spectra algorithm tried to shift the spectrum {distance} eV"
                    f" which is more than the programmed limit: {max_distance} eV.\n"
                    "If this is okay, you need to change this limit."
                )
                distance = 0
                if verbose:
                    print("***\nDismissing alignment algorithm !\n***")
                    print(msg)
                info.update(
                    {
                        scan: {
                            "x": None,
                            "y": None,
                            "xy": None,
                            "region": key,
                            "masses": masses,
                            "limits": limits,
                            "good": False,
                        }
                    }
                )
                continue
                # raise ValueError(msg)

            # Snap to nearest line
            shifted[key] = data.x - distance
            if plot_result:
                plt.plot(
                    shifted[key],
                    data.y,
                    "b-",
                    label="Aligned raw data",
                )
                plt.plot(
                    shifted[key],
                    ys,
                    "c-",
                    label="Aligned and smoothed",
                )
                data.add_mass_lines(masses)
                plt.legend()
            info.update(
                {
                    scan: {
                        "x": shifted[key],
                        "y": smoothed[key],
                        "xy": np.vstack((shifted[key], smoothed[key])).T,
                        "region": key,
                        "masses": masses,
                        "limits": limits,
                        "good": good,
                    }
                }
            )

        # Return new data
        return_data.append(info)
    if plot_result:
        plt.sca(old_ax)

    # Don't return as a list if it only contains a single item
    if len(iss_data) == 1:
        return return_data[0]
    else:
        return return_data


class DataIterator:
    """Iterate through datasets in Data class."""

    def __init__(self, data):
        self._data = data
        self._initial_index = data.default_scan
        self._index = 0

    def __next__(self):
        if self._index < self._data.scans:
            self._data.default_scan = self._index
            self._index += 1
            return self._index - 1
        # Restore original dataset before stopping iteration
        self._data.default_scan = self._initial_index
        raise StopIteration


class Data:
    """Load an ISS experiment exported as text or VAMAS file.

    Class loader copied from github.com/Ejler/DataTreatment/ISS.py
    Renamed Experiment() -> Data()

    Author: Jakob Ejler Sorensen
    Version: 5.2
    Date: 2021 July 21
    """

    def __init__(self, filename, mass=4, theta=146.7, E0=1000, default_scan=0):
        """Initialize the class"""
        # Constants
        self.settings = dict()
        self.settings["mass"] = mass
        self.settings["theta"] = theta
        self.settings["E0"] = E0
        self.default_scan = default_scan

        # Initialize variables
        self.energy = dict()
        self.cps = dict()
        self.dwell = dict()
        self.mode = dict()
        self.mode_value = dict()
        self.note = dict()
        self.date = ""
        filename = str(filename)
        self.filename = filename

        # Convenience function variables
        self.peak_positions = None
        self.peak_heights_raw = None
        self.peak_heights_bg = None
        self._background = None
        self.background_settings = {
            "type": None,
            "ranges": None,
            "on": False,
        }

        # ----------------------------------------------------------------------
        # Read data from old VAMAS block file
        if filename.endswith(".vms"):
            # Open filename with ISS data
            f = open(filename, "r")
            lines = f.readlines()
            f.close()
            # Old format:
            if lines[6].lower().startswith("experiment type"):
                self.setup = "omicron"
                self.format = "Old VAMAS"
                # print('Loading file: ' + filename)
                blocks_4 = [
                    i
                    for i, line in enumerate(lines)
                    if (line.strip() == "-1")
                    and (lines[i + 1].lower().strip() == "kinetic energy")
                ]
                blocks_2_ISS = [
                    i
                    for i, line in enumerate(lines)
                    if (line.strip() == "ISS") and (lines[i + 1].strip() == "")
                ]
                print(lines[9].rstrip())
                self.scans = len(blocks_4)
                if len(blocks_4) == int(lines[9].rstrip()) and len(blocks_4) == len(
                    blocks_2_ISS
                ):
                    self.scans = len(blocks_4)
                else:
                    msg = 'Error: Identified {} "Block 4", {} "Block 2", but "Block 1" says: {}'
                    msg = msg.format(
                        len(blocks_4), len(blocks_2_ISS), int(lines[9].rstrip())
                    )
                    raise ImportError(msg)

                # Copy data points
                self.note = dict()
                for counter, block in enumerate(blocks_4):
                    if not len(lines[blocks_2_ISS[counter] - 1]) == 5:
                        self.note[counter] = lines[blocks_2_ISS[counter] - 1].rstrip()
                    else:
                        self.note[counter] = ""
                    self.mode[counter] = lines[block - 11].rstrip()
                    self.mode_value[counter] = float(lines[block - 10].rstrip())
                    self.dwell[counter] = float(lines[block + 9].rstrip())
                    data_points = int(lines[block + 16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[block + 4].rstrip())
                    E_start = float(lines[block + 3].rstrip())
                    self.energy[counter] = np.arange(data_points) * E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = (
                            float(lines[block + 19 + counter_inner])
                            / self.dwell[counter]
                        )
                self.note[counter] = ""
                print(self.energy.keys())
                print("Comments: {}".format(self.note))
                print("Dwell time: {}".format(self.dwell))
                print("Modes: {}".format(self.mode))
                print("Mode values: {}".format(self.mode_value))
            # ----------------------------------------------------------------------
            # New format
            if lines[6].lower().startswith("created with"):
                self.setup = "omicron"
                self.format = "New VAMAS"
                ENDING = "_1-Detector_Region.vms"

                # Do a search to find all files with matching name structure
                filename = pathlib.Path(filename)
                path = filename.parent
                filename = filename.name
                filen = filename.split("--")[0]
                search_for = filen + "*.vms"
                list_of_files = list(path.rglob(search_for))
                # Make sure the list is properly sorted
                try:
                    keys = [
                        int(str(name).split("--")[1].split("_")[0])
                        for name in list_of_files
                    ]
                except IndexError:
                    for i in list_of_files:
                        print(i)
                    raise
                keys.sort()
                list_of_files = [f"{filen}--{key}{ENDING}" for key in keys]
                self.scans = len(list_of_files)
                for counter, filename in enumerate(list_of_files):
                    # Load contents
                    with open(path / filename, "r") as f:
                        lines = f.readlines()
                        f.close()

                    # Analyze contents
                    blocks_4 = [
                        i
                        for i, line in enumerate(lines)
                        if (line.rstrip() == "-1")
                        and (lines[i + 1].lower().rstrip() == "kinetic energy")
                    ]
                    if len(blocks_4) > 1:
                        print(
                            "*** Interesting! More than 1 scan has been detected in above file!"
                        )
                    # Copy data points
                    i = blocks_4[0]
                    ###
                    if counter == 0:
                        _counter = 0
                        while True:
                            if lines[_counter].startswith("CREATION COMMENT START"):
                                comment_start = _counter
                                break
                            else:
                                _counter += 1
                                if _counter > len(lines):
                                    break
                        _counter = 0
                        while True:
                            if lines[_counter].startswith("CREATION COMMENT END"):
                                comment_end = _counter
                                break
                            else:
                                _counter += 1
                                if _counter > len(lines):
                                    break
                    self.note = lines[comment_start + 1 : comment_end]
                    ###
                    self.mode[counter] = lines[i - 11].rstrip()
                    self.mode_value[counter] = float(lines[i - 10].rstrip())
                    self.dwell[counter] = float(lines[i + 9].rstrip())
                    data_points = int(lines[i + 16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[i + 4].rstrip())
                    E_start = float(lines[i + 3].rstrip())
                    self.energy[counter] = np.arange(data_points) * E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = (
                            float(lines[i + 19 + counter_inner]) / self.dwell[counter]
                        )
        # ----------------------------------------------------------------------
        # Import Thetaprobe .avg data
        elif filename.endswith(".avg"):
            self.setup = "thetaprobe"
            with open(filename, "r", encoding="latin-1") as f:
                lines = f.readlines()

            # Check for ISS
            info = {
                line.split(":")[0].strip(): line.split("=")[1].strip()
                for line in lines
                if line.startswith("DS_")
            }
            if info["DS_ANPROPID_LENS_MODE_NAME"] != "'ISS'":
                print(
                    "{} does not appear to be an ISS experiment!".format(self.filename)
                )
                print(
                    "Expected 'ISS', but encountered: {}".format(
                        info["DS_ANPROPID_LENS_MODE_NAME"]
                    )
                )
                raise ImportError("File not an ISS experiment!")
            if info["DS_EXT_SUPROPID_CREATED"] == info["DS_EXT_SUPROPID_SAVED"]:
                # print('Created and saved dates are identical - checking for empty dataset...')
                check_empty = True
            else:
                check_empty = False

            # Metadata
            self.note[0] = info["DS_EXT_SUPROPID_SUBJECT"]
            self.date = info["DS_EXT_SUPROPID_CREATED"]
            self.dwell[0] = float(info["DS_ACPROPID_ACQ_TIME"])
            self.mode[0] = int(info["DS_ANPROPID_MODE"])
            self.mode_value[0] = float(info["DS_ANPROPID_PASS"])
            if info["DS_GEPROPID_VALUE_LABEL"] == "'Counts'":
                normalize = True  # normalize to "counts per second"
            else:
                normalize = False

            # Data
            # data_info = {}
            line_number = [
                i for i, line in enumerate(lines) if line.startswith("$DATAAXES")
            ]
            self.scans = 1
            if len(line_number) > 1:
                print("Reading file: {}".format(self.filename))
                raise ImportError("Import of multiple dataaxes not implemented yet!")
            else:
                line_number = line_number[0]
            keys = [
                key.strip() for key in lines[line_number - 1].split("=")[1].split(",")
            ]
            values = [
                key.strip() for key in lines[line_number + 1].split("=")[1].split(",")
            ]
            data_info = {key: value for key, value in list(zip(keys, values))}

            start, end = float(data_info["start"]), float(data_info["end"])

            # space_info = {}
            line_number = [
                i for i, line in enumerate(lines) if line.startswith("$SPACEAXES")
            ]
            if len(line_number) > 1:
                print("Reading file: {}".format(self.filename))
                raise ImportError("Import of multiple dataaxes not implemented yet!")
            else:
                line_number = line_number[0]
            keys = [
                key.strip() for key in lines[line_number - 1].split("=")[1].split(",")
            ]
            values = [
                key.strip() for key in lines[line_number + 1].split("=")[1].split(",")
            ]
            space_info = {key: value for key, value in list(zip(keys, values))}

            num = int(space_info["numPoints"])
            if space_info["linear"] != "LINEAR":
                print("Reading file: {}".format(self.filename))
                raise ImportError("Check .avg file if energy axis is linear!")

            # Generate xy-data
            self.energy[0] = np.linspace(start, end, num)
            self.cps[0] = self.energy[0] * np.nan

            line_number = [
                i for i, line in enumerate(lines) if line.startswith("$DATA=")
            ]
            if len(line_number) > 1:
                msg = "Reading file: {}".format(self.filename)
                raise ImportError("Import of multiple dataaxes not implemented yet!")
            else:
                line_number = line_number[0]

            for j in range(num):
                if j % 4 == 0:  # values are grouped in chunks of 4
                    line_number += 1
                    line = lines[line_number].split("=")[1].split(",")
                try:
                    self.cps[0][j] = float(line[j % 4])
                except ValueError:
                    pass  # #empty# values
            if check_empty:
                if not np.any(np.isfinite(self.cps[0])):
                    raise ImportError("Dataset from {} is empty!".format(self.filename))
                else:
                    print(
                        "Dataset appeared to be empty from the saved timestamps, but is not empty."
                    )
            if normalize:
                self.cps[0] /= self.dwell[0]
        else:
            raise IOError(
                'File: "{}" not found or fileending not accepted.'.format(self.filename)
            )

        # Print loaded settings
        print("Successfully loaded file: {}".format(filename))
        string = "Used settings:\nProbing mass: {} amu\nScatter angle: {}\nPrimary energy: {} eV"
        # print(string.format(*[self.settings[key] for key in ['mass', 'theta', 'E0']]))

    def __iter__(self):
        return DataIterator(self)

    @property
    def x(self):
        return self.energy[self.default_scan]

    @x.setter
    def x(self, var):
        if not var in self.energy.keys():
            print('"{}" not an available key! {}'.format(var, self.energy.keys()))
        self.default_scan = var

    @property
    def y(self):
        return self.cps[self.default_scan]

    @y.setter
    def y(self, var):
        if not var in self.energy.keys():
            print('"{}" not an available key! {}'.format(var, self.energy.keys()))
        self.default_scan = var

    @property
    def xy(self):
        return np.vstack((self.x, self.y)).T

    def get_xy(self, index):
        return np.vstack((self.energy[index], self.cps[index])).T

    @property
    def background(self):
        if self._background is not None:
            return self._background[self.default_scan]
        else:
            return None

    def convert_energy(self, mass):
        """Converts a measured energy to mass of surface atom
        corresponding the settings stored in the experiment.
        """
        angle = self.settings["theta"] * np.pi / 180
        return (
            self.settings["E0"]
            * (
                (
                    self.settings["mass"] * np.cos(angle)
                    + np.sqrt(
                        mass**2 - self.settings["mass"] ** 2 * np.sin(angle) ** 2
                    )
                )
                / (mass + self.settings["mass"])
            )
            ** 2
        )

    def plot_all_scans(self, exclude=[None], color=None):
        """Plot all elements in file in single figure."""
        selection = [i for i in self.energy.keys() if not i in exclude]
        if not color:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i])
        else:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i], color=color)
        plt.xlabel("Kinetic energy (eV)")
        plt.ylabel("Counts per second")

    def normalize(self, interval="Max", exclude=[None], unit="Mass", delta_e=10):
        """Normalize to highest value in interval=[value1, value2]"""
        self.delta_e = delta_e
        if isinstance(interval, int):
            self.normalization_criteria = interval
        elif isinstance(interval, str):
            if interval == "Total":
                self.normalization_criteria = "all"
            elif interval.lower().startswith("max"):
                self.normalization_criteria = "max"
            elif interval == "Au":
                self.normalization_criteria = 196.0
        if not isinstance(interval, list):
            if self.normalization_criteria == "all":
                selection = [i for i in range(self.scans) if not i in exclude]
                for __counter in selection:
                    total = simps(self.cps[__counter], self.energy[__counter])
                    self.cps[__counter] /= total
            elif self.normalization_criteria == "max":
                selection = [i for i in range(self.scans) if not i in exclude]
                for __counter in selection:
                    ydata = ct.smooth(self.cps[__counter], width=2)
                    norm_value = max(ydata)
                    self.cps[__counter] /= norm_value
            else:
                interval = [0, 0]
                if unit.lower() == "mass":
                    interval[0] = (
                        self.convert_energy(self.normalization_criteria) - self.delta_e
                    )
                    interval[1] = (
                        self.convert_energy(self.normalization_criteria) + self.delta_e
                    )
                elif unit.lower() == "energy":
                    interval[0] = self.normalization_criteria - self.delta_e
                    interval[1] = self.normalization_criteria + self.delta_e
                selection = [
                    i
                    for i in range(self.scans)
                    if (not i in exclude)
                    and (not interval[0] > max(self.energy[i]))
                    and (not interval[1] < min(self.energy[i]))
                ]
                for __counter in selection:
                    range_1 = np.where(self.energy[__counter] < interval[1])[0]
                    range_2 = np.where(self.energy[__counter] > interval[0])[0]
                    energy_range = np.intersect1d(range_1, range_2)
                    value = max(self.cps[__counter][energy_range])
                    self.cps[__counter] = self.cps[__counter] / value

    def add_mass_lines(
        self,
        masses,
        ax=None,
        offset=0,
        color="k",
        labels=True,
        linestyle="dotted",
        **kwargs,
    ):
        """Add vertical lines for mass references."""
        energies = self.convert_energy(np.array(masses))
        if ax is None:
            ax = plt.gca()
        [x1, x2, y1, y2] = ax.axis()
        for energy, mass in zip(energies, masses):
            ax.axvline(
                x=energy - offset,
                ymin=0,
                ymax=1,
                linestyle=linestyle,
                color=color,
                **kwargs,
            )
            if labels:
                ax.text(
                    float(energy) / x2,
                    0.95,
                    "m-{}".format(mass),
                    transform=ax.transAxes,
                )

    def add_regions(self):
        """Add regions indicating the whereabouts of 3d, 4d, 5d metals and the
        lanthanides and actinides."""
        ax = plt.gca()
        d3 = [45, 65]
        d4 = [89, 112]
        d5 = [178, 201]
        lant = [139, 175]
        act = [227, 260]
        for i in [d3, d4, d5]:
            ax.axvspan(
                xmin=self.convert_energy(i[0]),
                xmax=self.convert_energy(i[1]),
                color="k",
                alpha=0.2,
            )
        for i in [lant, act]:
            ax.axvspan(
                xmin=self.convert_energy(i[0]),
                xmax=self.convert_energy(i[1]),
                color="y",
                alpha=0.2,
            )
