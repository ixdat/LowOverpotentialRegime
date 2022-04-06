"""This module groups TurnOverFrequency's for plotting and analysis"""

import json
import numpy as np
from .experiment import all_standard_experiments
from .tof import TurnOverFrequency


def nested_update(dict_1, dict_2):
    """Update dict_1 with dict_2, but update rather than replace values which are dicts.

    This is a nice tool to use later as we turn the datapoints into statistics
    """
    for key, value in dict_2.items():
        if key in dict_1 and isinstance(dict_2[key], dict):
            nested_update(dict_1[key], dict_2[key])
        else:
            dict_1[key] = value


def nested_update_with_layers(dict_1, dict_2, layers, **kwargs):
    """Update dict_1 with dict_2, filtering away some layers.

    Example:
        >>> nested_update_with_layers(
        ...     dict_1 = {},
        ...     dict_2 = {"a": {"ii": {"1": 121, "2": 122},
        ...                     "iii": {"3": 133, "4": 134}}},
        ...     layers = ("abc", "roman", "arabic"),
        ...     roman="ii"
        ... )
        OUT: {"a": {"1": 121, "2": 122}}

    Args:
        dict_1 (dict): The dictionary to update. Can be an empty dictionary.
        dict_2 (dict): The dictionary to update it with. This will in general be a
            homogeneous multi-layered string-keyed dictionary.
        layers (iter of str): The names of the layers of dict_2. Needs to have the
            same length as the number of times you would index dict_2.
        kwargs: The filters to apply to dict_2 while updating dict_1. Each key-word
            has to be the name of a layer. Those layers will not be recreated from
            dict_1 to dict_2. Instead, the key-word argument will be used as an index,
            selecting only one branch of the lower layers of dict_2 to update into
            dict_1.
    """
    for layer_name in kwargs:
        if layer_name not in layers:
            key = kwargs.pop(layer_name)
            print(
                f"Can't filter based in {layer_name}={key} "
                f"because {key} is not in layers={layers}. This will be ignored."
            )
    if layers[0] in kwargs:
        new_kwargs = kwargs.copy()
        key = new_kwargs.pop(layers[0])
        nested_update_with_layers(dict_1, dict_2[key], layers[1:], **new_kwargs)
    else:
        for key, value in dict_2:
            if isinstance(dict_2, dict):
                if key not in dict_1:
                    dict_1[key] = {}
                nested_update_with_layers(
                    dict_1[key], dict_2[key], layers[1:], **kwargs
                )
            else:
                dict_1[key] = value
    return dict_1


def get_sample_type(sample, sample_mapping):
    """Return the key in sample_mapping that has a match to the start of sample

    If two sample types include a sample name base in their sample name list that match
    the start of sample, then return the one matching sample more specificically, i.e.
    that with the longest match.
    """
    possibilities = {
        s: s_type
        for (s_type, s_list) in sample_mapping.items()
        for s in s_list
        if sample.startswith(s)
    }
    return possibilities[max(possibilities.keys(), key=len)]  # woa, funky code :D


def get_current_point(tof, current_point_mapping):
    """Return a string which is the nearest current in mA/cm2"""
    j = tof.current * 1e3 / tof.measurement.A_el  # current in uA/cm^2
    errs = list(zip(*[(s, np.abs(f - j)) for s, f in current_point_mapping.items()]))
    index = int(np.argmin(errs[1]))
    current_point = errs[0][index]
    return current_point


class StabilityResultsCollection:
    def __init__(
        self,
        tof_collection=None,
        layers=("sample_type", "current_point", "result_time", "result_type"),
        sample_mapping=None,
        current_point_mapping=None,
    ):
        # okay, hold your horses, four-layer nested dictionary here.
        # Layers are (in the order that you'd index the dictionary):
        #   1. sample_type (RT-RuO2, ..., IrOx/Ir, ...).
        #   2. current_point (0.5 mA/cm^2, ... 0.05 mA/cm^2)
        #   3. Timespan: start of electrolysis or steady-state
        #   4. rate_type (activity, dissolution, or exchange)
        # Having all of this connected in a relational way with the raw data is why
        # pyOER is such a big package. But it will be done better with ixdat.

        # tof_collection will have a list of tof_id's (remember tof is the
        # unfortunate name given to the main results table in pyOER)
        self.tof_collection = tof_collection
        self.layers = layers
        if not tof_collection:
            self.generate_tof_collection(sample_mapping, current_point_mapping)

    def generate_tof_collection(self, sample_mapping, current_point_mapping):
        tof_collection = {
            sample_type: {
                current_point: {
                    tof_time: {
                        rate_type: []
                        for rate_type in ["activity", "dissolution", "exchange"]
                    }
                    for tof_time in ["start", "steady"]
                }
                for current_point in current_point_mapping.keys()
            }
            for sample_type in sample_mapping.keys()
        }

        for e in all_standard_experiments():
            try:
                sample_type = get_sample_type(e.sample_name, sample_mapping)
            except ValueError:
                print(f"couldn't find a sample type match for {e}. Skipping.")
                continue
            tofs = e.get_tofs()
            for tof in tofs:
                rate_type = tof.tof_type
                current_point = get_current_point(tof, current_point_mapping)
                if "steady" in tof.description or "composite" in tof.description:
                    rate_time = "steady"
                elif "start" in tof.description or "first" in tof.description:
                    rate_time = "start"
                else:
                    print(
                        f"{tof} with description = '{tof.description}' as "
                        f"it seems to be neither start nor steady"
                    )
                    continue
                tof_collection[sample_type][current_point][rate_time][rate_type].append(
                    tof.id
                )

        self.tof_collection = tof_collection

    def __getitem__(self, key):
        return StabilityResultsCollection(
            tof_collection=self.tof_collection[key], layers=self.layers[1:]
        )

    def __iter__(self):
        yield from self.tof_collection

    def keys(self):
        yield from self

    def items(self):
        for key in self:
            yield key, self[key]

    def get_sub_collection(self, **kwargs):
        new_tof_collection = {}
        old_tof_collection = self.tof_collection.copy()

        new_layers = []
        for i, layer in enumerate(self.layers):
            if layer not in kwargs:
                new_layers.append(layer)

        nested_update_with_layers(
            new_tof_collection, old_tof_collection, self.layers, **kwargs
        )
        return StabilityResultsCollection(
            tof_collection=new_tof_collection, layers=tuple(new_layers)
        )

    def as_dict(self):
        return {"tof_collection": self.tof_collection, "layers": self.layers}

    def save(self, path_to_file):
        with open(path_to_file, "w") as f:
            json.dump(self.as_dict(), f, indent=4)

    @classmethod
    def open(cls, path_to_file):
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        return cls(**self_as_dict)

    def get_coherent_results(self, sample_type, current_point, tof_time):
        act = []
        diss = []
        exc = []
        for t_id in self.tof_collection[sample_type][current_point][tof_time][
            "activity"
        ]:
            act_tof, diss_tof, exc_tof = TurnOverFrequency.open(t_id).get_tof_triplet()
            if not (act_tof and diss_tof and exc_tof):
                continue
            act.append(act_tof.rate)
            diss.append(diss_tof.rate)
            exc.append(exc_tof.rate)
        return {
            "activity": np.array(act),
            "dissolution": np.array(diss),
            "exchange": np.array(exc),
        }

    def get_results(self, sample_type, current_point, tof_time):

        diss = np.array(
            [
                TurnOverFrequency.open(t_id).rate
                for t_id in self.tof_collection[sample_type][current_point][tof_time][
                    "dissolution"
                ]
            ]
        )
        exc = np.array(
            [
                TurnOverFrequency.open(t_id).rate
                for t_id in self.tof_collection[sample_type][current_point][tof_time][
                    "exchange"
                ]
            ]
        )
        act = np.array(
            [
                TurnOverFrequency.open(t_id).rate
                for t_id in self.tof_collection[sample_type][current_point][tof_time][
                    "activity"
                ]
            ]
        )
        return {"activity": act, "dissolution": diss, "exchange": exc}

    def get_stats(self, sample_type, current_point, tof_time):
        """Return {type: [mean, standard deviation]} for types of stability numbers.

        The types of stability numbers are "S_number" for OER / metal dissolution and
        "S_number_lattice for OER / lattice oxygen evolution

        Args:
            sample_type (str): The type of sample, e.g. "RT-RuO2", "IrOx/Ir", etc.
            current_point (str): The current in [mA/cm^2], e.g. "0.5", "0.15", etc.
            tof_time (str): The part of the measurement. "start" or "steady".
        """
        stats = {sample_type: {current_point: {"S_number": {}, "S_number_lattice": {}}}}

        results = self.get_results(sample_type, current_point, tof_time)
        diss = results["dissolution"]
        exc = results["exchange"]
        act = results["activity"]
        # get rid of all nan's and inf's!
        diss = diss[~np.isnan(diss)]
        diss = diss[~np.isinf(diss)]
        exc = exc[~np.isnan(exc)]
        exc = exc[~np.isinf(exc)]
        act = act[~np.isnan(act)]
        act = act[~np.isinf(act)]

        S_number = np.mean(act) / np.mean(diss)
        S_number_lattice = np.mean(act) / np.mean(exc)
        if len(diss) > 1:
            sigma_S = S_number * np.sqrt(
                (np.std(act) / np.mean(act)) ** 2 + (np.std(diss) / np.mean(diss)) ** 2
            )  # propagate uncertainty in activity and uncertainty in dissolution.
            stats["S_number"] = [S_number, sigma_S]
            # [average, uncertainty]
        else:
            stats["S_number"] = [S_number, None]
            # no uncertainty since there's only one.
        if len(exc) > 1:
            sigma_S_lattice = S_number_lattice * np.sqrt(
                (np.std(act) / np.mean(act)) ** 2 + (np.std(exc) / np.mean(exc)) ** 2
            )  # propagate uncertainty in activity and uncertainty in exchange.
            stats["S_number_lattice"] = [
                S_number_lattice,  # average
                sigma_S_lattice,  # uncertainty (standard deviation)
            ]
        else:
            stats["S_number_lattice"] = [
                S_number_lattice,  # average (since there's only one
                None,  # can't quantify uncertainty if there's only one
            ]

        return stats
