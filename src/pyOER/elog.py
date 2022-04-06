"""This module defines the ElogEntry and a function for parsing elog html
See: https://elog.psi.ch/elog/

Made for DTU SurfCat's cinfelog by Soren B. Scott on July 29, 2020
"""
from pathlib import Path
import re
import json
from .constants import ELOG_DIR

SETUP = "ECMS"

FIELD_NAME_MATCHER = r"""<th class="listtitle"><a.*>(.*)</a></th>"""
ENTRY_NUMBER_MATCHER = (
    r"""<tr><td class="list1"><input type="checkbox" name="s[0-9]+" value="([0-9]+)">"""
)
FIELD_VALUES_START = r"""<td class="list1"""
FIELD_VALUE_MATCHER = r"""<td class="list1".*">([^<>]*)<"""
ENTRY_END_MATCHER = r"</pre></td></tr>"


def read_elog_html(
    path_to_elog_html,
    setup=SETUP,
):
    """Return a list of ElogEntry's with data from the html file

    TODO: replace this with something using BeautifulSoup
    """
    elog_entries = []  # this is what we will return
    field_names = []  # this will list the metadata field_names specified in the elog
    n_elog = None  # this will be the elog entry number
    notes_lines = []
    notestext = False  # a boolean to indicate whether or not we're reading notes
    gotvalues = False  # a boolean to indicate whether or not we got the field values

    with open(path_to_elog_html) as f:
        line = True  # need an initial value
        while line:
            # I prefer reading at start of while so it doesn't get skipped:
            try:
                line = f.readline()
            except UnicodeDecodeError:
                print(f"Error on line after line = {line}")
                continue
            elog_field_match = re.search(FIELD_NAME_MATCHER, line)
            if elog_field_match:
                # then this line names a field in the elog metadata
                field_names += [elog_field_match.group(1)]
                print(f"got field names = {field_names}")
                # don't need to waste time on the rest. Move on to the next line:
                continue

            entry_number_match = re.search(ENTRY_NUMBER_MATCHER, line)
            if entry_number_match:
                # then this line starts a new entry and specifies its number
                n_elog = int(entry_number_match.group(1))
                field_values = [n_elog]  # the first field value is ID
                print(f"working on entry number {n_elog}")

            if line.startswith(FIELD_VALUES_START):
                # ... then this specifies the values of the fields. great.
                # Dates get their own lines, but follow the same structure.
                # This will not match the ID, since that line is broken in the html.
                #     (html is disgusting)
                field_strings = line.split("</td>")
                if len(field_strings) > 1:  # will be the case if the line is real.
                    # The last entry is just noise and I want to toss it:
                    field_strings = field_strings[:-1]
                    # print(f"getting values from = f{field_strings}")  # debugging
                    # Now go through and use re to extract the value
                    for field_string in field_strings:
                        field_value_match = re.search(FIELD_VALUE_MATCHER, field_string)
                        if field_value_match:
                            field_values += list(field_value_match.groups())
                        else:
                            field_values += [None]

                    gotvalues = True  # the next non-blank non-field line is notes

            elif gotvalues and line.strip():  # then we're done with the field values
                # so we can package the field values
                field_data = dict(zip(field_names, field_values))
                # and get ready for notes
                notestext = True

            if notestext:  # then add the line to the notes of the entry!
                # notes = notes.join(line)  # <-- this leads to MemoryErrors :(
                notes_lines += [line]

            entry_end_match = re.search(ENTRY_END_MATCHER, line)
            if entry_end_match:
                # Then we just put it together:
                elog_entry = ElogEntry(
                    setup=setup,
                    number=n_elog,
                    field_data=field_data,
                    notes="".join(notes_lines),
                )
                # add the result to the list
                elog_entries += [elog_entry]
                # clear stored values and be fresh for next elog entry
                notes_lines = []
                field_values = []
                notestext = False
                gotvalues = False

    return elog_entries


def all_elog_entries(elog_dir=ELOG_DIR):
    """iterate over elog entries."""
    for file in elog_dir.iterdir():
        if file.suffix == ".json":
            with open(file, "r") as f:
                log_as_dict = json.load(f)
            try:
                yield ElogEntry(**log_as_dict)
            except TypeError as e:  # raised if __init__ doesn't get right args
                print(f"{e}... maybe {file} shouldn't be in {elog_dir}?")
                raise


class ElogEntry:
    """simple data class summarizing contents of an elog entry"""

    def __init__(
        self,
        setup,
        number,
        date=None,
        field_data=None,
        sample_measurements=None,
        measurement_EC_tags=None,
        notes=None,
    ):
        self.setup = setup
        self.number = number
        self.date = date
        self.field_data = field_data
        self.sample_measurements = sample_measurements
        self.measurement_EC_tags = measurement_EC_tags
        self.notes = notes

    @classmethod
    def load(cls, file_name, measurement_dir=ELOG_DIR):
        """Load the elog entry given its file name and dir"""
        path_to_file = Path(measurement_dir) / file_name
        with open(path_to_file, "r") as f:
            self_as_dict = json.load(f)
        return cls(**self_as_dict)

    @classmethod
    def open(cls, e_id, setup=SETUP, elog_dir=ELOG_DIR):
        """Open the elog entry given its id"""
        try:
            path_to_file = next(
                path
                for path in Path(elog_dir).iterdir()
                if path.stem.startswith(f"{setup} {e_id}")
            )
        except StopIteration:
            raise FileNotFoundError(f"no elog with number={e_id}")
        return cls.load(path_to_file)

    def save(self, file_name=None, elog_dir=ELOG_DIR):
        """Save the elog entry, uses self.get_name for file name by default"""
        self_as_dict = self.as_dict()
        if not file_name:
            file_name = self.get_name()
        path_to_json = (Path(elog_dir) / file_name).with_suffix(".json")
        with open(path_to_json, "w") as f:
            json.dump(self_as_dict, f, indent=4)

    def as_dict(self):
        """Return a dictionary representation of the elog entry"""
        self_as_dict = dict(
            setup=self.setup,
            number=self.number,
            date=self.date,
            field_data=self.field_data,
            sample_measurements=self.sample_measurements,
            measurement_EC_tags=self.measurement_EC_tags,
            notes=self.notes,
        )
        return self_as_dict

    def get_name(self):
        return f"{self.setup} {self.number} {self.date}"

    def __repr__(self):
        return f"{self.setup} {self.number} {self.date}"

    def update_with(self, elog_2):
        """Put non-empty attributes from elog_2 into self, update not replace dicts"""
        new_attrs = elog_2.as_dict()
        for attr, value in new_attrs.items():
            if (
                hasattr(self, attr)
                and isinstance(getattr(self, attr), dict)
                and isinstance(value, dict)
            ):
                getattr(self, attr).update(value)
            elif value:
                print(f"update is setting {attr} to {value}")
                setattr(self, attr, value)
