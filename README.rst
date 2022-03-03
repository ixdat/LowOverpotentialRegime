LowOverpotentialRegime
======================

Data and analysis for _________________
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This repository has gotten the data analysis etc for the DTU-MIT RuO2 project organized and under version control.
It will also be used to open-source the raw data when we publish.

The core code is organized int a python package ``pyOER``, the modules of which are in src/pyOER.
This package must be added to your python path for other parts to run.

The package uses `ixdat <https://ixdat.readthedocs.io/>`_ to import and handle electrochemistry - mass spectrometry 
(EC-MS) data, but defines its own internal relational data structure for higher-level organization.

This repository may eventually be redone using ``ixdat`` for higher data organization if there is future interest.

Setup
-----

Install ``ixdat`` v 0.1.6 or above using::

  pip install --upgrade ixdat

At present, a settings.py file is required to indicate the location of raw data (which is not yet publically available)
Create **src/pyOER/settings.py** and copy the following to it::

  from pathlib import Path

  # Change below path to match the path to the shared folder of your project.
  # DATA_DIR = Path("/absolute/path/to/your/shared/project-folder/")
  DATA_DIR = Path(r"C:\Users\scott\Dropbox\WORKSPACES\DTU-MIT RuO2")

Change the DATA_DIR to the location of the project with the data files as it is called on your computer.

Figures
-------
This folder includes the scripts needed to make the figures in the manuscripts, which are in preparation.
Each script uses data in the database (see below) and the pyOER package, and produces as .png and .svg
the figures as they are in the publication (minus some annotations which are added later in Inkscape)

These scripts are full of comments in order to serve as tutorials on the use of pyOER and the database.

The figures themselves will be added when the works are published, and this README will be updated
with a brief description of each figure script, like here: https://github.com/ScottSoren/pyCOox_public.

Database
--------

The relational database is represented in tables/ . Each subfolder of table/ acts as a database table, with
each .json file acts as a row. The .json files can be opened with any text editor (e.g. Notepad) for ease of
use. Each table corresponds to a class in pyOER, and each row corresponds to an object of that class.

The tables include but are not limited to:

- **measurements**. A measurement is a collection of metadata about a single EC, MS or EC-MS measurement.
It includes a pointer to any notes taken during the measurement (elog table), as well as the raw data file.
pyOER is lazy in that raw data is not loaded by the measurement object until it is needed, saving RAM.

- **icpms**. An row in this table represents a single ICPMS measurement, with the counts and metadata
including which element was measured, which ICPMS calibration it corresponds to, which measurement the
electrolyte sample was taken during (measurement table) and when in that measurement.

- **experiments**. An experiment contains the additional metadata needed to do a standard analysis of
a measurement (measurements table) to extract results such as activity, faradaic efficiency, lattice oxygen evolution rate,
and dissolution rates. The additional metadata can include references to where in the measurement or elsewhere
to read the calibration of the O2 signal at m/z=32 and the 16O/18O isotopic ratio in the electrolyte

- **tofs**. A row in this table represents a result, which is a rate of oxygen evolution, lattice oxygen evolution,
or dissolution. It also includes a capacitance value for normalization to number of sites (thus the name of the
table, for turn-over frequency). Each row also includes all the metadata used to derive the result - most importantly a pointer
to the corresponding experiment (experiment table) and the timespan during that experiment for which the rate applies.

- **elog**. A row in this table includes an entry in the electronic lab notebook. They contain metadata and the
some of the experimenter's thoughts during the measurement, and can be useful if a specific measurement
seems hard to interpret.

The tables were populated semi-automatically using the scripts in data_processing/

Raw data is at present available only to the authors via a dropbox folder.
It will be made publically available upon publication.

