"""This moves the LEIS raw data to the dropbox folder from which it will be shared
It also converts it from pickle-5 to standard (as of python 3.7) pickle 4.0.
"""

import shutil
from pathlib import Path
from pyOER.iss import ISS
import pickle5
import pickle

new_data_dir_str = "~/Dropbox/DATA/LowOverpotentialRegime/LEIS"
new_data_dir = Path(new_data_dir_str).expanduser()

old_data_dir = Path("~/Dropbox/WORKSPACES/DTU-MIT RuO2/Data/ISS").expanduser()

main_iss = ISS()

for sample_name in main_iss.all_leis():
    main_iss.get_sample(sample_name)
    keys = [key for key in main_iss._meta['data']]
    file_names = [main_iss._meta['data'][key]['pickle_name'] for key in keys]
    # print(sample_name)
    # print(filenames)
    for i, file_name in enumerate(file_names):
        old_path_to_file = old_data_dir / "organized_pickles" / file_name
        new_path_to_file = new_data_dir / "spectra" / (file_name.split(".")[0] + ".pkl")
        with open(old_path_to_file, "rb") as old:
            data = pickle5.load(old)
        with open(new_path_to_file, "wb") as new:
            pickle.dump(data, new)
        main_iss._meta["data"][i]["pickle_name"] = new_path_to_file.name

    old_path_to_extra_file = old_data_dir / "pickled_pickles" / (sample_name + ".pickle")
    new_path_to_extra_file = new_data_dir / "extras" / (sample_name + ".pkl")
    with open(old_path_to_extra_file, "rb") as old:
        extra_data = pickle5.load(old)
    with open(new_path_to_extra_file, "wb") as new:
        pickle.dump(extra_data, new)

    main_iss.save_json()

