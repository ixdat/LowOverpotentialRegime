"""This moves the EC-MS raw data to the dropbox folder from which it can be shared"""

import shutil
from pathlib import Path
from pyOER import all_measurements


new_data_dir_str = "~/Dropbox/DATA/LowOverpotentialRegime/EC_MS"
new_data_dir = Path(new_data_dir_str).expanduser()

continue_from = 217
for m in all_measurements():
    if True:  # make sure it worked
        print(m)
        if "broken" not in m.name:
            try:
                print("\t" + m.meas.name)
            except Exception as e:
                print("\t" + str(e))
        continue
    if m.id < continue_from:
        continue
    old_data_path = m.raw_data_path
    later_part = str(old_data_path).split("DTU-MIT RuO2")[-1]
    later_part = later_part.replace("\\", "/").strip("/")
    new_data_path = new_data_dir / later_part
    new_data_path_str = new_data_dir_str + "/" + later_part
    new_data_folder = new_data_path.parent
    if not new_data_folder.exists():
        new_data_folder.mkdir(parents=True)
    shutil.copyfile(old_data_path, new_data_path)
    m.raw_data_path = new_data_path_str
    m.save()
