"""This moves the EC-MS raw data to the dropbox folder from which it can be shared"""

import shutil
import re
from pathlib import Path
from pyOER import all_measurements


new_data_dir_str = "~/Dropbox/DATA/LowOverpotentialRegime/EC_MS"
new_data_dir = Path(new_data_dir_str).expanduser()

continue_from = 0
for m in all_measurements():
    if m.id < continue_from:
        continue

    if True:  # just make sure it worked
        print(m)
        if "broken" not in m.name:
            try:
                print("\t" + m.meas.name)
            except Exception as e:
                print("\t" + str(e))
        continue
    if False:  # do the move!
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
    if False:  # Simplify the path by removing levels between 'EC_MS' and dated folder
        # The pickle files were moved manually, so here is just updating in the jsons's
        old_data_path = m.raw_data_path
        later_part = str(old_data_path).split("EC_MS")[-1]
        later_part = later_part.replace("\\", "/").strip("/")
        parts = later_part.split("/")

        skipping = True
        new_data_path = new_data_dir_str
        for part in parts:
            if re.search("^[0-9]{2}[A-L][0-9]{2}", part):
                skipping = False
            if not skipping:
                new_data_path = new_data_path + "/" + part
        m.raw_data_path = new_data_path
        m.save()


