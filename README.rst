LowOverpotentialRegime
======================

Data and analysis for:

Soren B. Scott and Reshma R. Rao, Choongman Moon, Jakob Ejler Sørensen, Jakob Kibsgaard, Yang Shao-Horn, and Ib Chorkendorff. **The low overpotential regime of acidic water oxidation part I: The importance of O2 detection**. `Energy & Environmental Science. In Press. <https://doi.org/10.1039/D1EE03914H>`_

and


Soren B. Scott, Jakob Ejler Sørensen,  Reshma R. Rao, Choongman Moon, Jakob Kibsgaard, Yang Shao-Horn, and Ib Chorkendorff. **The low overpotential regime of acidic water oxidation part II: Trends in metal and oxygen stability numbers**.  `Energy & Environmental Science. In Press. <https://doi.org/10.1039/D1EE03915F>`_

Setup
-----

1. Clone or download the repository

2. Install ``ixdat`` version 0.2.2 or above using::

     pip install --upgrade ixdat

3. To enable full analysis, download the raw data using this dropbox link:
   https://www.dropbox.com/sh/u0ttselmhd08ntq/AAA43jVl42MB--sV5pyf_6CPa?dl=0

   The size is 3 GB.

   a. By default, the code will look for raw data in ``~/Dropbox/DATA/LowOverpotentialRegime``.
      If you download the data to another location, you must create the file ``src/pyOER/settings.py``
      and in it, define the ``DATA_DIR`` variable. An example ``settings.py``::

          from pathlib import Path

          # Change below path to match the path to the shared folder of your project.
          DATA_DIR = Path(r"C:\DATA\other_peoples_data\LowOverpotentialRegime")


Figures
-------

`Part I: The importance of O2 detection. <https://doi.org/10.1039/D1EE03914H>`_
...............................................................................

Figure 1
^^^^^^^^
(a), (b), and (c) are diagrams. 

(d) and (e):

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_main_figs/part_I_fig1.py

Requires downloading the raw data files.

These panels compare working electrode current and oxygen production during cyclic voltammatry on RuOx. 

Note also that we added annotations to most figures in Inkscape after making the panels with these scripts.

Figure 2
^^^^^^^^

(a), (b), (c), and (d):

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_main_figs/part_I_fig2.py

(a), (b), and (d) show plots of activity measurements. (c) shows an analysis of the Faradaic Efficiency during the activity measurement in (b).

Requires downloading the raw data files.

Figure 3
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_main_figs/part_I_fig3.py

This script defines a function ``plot_all_activity_results`` which is imported by fig5.py amoung others. 
This function is used to make a figure in the lower part (under ``if __name__ == "__main__"``

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 4
^^^^^^^^

This was made in Origin. Contact R. R. R.

Figure 5
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_main_figs/part_I_fig5.py

In addition to making the figure (using an import from fig3.py, this script also solves the model described in `the ESI <https://www.rsc.org/suppdata/d1/ee/d1ee03914h/d1ee03914h1.pdf>`_ for j0, dG1 and dG2.

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure S1
^^^^^^^^^

On its way.

Figure S2
^^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_SI_figs/part_I_figS2.py

Requires downloading the raw data files.

Figure S3
^^^^^^^^^

On its way.

Figure S4
^^^^^^^^^

Made in Origin

Figure S5
^^^^^^^^^

Made in Origin

Figure S6
^^^^^^^^^

Made in Origin

Figure S7
^^^^^^^^^

Made in Origin

Figure S8
^^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_I_SI_figs/part_I_figS8.py

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

`Part II: Trends in metal and oxygen stability numbers. <https://doi.org/10.1039/D1EE03915F>`_
.................................................................................................

Figure 1
^^^^^^^^

This is a diagram made with Inkscape.

Figure 2
^^^^^^^^

(a) is a diagram.

(b):

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_main_figs/part_II_fig2.py

Requires downloading the raw data files.

Figure 3
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_main_figs/part_II_fig3.py

Requires downloading the raw data files.

Figure 4
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_main_figs/part_II_fig4.py

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 5
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_main_figs/part_II_fig5.py

This uses the plotting function from Figures 3 and 5 of Paper 1.

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 6
^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_main_figs/part_II_fig6.py

Figure 7
^^^^^^^^

This is a diagram.

Figure S1
^^^^^^^^^

This is a diagram.

Figure S2
^^^^^^^^^

On its way.

Figure S3
^^^^^^^^^

On its way.

Figure S4
^^^^^^^^^

On its way.

Figure S5
^^^^^^^^^

https://github.com/ixdat/LowOverpotentialRegime/blob/move_leis/figures/part_II_SI_figs/part_II_figS5.py
