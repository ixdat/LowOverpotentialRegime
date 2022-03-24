LowOverpotentialRegime
======================

Data and analysis for:

Soren B. Scott and Reshma R. Rao, Choongman Moon, Jakob Ejler Sørensen, Jakob Kibsgaard, Yang Shao-Horn, and Ib Chorkendorff. **The low overpotential regime of acidic water oxidation part I: The importance of O2 detection**. `Energy & Environmental Science. In Press. <https://doi.org/10.1039/D1EE03914H>`_

and


Soren B. Scott, Jakob Ejler Sørensen,  Reshma R. Rao, Choongman Moon, Jakob Kibsgaard, Yang Shao-Horn, and Ib Chorkendorff. **The low overpotential regime of acidic water oxidation part II: Trends in metal and oxygen stability numbers**.  `Energy & Environmental Science. In Press. <https://doi.org/10.1039/D1EE03915F>`_

Setup
-----

We have just started moving things over to this repository from an old repository, where we developed the scripts before the articles were submitted. 

Please have patience!

For now you can can see the old repository at https://github.com/ScottSoren/pyOER20
See its README for instructions. 


Figures
-------
Links for now point towards the old repository. We will update these when we move everything here. 

`Part I: The importance of O2 detection. <https://doi.org/10.1039/D1EE03914H>`_
...............................................................................

Figure 1
^^^^^^^^
(a), (b), and (c) are diagrams. 

(d) and (e):

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/paper_I_v6_fig1.py

Requires downloading the raw data files.

These panels compare working electrode current and oxygen production during cyclic voltammatry on RuOx. 

Note also that we added annotations to most figures in Inkscape after making the panels with these scripts.

Figure 2
^^^^^^^^

(a), (b), (c), and (d):

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/paper_I_v6_fig2.py

(a), (b), and (d) show plots of activity measurements. (c) shows an analysis of the Faradaic Efficiency during the activity measurement in (b).

Requires downloading the raw data files.

Figure 3
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/paper_I_v6_fig3.py

This script defines a function ``plot_all_activity_results`` which is imported by fig5.py amoung others. 
This function is used to make a figure in the lower part (under ``if __name__ == "__main__"``

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 4
^^^^^^^^

This was made in Origin.

Figure 5
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/paper_I_v6_fig5.py

In addition to making the figure (using an import from fig3.py, this script also solves the model described in `the ESI <https://www.rsc.org/suppdata/d1/ee/d1ee03914h/d1ee03914h1.pdf>`_ for j0, dG1 and dG2.

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure S1
^^^^^^^^^

On its way.

Figure S2
^^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/SI_paper_I_v6_S2_O2.py

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

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_I_v6_figs/SI_paper_I_v6_fig_S8_equilibrium.py

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

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_II_v4_figs/paper_II_v4_fig2.py

Requires downloading the raw data files.

Figure 3
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_II_v4_figs/paper_II_v4_fig3.py

Requires downloading the raw data files.

Figure 4
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_II_v4_figs/paper_II_v4_fig4.py

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 5
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/figures/paper_II_v4_figs/paper_II_v4_fig5.py

This uses the plotting function from Figures 3 and 5 of Paper 1.

Does not require downloading the raw data, as it uses the results stored in the tables of the repository.

Figure 6
^^^^^^^^

https://github.com/ScottSoren/pyOER20/blob/master/src/plot_ruo2_leis_figures.py

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

On its way.
