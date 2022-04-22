from pyOER import Experiment
from matplotlib import pyplot as plt

forpublication = True
if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    # plt.rc('text', usetex=True)  # crashingly slow
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.6)
else:
    plt.style.use("default")

e1 = Experiment.open(33)

e1.measurement.meas.calibrate(
    RE_vs_RHE=e1.measurement.RE_vs_RHE,
    A_el=0.196,
)

axes_b = e1.measurement.plot(
    tspan=None,
    mass_list=[
        # "M4", "M18",
        "M28",
        "M32",
        "M34",
        "M36",
        # "M40",  # yup, moves with M36.
        "M44",
        "M46",
        "M48",
    ],
    unit="A",
)

axes_b[0].set_ylim([1e-13, 1e-7])
axes_b[0].set_ylabel("raw MS signal / (A)")
axes_b[1].set_ylabel("E vs RHE / (V)")
axes_b[2].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
axes_b[0].set_xlabel("time / (s)")
axes_b[1].set_xlabel("time / (s)")

fig_b = axes_b[0].get_figure()
fig_b.set_figheight(fig_b.get_figwidth())
if forpublication:
    fig_b.savefig("paper_II_v4_fig1b.png")
    fig_b.savefig("paper_II_v4_fig1b.svg")
