from pyOER import StandardExperiment
import numpy as np
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

lines = []

for panel, e_id in {"a": 33, "b": 41}.items():  # 5, 39

    e = StandardExperiment.open(e_id)

    axes = e.plot_EC_MS_ICPMS()
    ylim_left = np.array([-1, 10])
    axes[0].set_ylim(ylim_left)
    axes[3].set_ylim(ylim_left / e.beta)

    if e_id in [33, 5]:
        axes[-1].set_ylim([0, 25])
    if e_id in [39, 41]:
        axes[-1].set_ylim([0, 2.5])

    axes[0].set_ylabel(
        "$^{16}$O$^{18}$O and $^{18}$O$_2$ / \n" "(pmol s$^{-1}$cm$^{-2}_{geo})$"
    )
    axes[1].set_ylabel("E vs RHE / (V)")
    axes[2].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
    axes[3].set_ylabel("$^{16}$O$_2$ / \n" "(pmol s$^{-1}$cm$^{-2}_{geo})$")
    axes[4].set_ylabel(
        f"{e.sample.element} dissolution / \n" "(pmol s$^{-1}$cm$^{-2}_{geo})$"
    )
    fig = axes[0].get_figure()
    if forpublication:
        fig.set_figwidth(3.25)
        fig.set_figheight(3.25)
        fig.savefig(f"paper_II_v4_fig2{panel}.png")
        fig.savefig(f"paper_II_v4_fig2{panel}.svg")

    lines += ["\n", f"----- {panel} = {e}:\n"]
    tofs = e.get_tofs()

    for tof in tofs:
        lines += [
            f"{tof.tof_type} rate at {tof.description} = {tof.rate} [mol/s] \n",
            (
                f"\t= {tof.amount * 1e12 / 0.196} [pmol/cm^2] "
                + f"in the {tof.t_interval} [s] interval\n"
            ),
            f"\twith average potential = {tof.potential} [V] vs RHE.\n",
        ]


with open("paper_II_v4_fig2_annotation.txt", "w") as f:
    f.writelines(lines)
