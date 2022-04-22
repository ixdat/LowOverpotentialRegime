"""Generates the RuO2 LEIS figures for my phd thesis."""
import numpy as np
import matplotlib.pyplot as plt

import pyOER
from pyOER.iss import get_common_y
from pyOER.tools import get_range

forpublication = True
tex = False
if forpublication:  # for the publication figure
    import matplotlib as mpl

    mpl.rcParams["figure.figsize"] = (3.25, 2.75)
    if tex:
        plt.rc("text", usetex=tex)  # crashingly slow
    plt.rc("font", family="sans-serif")
    plt.rc("font", size=8)
    plt.rc("lines", linewidth=0.6)
    plt.rc("lines", markersize=3)  # extra
else:
    plt.style.use("default")

save = False
save_path = "../figures/extras/leis_tmp/"

colors = {
    "labeled": (22, 135, 52),
    "OER": (0, 147, 176),
    "sputter": (176, 120, 0),
}
description = [
    ("As-prepared", 640),
    ("After OER", 640),
    (r"OER \& Sputtered", 620) if tex else ("OER & Sputtered", 570),
]
for key, color in colors.items():
    r, g, b = color
    color = (r / 255, g / 255, b / 255)
    colors[key] = color
colors = list(colors.values())
print(colors)

# Main ISS handle
handle = pyOER.ISS()
samples = ["Taiwan1C", "Decade1G", "John4A", "John4C"]
selection = {
    "Taiwan1C": [0, 2, 3],
    "Decade1G": [2, 3, 4, 6],
    "John4A": [0, 2, 3],
    "John4C": [0, 2],
}
offsets = {
    "Taiwan1C": 15,
    "Decade1G": 10,
    "John4A": 20,
    "John4C": 20,
}
ylim = {
    "Taiwan1C": 50,
    "Decade1G": 40,
    "John4A": 65,
    "John4C": 45,
}
y0 = 1000
################################
########## Figures #############
# sample = 'Taiwan1C'
for sample in samples:
    print(f"\n{sample}")
    handle.get_sample(sample)
    handle.region = "oxygen"
    # ratios, coeffs = handle.fit_with_reference(peaks=[[16, 18]], plot_result=False)
    # TODO: add the background subtraction to data, so it is remembered

    fig = plt.figure(sample)
    ax = fig.add_subplot(111)
    handle.data.add_mass_lines(
        [16, 18, 101], labels=False, color="lightgray", linewidth=1.2
    )

    offset = offsets[sample]
    i = 0
    for index, data in enumerate(handle):
        if index not in selection[sample]:
            continue
        if sample == "Decade1G":
            if index in [2, 3] and data.comment == "spot":
                # wrong dataset
                continue
        j = handle.active
        print(f'{i}\t- O16: {handle.meta("O16")}')
        print(f' \t- O18: {handle.meta("O18")}\n')

        # Plot components
        setup = data.setup
        ref1 = handle._ref[setup][16]
        ref2 = handle._ref[setup][18]
        x = ref1["xy"][:, 0]
        mask = get_range(x, 300, 600)

        # Sum of components
        mask = get_range(handle.data.x, 300, 600)
        plt.plot(
            handle.data.x[mask],
            (
                i * offset
                + get_common_y(
                    handle.data.x,
                    ref1["xy"][:, 0],
                    ref1["peak"],
                )[mask]
                * handle.meta("c_16")
                / y0
                + get_common_y(
                    handle.data.x,
                    ref2["xy"][:, 0],
                    ref2["peak"],
                )[mask]
                * handle.meta("c_18")
                / y0
                + handle.background("y")[mask] / y0
            ),
            "k",
            label="Sum of components",
            # linewidth=0.5,
        )

        # Individual components:
        plt.plot(
            handle.data.x[mask],
            (
                i * offset
                + handle.background("y")[mask] / y0
                + get_common_y(
                    handle.data.x,
                    ref1["xy"][:, 0],
                    ref1["peak"],
                )[mask]
                * handle.meta("c_16")
                / y0
            ),
            "r-",
            label="O-16 component",
        )
        plt.fill_between(
            handle.data.x[mask],
            (
                i * offset
                + handle.background("y")[mask] / y0
                + get_common_y(
                    handle.data.x,
                    ref1["xy"][:, 0],
                    ref1["peak"],
                )[mask]
                * handle.meta("c_16")
                / y0
            ),
            handle.background("y")[mask] / y0 + i * offset,
            color=(250 / 255, 210 / 255, 210 / 255),
        )
        plt.plot(
            handle.data.x[mask],
            (
                i * offset
                + handle.background("y")[mask] / y0
                + get_common_y(
                    handle.data.x,
                    ref2["xy"][:, 0],
                    ref2["peak"],
                )[mask]
                * handle.meta("c_18")
                / y0
            ),
            "g-",
            label="O-18 component",
        )
        plt.fill_between(
            handle.data.x[mask],
            (
                i * offset
                + handle.background("y")[mask] / y0
                + get_common_y(
                    handle.data.x,
                    ref2["xy"][:, 0],
                    ref2["peak"],
                )[mask]
                * handle.meta("c_18")
                / y0
            ),
            handle.background("y")[mask] / y0 + i * offset,
            color=(217 / 255, 255 / 255, 232 / 255),
        )

        # Plots data and background
        ax.plot(
            handle.background("x"),
            handle.background("y") / y0 + i * offset,
            ":",
            color=colors[i],
            label="Background",
        )
        ax.plot(
            handle.shifted("x"),
            data.y / y0 + i * offset,
            "-",
            color=colors[i],
            label="Raw data (realigned)",
        )

        # Annotate graphs
        val = (
            handle.background("y")[np.isfinite(handle.background("y"))] / y0
            + i * offset
        )
        val = val[get_range(handle.background("x"), 300, 550)]
        val1 = max(val) + 1

        ax.text(
            380,
            val1,
            (
                f'{handle.meta("O16"):.1f}' + r"\,\%"
                if tex
                else f'{handle.meta("O16")/100:.1%}'
            ),
            horizontalalignment="right",
            verticalalignment="bottom",
        )
        ax.text(
            495,
            val1,
            (
                f'{handle.meta("O18"):.1f}' + r"\,\%"
                if tex
                else f'{handle.meta("O18")/100:.1%}'
            ),
            horizontalalignment="left",
            verticalalignment="bottom",
        )
        label, xpos = description[i]
        ax.text(
            xpos,
            min(i * offset + handle.data.y[get_range(handle.data.x, xpos, 800)] / y0)
            - 1,
            label,
            horizontalalignment="left",
            verticalalignment="top",
        )

        i += 1

    ax.set_xlabel("kinetic energy / (eV)")
    ax.set_ylabel("counts / ($10^3$ s$^{-1}$)")
    ax.set_xlim(300, 800)
    ax.set_ylim(0, ylim[sample])

    plt.tight_layout()
    if save:
        plt.savefig(save_path + f"{sample}_custom.png")
        plt.savefig(save_path + f"{sample}_custom.svg")
plt.show()
