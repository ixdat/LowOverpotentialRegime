import re
from pyOER.tof import all_tofs, TurnOverFrequency
import numpy as np
from matplotlib import pyplot as plt

from pyOER.constants import (
    STANDARD_SPECIFIC_CAPACITANCE,
    STANDARD_SITE_DENSITY,
    FARADAY_CONSTANT,
)


colors = {  # from Reshma's Figure 1.
    "1": "#54bdebff",
    "2": "#308ab2ff",
    "3": "#2e7094ff",
    "4": "#163854ff",
}
colors_soren = {"1": "k", "2": "b", "3": "g"}
# colors = colors_soren
markers = {"16": "o", "18": "s"}
number_finder = "([0-4])"


def get_color(sample_name):
    if "Reshma" in sample_name:
        try:
            T_number = re.search(number_finder, sample_name).group(1)
        except AttributeError:
            return "y"
        else:
            return colors[T_number]
    elif "Evans" in sample_name:
        return "m"
    # elif "Melih" in sample_name:   # not using this anymore
    #     return "m"
    # elif "Rao" in sample_name:   # not using this anymore
    #     return "orange"
    else:
        return "y"  # a sign of a sample that shouldn't be included


def get_sample_type(sample_name):
    if "Reshma" in sample_name:
        return sample_name[:7]
    elif "Evans" in sample_name:
        return "Evans"


def plot_all_activity_results(
    ax=None,
    result="rate",
    factor=1,
    takelog=False,
    for_model=True,
):
    print(f"plotting from {__file__}")
    potential_list = []
    result_list = []
    for tof in all_tofs():
        sample_name = tof.sample_name
        if not (
            tof.tof_type in ("activity", "ec_activity")
            and tof.id > 239
            and (
                "Reshma" in sample_name
                # or "Rao" in sample_name  # not anymore
                or "Evans" in sample_name
                # or "Melih" in sample_name  # not anymore
            )
        ):
            continue
        rate = tof.rate
        potential = tof.potential
        f = tof.tof

        if rate * 1e9 < 3e-5:
            # This is the detection limit.
            continue

        if for_model:
            if (
                (potential > 1.45 and not tof.tof_type == "ec_activity")
                or "Rao" in tof.sample_name
                or "Melih" in tof.sample_name
            ):
                continue
        else:
            if tof.tof_type == "ec_activity":
                continue

        if tof.tof_type == "ec_activity":
            marker = "^"
        else:
            marker = markers.get(tof.measurement.isotope, "o")
        # marker = "o"

        to_plot = None
        if result == "rate":
            result_list.append(rate)
            potential_list.append(potential)
            to_plot = rate * 1e9 * factor
        elif result == "tof":
            result_list.append(f)
            potential_list.append(potential)
            to_plot = f * factor

        if not to_plot:
            print(f"Don't know what {result} is. Not plotting.")
        else:
            if takelog:
                if to_plot <= 0:
                    print(f"Waring! Negative value encountered in {tof}")
                    continue
                to_plot = np.log10(to_plot)
            if ax:
                color = get_color(sample_name)
                ax.plot(
                    potential, to_plot, color=color, marker=marker, fillstyle="none"
                )
    return np.array(potential_list), np.array(result_list)


def plot_aggregated_activity_results(
    ax=None, result="rate", factor=1, takelog=False, for_model=True, separate_ec=True
):
    print(f"plotting from {__file__}")
    aggregated_results = {}
    for tof in all_tofs():
        sample_name = tof.sample_name
        if not (
            tof.tof_type in ("activity", "ec_activity")
            and tof.id > 239
            and (
                "Reshma" in sample_name
                # or "Rao" in sample_name  # not anymore
                or "Evans" in sample_name
                # or "Melih" in sample_name  # not anymore
            )
        ):
            continue
        rate = tof.rate
        potential = np.round(tof.potential, 2)
        f = tof.tof
        sample_type = get_sample_type(sample_name)
        if separate_ec and tof.tof_type == "ec_activity":
            sample_type = "ec_" + sample_type

        if rate * 1e9 < 3e-5:
            # This is the detection limit.
            continue

        if for_model:
            if (
                (potential > 1.45 and not tof.tof_type == "ec_activity")
                or "Rao" in tof.sample_name
                or "Melih" in tof.sample_name
            ):
                continue
        else:
            if tof.tof_type == "ec_activity":
                continue

        if sample_type not in aggregated_results:
            aggregated_results[sample_type] = {}
        if potential not in aggregated_results[sample_type]:
            aggregated_results[sample_type][potential] = []
        aggregated_results[sample_type][potential].append((rate, f))

    for sample_type, sample_results in aggregated_results.items():

        for potential, points in sample_results.items():

            rates, fs = zip(*points)

            rate = np.mean(rates)
            rate_error = np.std(rates)
            f = np.mean(fs)
            f_error = np.std(fs)

            to_plot = None
            error_bars = None
            if result == "rate":
                to_plot = rate * 1e9 * factor
                if len(rates) > 1:
                    error_bars = [
                        to_plot - rate_error * 1e9 * factor,
                        to_plot + rate_error * 1e9 * factor,
                    ]
            elif result == "tof":
                to_plot = f * factor
                if len(fs) > 1:
                    error_bars = [
                        to_plot - f_error * factor,
                        to_plot + f_error * factor,
                    ]

            if sample_type.startswith("ec_"):
                marker = "^"
                sample_type = sample_type[3:]
            else:
                marker = "o"

            if not to_plot:
                print(f"Don't know what {result} is. Not plotting.")
            else:
                if takelog:
                    if to_plot <= 0:
                        print(f"Waring! Negative value encountered in {sample_type}")
                        continue
                    to_plot = np.log10(to_plot)
                if ax:
                    color = get_color(sample_type)
                    ax.plot(
                        potential, to_plot, color=color, marker=marker, fillstyle="none"
                    )
                    if error_bars:
                        ax.plot(
                            [potential, potential], error_bars, color=color, marker="_"
                        )

    return aggregated_results


if __name__ == "__main__":

    # plt.interactive(False)  # show the plot when I tell you to show() it!

    forpublication = True
    if forpublication:  # for the publication figure
        import matplotlib as mpl

        mpl.rcParams["figure.figsize"] = (3.25, 2.75)
        # plt.rc('text', usetex=True)  # crashingly slow
        plt.rc("font", family="sans-serif")
        plt.rc("font", size=8)
        plt.rc("lines", linewidth=0.4)
        plt.rc("lines", markersize=3)
    else:
        plt.style.use("default")

    # ------------ geometric activity ------------ #

    fig1, ax1 = plt.subplots()

    plot_all_activity_results(ax=ax1, result="rate", for_model=False)
    ax1.set_xlabel("E vs RHE / (V)")
    ax1.set_ylabel("rate / (nmol s$^{-1}$)")
    ax1.set_yscale("log")

    if False:  # axis to indicate geometric flux density
        ax1b = ax1.twinx()
        ax1b.set_ylim([lim / 0.196 for lim in ax1.get_ylim()])
        ax1b.set_yscale("log")
        ax1b.set_ylabel("rate / (nmol s$^{-1}$cm$^{-2}_{geo}$)")

    if True:  # a partial current density axis
        ax1b = ax1.twinx()
        j_O2_lim = [lim / 0.196 * 4 * FARADAY_CONSTANT * 1e-6 for lim in ax1.get_ylim()]
        # ^ [nmol/s] * [1/cm^2] * [C/mol] * [mA/nA] -> [mA/cm^2]
        ax1b.set_ylim(j_O2_lim)
        ax1b.set_yscale("log")
        ax1b.set_ylabel("j$_{O2}$ / (mA cm$^{-2}_{geo}$)")

    if forpublication:
        fig1.subplots_adjust(left=0.15, right=0.85)
        fig1.savefig("paper_I_v6_fig3.png")
        fig1.savefig("paper_I_v6_fig3.svg")

    if False:  # other versions of the fig

        # --------- capacitance-normalized activity ------------ #

        fig2, ax2b = plt.subplots()
        ax2 = ax2b.twinx()
        plot_all_activity_results(ax=ax2, result="tof", for_model=False)
        ax2b.set_xlabel("E vs RHE / (V)")

        ax2.set_xlabel("E vs RHE / (V)")
        ax2.set_yscale("log")

        tof_lim = ax2.get_ylim()
        norm_flux_lim = [
            # FIXME: Should get capacitance-normalized current
            #  and then make a TOF axis, not vice versa...
            lim
            * STANDARD_SITE_DENSITY
            / STANDARD_SPECIFIC_CAPACITANCE
            * 4
            * FARADAY_CONSTANT
            for lim in tof_lim
        ]
        #  [1/s] * [mol/cm^2] / [F/cm^2] * [C/mol] = [A/F]
        ax2b.set_ylim(norm_flux_lim)
        ax2b.set_yscale("log")
        ax2b.set_ylabel("j$_{O2, norm}$ / (A F$^{-1}$)")

        if True:  # a TOF axis
            ax2.set_ylabel("TOF / (s$^{-1}$)")
        else:  # no TOF axis
            ax2.set_yticks([])
            ax2.yaxis.set_tick_params(which="both", right=False)
            ax2.set_ylabel("")

        # --------- aggregated activity ------------ #

        fig3, ax3 = plt.subplots()
        plot_aggregated_activity_results(ax=ax3, result="rate", for_model=False)
        ax3.set_yscale("log")
        ax3.set_ylabel("j$_{O2, norm}$ / (A F$^{-1}$)")

        if True:  # a partial current density axis
            ax3b = ax3.twinx()
            j_O2_lim = [
                lim / 0.196 * 4 * FARADAY_CONSTANT * 1e-6 for lim in ax3.get_ylim()
            ]
            # ^ [nmol/s] * [1/cm^2] * [C/mol] * [mA/nA] -> [mA/cm^2]
            ax3b.set_ylim(j_O2_lim)
            ax3b.set_yscale("log")
            ax3b.set_ylabel("j$_{O2}$ / (mA cm$^{-2}_{geo}$)")

        # --------- aggregated capacitance-normalized activity ------------ #

        fig4, ax4 = plt.subplots()
        ax4b = ax4.twinx()
        plot_aggregated_activity_results(ax=ax4b, result="tof", for_model=False)
        ax4b.set_yscale("log")
        ax4b.set_ylim([1e-7, 1])
        ax4b.set_xlabel("E vs RHE / (V)")
        ax4b.set_ylabel("TOF / (s$^{-1}$)")

        tof_lim = ax4b.get_ylim()
        norm_flux_lim = [
            # FIXME: Should get capacitance-normalized current
            #  and then make a TOF axis, not vice versa...
            lim
            * STANDARD_SITE_DENSITY
            / STANDARD_SPECIFIC_CAPACITANCE
            * 4
            * FARADAY_CONSTANT
            for lim in tof_lim
        ]
        #  [1/s] * [mol/cm^2] / [F/cm^2] * [C/mol] = [A/F]
        ax4.set_ylim(norm_flux_lim)
        ax4.set_yscale("log")
        ax4.set_ylabel("j$_{O2, norm}$ / (A F$^{-1}$)")

        if forpublication:
            fig2.subplots_adjust(left=0.15, right=0.85)
            fig2.savefig("paper_I_v6_fig3_norm.png")
            fig2.savefig("paper_I_v6_fig3_norm.svg")

            fig3.subplots_adjust(left=0.15, right=0.85)
            fig3.savefig("paper_I_v6_fig3_agg.png")
            fig3.savefig("paper_I_v6_fig3_agg.svg")

            fig4.subplots_adjust(left=0.15, right=0.85)
            fig4.savefig("paper_I_v6_fig3_agg_norm.png")
            fig4.savefig("paper_I_v6_fig3_agg_norm.svg")
