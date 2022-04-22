import json
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy.optimize import minimize
import pandas as pd

from pyOER.constants import (
    FARADAY_CONSTANT as F,
    STANDARD_TEMPERATURE as T0,
    GAS_CONSTANT as R,
    STANDARD_SYMMETRY_FACTOR as alpha_0,  # 0.5
    STANDARD_SITE_DENSITY,
    STANDARD_SPECIFIC_CAPACITANCE,
)
from part_I_fig3 import plot_all_activity_results
from pyOER.modelling import get_states


try:  # See if the results to fit have already been gathered and saved
    df = pd.read_csv("results_for_fitting.csv")
    u_exp = df["E v RHE / (V)"]
    j_norm_exp = df["norm. activity / (A/F)"]
except FileNotFoundError:  # if not, then gather and save them.
    u_exp, tof_exp = plot_all_activity_results(ax=None, result="tof", for_model=True)
    # clean up the data a bit (remove nan's and negatives if they appear):
    mask = np.logical_not(np.logical_or(np.isnan(tof_exp), tof_exp < 0))
    u_exp = u_exp[mask]
    tof_exp = tof_exp[mask]
    # get the capacitance-normalized current from the tof returned above:
    # ^ FIXME: Should work in capacitance-normalized currents by default not TOF.
    j_norm_exp = tof_exp * STANDARD_SITE_DENSITY / STANDARD_SPECIFIC_CAPACITANCE * 4 * F
    # save it to a .csv file via a pandas dataframe:
    df = pd.DataFrame({"E v RHE / (V)": u_exp, "norm. activity / (A/F)": j_norm_exp})
    df.to_csv("results_for_fitting.csv")


U0 = 1.23  # reference potential vs RHE in [V]


def get_model_j_over_j0(G1, G2, u, alpha=alpha_0):
    """Return the modeled activity normalized to the active state exchange activity

    Args:
        G1 (float): Energy wrt the active state of the state two electrons prior [eV]
        G2 (float): Energy wrt the active state of the state two electrons prior [eV]
        u (Array): The potential data
        alpha (float): The symmetry factor of the rate-determining step

    Returns Array: The model-calculated activity normalized to exchange activity
    """

    # ---- and lets do some math! :D
    states = get_states(G1, G2)

    # This matrix is the coverage of each state (outer dimension) relative to the RDS
    #   state as a function of potential (outer dimension):
    theta_over_theta_i = np.array(
        [
            1 / state.K_rds * np.exp(-state.n_to_rds * F * (u - U0) / (R * T0))
            for state in states
        ]
    )
    # The coverage of the RDS state is the first of these over their sum:
    theta_i = theta_over_theta_i[0] / np.sum(theta_over_theta_i, axis=0)

    # From here it is easy to calculate the current relative to j0 of the RDS
    j_over_j0 = theta_i * np.exp(alpha * F * (u - U0) / (R * T0))
    return j_over_j0


def get_model_j_norm(G1, G2, j_norm_0, u=u_exp, alpha=alpha_0):
    """Return the modeled capacitance-normalized activity in [A/F]

    Args:
        G1 (float): Energy wrt the active state of the state two electrons prior [eV]
        G2 (float): Energy wrt the active state of the state two electrons prior [eV]
        j_norm_0 (float): The capacitance-normalized OER exchange current density of the
            active state in [A/F]
        u (Array): The potential data in [V] vs RHE. Defaults the experimental set.
        alpha (float): The symmetry factor of the rate-determining step. Defaults to 0.5
    Returns Array: The model-calculated activity normalized to capacitance in [A/F]
    """
    j_over_j0 = get_model_j_over_j0(G1=G1, G2=G2, u=u, alpha=alpha)

    return j_over_j0 * j_norm_0


def square_error(params):
    """Return the sum of sq. error on log(j_norm) calc'd with params=(G1, G2, j_norm_0)"""
    G1, G2, j_norm_0 = params
    j_model = get_model_j_norm(G1, G2, j_norm_0, u_exp, alpha_0)
    error = np.log10(j_model) - np.log10(j_norm_exp)
    return error.dot(error)


if __name__ == "__main__":

    forpublication = True
    # plt.interactive(False)  # show the plot when I tell you to show() it!

    if forpublication:  # for the publication figure
        import matplotlib as mpl

        mpl.rcParams["figure.figsize"] = (3.25, 2.75)
        # plt.rc('text', usetex=True)  # crashingly slow
        plt.rc("font", family="sans-serif")
        plt.rc("font", size=8)
        plt.rc("lines", linewidth=0.5)
        plt.rc("lines", markersize=3)
    else:
        plt.style.use("default")

    result = minimize(
        square_error,
        np.array([-0.1, -0.3, 2e-3]),
        bounds=[(None, None), (None, None), (1e-10, None)],
    )

    G1, G2, j_norm_0 = result.x

    with open("fit_results.json", "w") as f:
        json.dump({"G1": G1, "G2": G2, "j_norm_0": j_norm_0}, f)

    u = np.linspace(1.23, 1.55, 100)
    j_over_j0 = get_model_j_over_j0(G1, G2, u)

    gs = gridspec.GridSpec(8, 1)
    # gs.update(hspace=0.025)
    ax1 = plt.subplot(gs[0:4, 0])
    ax2 = plt.subplot(gs[4:6, 0])
    ax3 = ax2.twinx()
    ax4 = plt.subplot(gs[6:8, 0])

    ax1.xaxis.set_label_position("top")
    ax1.tick_params(axis="x", top=True, bottom=True, labeltop=True, labelbottom=False)
    ax2.tick_params(axis="x", top=False, bottom=True, labeltop=False, labelbottom=False)
    ax4.tick_params(axis="x", top=False, bottom=True, labeltop=False, labelbottom=True)
    ax3.spines["right"].set_color("r")
    ax3.tick_params(axis="y", color="r")
    ax3.tick_params(axis="y", labelcolor="r")
    ax3.yaxis.label.set_color("r")

    fig = ax1.get_figure()

    if True:  # plot the experimental results with the colors (a bit slow):

        tof_to_j_norm = STANDARD_SITE_DENSITY * (4 * F) / STANDARD_SPECIFIC_CAPACITANCE
        tof_0 = j_norm_0 / tof_to_j_norm

        plot_all_activity_results(ax=ax1, result="tof", factor=1 / tof_0, takelog=True)
    elif True:  # plot all the experimental results black and white:
        ax1.plot(
            u_exp,
            np.log10(j_norm_exp / j_norm_0),
            "k",
            marker="o",
            linestyle="none",
            fillstyle="none",
        )

    ax1.plot(u, np.log10(j_over_j0), linewidth=0.6, color="k")

    ax1.set_xlabel("E vs RHE / (V)")
    # ax1.set_ylabel("log(j$_{O2, norm}$ / j$_{RDS}^0$)")
    ax1.set_ylabel("log(j / j$^0$)")

    if True:  # explanatory fig

        states = get_states(G1, G2)
        theta_over_theta_i = np.array(
            [
                1 / state.K_rds * np.exp(-state.n_to_rds * F * (u - U0) / (R * T0))
                for state in states
            ]
        )
        # ... and this matrix is the coverage of each individual state:
        thetas = theta_over_theta_i / np.sum(theta_over_theta_i, axis=0)

        # ...and also the expected Tafel slope
        # # from the average number of electron transfers to get to the RDS:
        n_to_rds_vec = np.sum(
            [thetas[j] * state.n_to_rds for j, state in enumerate(states)], axis=0
        )
        tafel_vec = R * T0 * np.log(10) / F / (n_to_rds_vec + 1 / 2) * 1e3

        for j, state in enumerate(states):
            ax2.plot(u, thetas[j], state.color)
        ax3.plot(u, tafel_vec, "r--")
        ax3.set_ylim(bottom=0)
        ax3.set_yticks([0, 30, 60, 90, 120])
        ax3.set_yticklabels(["0", "", "60", "", "120"])

        # ax2.set_xlabel("U$_{RHE}$ / [V]")
        ax2.set_ylabel("coverage")
        ax3.set_ylabel("Tafel / (mV dec$^{-1}$)")

    if True:
        # delta G plot
        u_vec = np.array([1.23, 1.55])
        for state in states:
            n = state.n_to_rds
            G0 = state.eV_1p23_vs_rds
            G_vec = G0 + (u_vec - U0) * n
            ax4.plot(u_vec, G_vec, state.color)

        ax4.set_ylim(top=0.1)

        ax4.set_xlabel("E vs RHE / (V)")
        ax4.set_ylabel("$G^0 - G^0_{RDS}$ / (eV)")

    if True:  # wrong models
        # if there was constant coverage of the RDS state,
        #   the current relative to j0 would simply be:
        j_over_j0_simple = np.exp(alpha_0 * F * (u - U0) / (R * T0))

        # if there were two sites with diff. symmetry factors, we might get something like:
        j_over_j0_early_site = 0.99 * np.exp(0.1 * F * (u - U0) / (R * T0))
        j_over_j0_late_site = 0.01 * np.exp(0.9 * F * (u - U0) / (R * T0))
        j_over_j0_2_simple = j_over_j0_early_site + j_over_j0_late_site

        # lets plot the other two models:
        ax1.plot(u, np.log10(j_over_j0_simple), "k--")
        ax1.plot(u, np.log10(j_over_j0_2_simple * 1e-4), "k:")

    fig.set_figheight(fig.get_figwidth() * 1.5)
    fig.savefig("paper_I_v6_fig5.png")
    fig.savefig("paper_I_v6_fig5.svg")
