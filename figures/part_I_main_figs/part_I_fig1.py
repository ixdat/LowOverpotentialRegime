from matplotlib import pyplot as plt
import numpy as np

from pyOER import Measurement
from ixdat.constants import FARADAY_CONSTANT

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

plt.close("all")

T = 298.15
A_el = 0.196

measurement = Measurement.open(10)

meas = measurement.meas

meas.plot()

V_str = meas.calibrate_RE(RE_vs_RHE=0.715)
J_str = meas.normalize_current(A_el=0.196)

O2 = meas.ecms_calibration(mol="O2", mass="M32", tspan=[300, 350], n_el=4)
cv = meas.cut(tspan=[6370, 6700], t_zero="start").as_cv()

axes = cv.plot_measurement(
    mol_list=[O2],
    logplot=False,
    tspan_bg=[0, 20],
    # J_str="cycle"
)
axes[0].set_ylabel("O$_2$ / (pmol s$^{-1}$cm$^{-2}_{geo})$")
axes[1].set_ylabel("E vs RHE / (V)")
axes[3].set_ylabel("J / (mA cm$^{-2}_{geo}$)")
axes[0].set_xlabel("time / (s)")
axes[1].set_xlabel("time / (s)")

if forpublication:
    axes[0].get_figure().savefig("paper_I_v5_fig1d.png")
    axes[0].get_figure().savefig("paper_I_v5_fig1d.svg")


cv1 = cv.cut(tspan=[86, 172])

t, v = cv1.grab(V_str)
j = cv1.grab_for_t(J_str, t)

fig, ax = plt.subplots()
ax.plot(v, j, "b")

x, y = cv1.grab_flux(O2, tspan_bg=[90, 100])

if False:  # use the integral and a modeled shape
    n_dot_O2 = np.trapz(y, x)  # in mol
    Q_O2 = n_dot_O2 * 4 * Chem.Far  # in C
    q_O2 = Q_O2 * 1e3 / A_el  # in mC/cm^2
    V_times_J_O2 = 20e-3 * q_O2 / 2  # in V/s * mC/cm^2 = V * mA/cm^2
    v_max = max(v)
    nernst = Chem.R * T / Chem.Far  # 26 mV per factor e  (60 mV per decade)
    if False:  # put the highlight at the top of the CV
        mask_increasing = np.logical_and(v > np.append(v[0] - 1, v[:-1]), v > 1)
        v_O2_model = v[mask_increasing]
        j_ciel = j[mask_increasing]
        j_O2_model = V_times_J_O2 / nernst * np.exp((v_O2_model - v_max) / nernst)
        ax.fill_between(v_O2_model, j_ciel, j_ciel - j_O2_model, color="k", alpha=0.5)
    else:  # put the highlight in the middle of the CV
        v_O2_model = np.linspace(1.23, v_max, 100)
        j_O2_model = V_times_J_O2 / nernst * np.exp((v_O2_model - v_max) / nernst)
        zero_O2_model = np.zeros(v_O2_model.shape)
        ax.fill_between(v_O2_model, j_O2_model, zero_O2_model, color="k", alpha=0.5)
else:  # use the measured signal
    j_O2_model = y * 4 * FARADAY_CONSTANT * 1e3 / A_el  # in mA/cm^2
    dt_O2 = 5  # oxygen delay in s
    v_O2_model = np.interp(x - dt_O2, t, v)
    zero_O2_model = np.zeros(v_O2_model.shape)
    mask_increasing = np.logical_and(
        v_O2_model > np.append(v_O2_model[1:], v_O2_model[-1] - 1), v_O2_model > 1
    )
    mask_decreasing = np.logical_and(
        v_O2_model > np.append(v_O2_model[0] - 1, v_O2_model[:-1]), v_O2_model > 1
    )
    ax.plot(v_O2_model, j_O2_model, "k")
    for mask in mask_increasing, mask_decreasing:
        ax.fill_between(
            v_O2_model[mask],
            j_O2_model[mask],
            zero_O2_model[mask],
            # where=j_O2_model > zero_O2_model,
            color="k",
            alpha=0.2,
        )


ax.set_xlabel("E vs RHE / (V)")
ax.set_ylabel("J / (mA cm$^{-2}_{geo}$)")

if forpublication:
    ax.get_figure().savefig("paper_I_v6_fig1e.png")
    ax.get_figure().savefig("paper_I_v6_fig1e.svg")
