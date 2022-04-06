import numpy as np
from matplotlib import pyplot as plt
from pyOER import Measurement, CalibrationSeries, FARADAY_CONSTANT

m = Measurement.open(53)

m.meas.calibrate_RE(m.RE_vs_RHE)
m.meas.normalize_current(m.A_el)
m.meas.correct_ohmic_drop(m.R_Ohm)

axes_a = m.plot(mass_list=["M4", "M32", "M34", "M36"], tspan=[5, 6500])
fig_a = axes_a[0].get_figure()
fig_a.set_figwidth(fig_a.get_figheight() * 2.5)
fig_a.savefig("fig_S2a.png")

sensitivity_factor = CalibrationSeries.load().F_of_tstamp(m.tstamp)

tspans = {
    1.45: {"He": [3000, 3030], "O2": [5900, 5930]},
    1.44: {"He": [2700, 2730], "O2": [5600, 5630]},
    1.42: {"He": [2400, 2430], "O2": [5300, 5330]},
    1.40: {"He": [2150, 2180], "O2": [5050, 5080]},
    1.38: {"He": [1900, 1930], "O2": [4800, 4830]},
    1.36: {"He": [1650, 1680], "O2": [4550, 4580]},
    1.34: {"He": [1400, 1430], "O2": [4300, 4330]},
}
tspan_bgs = {"He": [1500, 1530], "O2": [4400, 4430]}

fig_b, ax_b = plt.subplots()
ax_b.set_xlabel("U vs RHE / [V]")
ax_b.set_ylabel("(partial) J / [mA cm$^{-2}]$")

fig_c, ax_c = plt.subplots()
ax_c.set_xlabel("U vs RHE / [V]")
ax_c.set_ylabel("Faradaic Eficiency / [%]")

V_list = []
I_He_list = []
I_M36_He_list = []
I_O2_list = []
I_M36_O2_list = []

width = 0.002

for V, tspan_pair in tspans.items():
    tspan_He = tspan_pair["He"]
    tspan_O2 = tspan_pair["O2"]

    n_dot_M34_He = np.mean(m.meas.grab("M34", tspan=tspan_He, tspan_bg=tspan_bgs["He"])[1]) / sensitivity_factor  # average 16O18O flux in [mol/s]
    n_dot_M36_He = np.mean(m.meas.grab("M36", tspan=tspan_He, tspan_bg=tspan_bgs["He"])[1]) / sensitivity_factor  # average 18O2 flux in [mol/s]
    I_He = np.mean(m.meas.grab("raw_current", tspan=tspan_He)[1]) * 1e-3   # average current in [A]
    I_M34_He = n_dot_M34_He * 4 * FARADAY_CONSTANT
    I_M36_He = n_dot_M36_He * 4 * FARADAY_CONSTANT

    n_dot_M34_O2 = np.mean(m.meas.grab("M34", tspan=tspan_O2, tspan_bg=tspan_bgs["O2"])[1]) / sensitivity_factor  # average 16O18O flux in [mol/s]
    n_dot_M36_O2 = np.mean(m.meas.grab("M36", tspan=tspan_O2, tspan_bg=tspan_bgs["O2"])[1]) / sensitivity_factor  # average 18O2 flux in [mol/s]
    I_O2 = np.mean(m.meas.grab("raw_current", tspan=tspan_O2)[1]) * 1e-3   # average current in [A]
    I_M34_O2 = n_dot_M34_O2 * 4 * FARADAY_CONSTANT
    I_M36_O2 = n_dot_M36_O2 * 4 * FARADAY_CONSTANT

    V_list.append(V)
    I_He_list.append(I_He)
    I_O2_list.append(I_O2)
    I_M36_He_list.append(I_M36_He)
    I_M36_O2_list.append(I_M36_O2)

    FE_M36_He = I_M36_He / I_He * 100
    FE_M34_He = I_M34_He / I_He * 100
    ax_c.bar(V, FE_M36_He, width=width, color="g")
    ax_c.bar(V, FE_M34_He, bottom=FE_M36_He, width=width, color="r", hatch="//")

    FE_M36_O2 = I_M36_O2 / I_O2 * 100
    FE_M34_O2 = I_M34_O2 / I_O2 * 100
    ax_c.bar(V + width, FE_M36_O2, width=width, color="g", hatch="//")
    ax_c.bar(V + width, FE_M34_O2, bottom=FE_M36_O2, width=width, color="r", hatch="//")
