import numpy as np


def calc_OER_rate(experiment, tspan, mol=None):
    """Return the total average flux of O2 in [mol/s] in the experiment over tspan"""
    rate = 0
    mol_list = [mol] if mol else experiment.mol_list
    for mol in mol_list:
        x, y = experiment.calc_flux(mol, tspan=tspan)
        rate += np.mean(y)  # mean flux in [mol/s]

    return rate


def calc_dissolution_rate(experiment, tspan, t_electrolysis=None):
    """Return the average dissolution rate during tspan in [mol/s]

    Args:
        experiment (StandardExperiment): the standard experiment
        tspan (timespan): The time interval for which to get average dissolution rate
        t_electrolysis (float): if given, then assume dissolution only occurs over that
            length of time. I.e., divide the icpms_point amount by t_electrolysis.
    """
    t_vec, n_vec = experiment.get_dissolution_points()

    i_before = int(np.argmax(t_vec > tspan[0])) - 1
    i_after = int(np.argmax(t_vec > tspan[-1]))

    t_interval = tspan[-1] - tspan[0]
    if t_electrolysis:
        if i_after > i_before + 1:
            raise TypeError(
                f"Can't adjust '{experiment}' at tspan={tspan} for t_electrolysis"
                f"because the experiment has icpms samples taken during tspan"
            )
        n_during_interval = n_vec[i_after] * t_interval / t_electrolysis
    else:
        if i_after > i_before + 1:
            n_during_interval = (
                n_vec[i_before + 1]
                * (t_vec[i_before + 1] - tspan[0])
                / (t_vec[i_before + 1] - t_vec[i_before])
            )
            for i in range(i_before + 2, i_after):
                n_during_interval += n_vec[i]
            n_during_interval += (
                n_vec[i_after]
                * (tspan[-1] - t_vec[i_after - 1])
                / (t_vec[i_after] - t_vec[i_after - 1])
            )
        else:
            n_during_interval = (
                n_vec[i_after] * t_interval / (t_vec[i_after] - t_vec[i_after - 1])
            )

    return n_during_interval / t_interval


def calc_exchange_rate(experiment, tspan):
    """Return the average rate of lattice O incorporation in O2 in [mol/s] over tspan"""
    beta = experiment.beta
    x_32, y_32 = experiment.calc_flux("O2_M32", tspan=tspan)
    x_34, y_34 = experiment.calc_flux("O2_M34", tspan=tspan)
    return np.mean(y_34) - np.mean(y_32) * beta


def calc_potential(experiment, tspan):
    """Return the average potential vs RHE in [V] during the experiment over tspan"""
    t, U = experiment.meas.grab("potential", tspan=tspan)
    return np.mean(U)


def calc_current(experiment, tspan):
    """Return the average current in [A] during the experiment over tspan"""
    t, I = experiment.meas.grab(experiment.meas.I_name, tspan=tspan)
    I *= 1e-3  # [mA] -> [A]
    tspan_bg = experiment.tspan_bg_current
    if tspan_bg:
        t_bg, I_bg = experiment.meas.grab(experiment.meas.I_name, tspan=tspan_bg)
        I_bg *= 1e-3  # [mA] -> [A]
        I_bg = np.mean(I_bg)
    elif experiment.tspan_cap:
        I_bg = experiment.calc_background_current()
    else:
        I_bg = 0
    return np.mean(I) - I_bg
