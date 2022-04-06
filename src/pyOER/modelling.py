import numpy as np
from .constants import (
    FARADAY_CONSTANT as F,
    STANDARD_TEMPERATURE as T0,
    GAS_CONSTANT as R,
    STANDARD_SYMMETRY_FACTOR,
)


U0 = 1.23  # sort of arbitrary reference potential vs RHE in [V]


class State:
    def __init__(self, n_to_rds, eV_1p23_vs_rds, color):
        """Nice little class to represent a surface state

        States are described relative to the RDS state, which is the state from which
        the rate-limiting step of OER takes place.
        States are described with reference to a standard potential U0 = 1.23 V_RHE

        Args:
            n_to_rds (int): ox. state of RDS intermediate relative to self
            eV_1p23_vs_rds (float): energy relative to RDS at U0 / [eV]
            color (str): default color for plotting
        """
        self.n_to_rds = n_to_rds
        self.eV_1p23_vs_rds = eV_1p23_vs_rds
        self.color = color

    @property
    def G_1p23_vs_rds(self):
        """Energy relative to RDS at U0 / [J/mol]"""
        return self.eV_1p23_vs_rds * F

    @property
    def K_rds(self):
        """equilibrium constant = theta_RDS / theta_self"""
        return np.exp(self.G_1p23_vs_rds / (R * T0))


def get_states(G1, G2):
    # ---------- state definitions --------------- #
    return (
        State(0, 0, "k"),  # the RDS state
        State(1, G1, "b"),  # the state one electron transfer behind the RDS state
        State(2, G2, "g"),  # the state two electron transfers behind the RDS state
    )
