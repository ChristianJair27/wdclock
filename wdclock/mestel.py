"""Mestel's analytic cooling law (Mestel 1952, MNRAS 112, 583).

The simplest white-dwarf clock: the ions of the core hold the heat (ideal-gas
heat capacity, 3/2 k per ion), the degenerate electrons conduct it to the base of
a thin radiative envelope, and Kramers opacity in that envelope sets the leak
rate, L proportional to M T_c^3.5. Integrating dE/dt = -L gives the cooling time
as a function of the present luminosity. In the form used in Michael Zingale's
stellar-physics notes (zingale.github.io/stars, "WD Cooling"):

    t_cool = 6e6 yr * (A/12)^-1 * (M/Msun)^(5/7) * (mu/2)^(-2/7)
             * [ (L/Lsun)^(-5/7) - (L0/Lsun)^(-5/7) ]

with A the mean ion mass number (mu_I), mu the mean molecular weight of the
envelope-ish plasma (ions plus electrons) and L0 the starting luminosity. It
ignores neutrino cooling (dominant above ~25,000 K), the Debye drop of the heat
capacity in old stars, crystallisation (latent heat, C/O phase separation) and the
hydrogen-layer thickness, so it is expected to run fast for young, hot stars and to
disagree by a factor of order unity below ~1e-4 Lsun. Comparing it with the Montreal
sequences is the point of ``wdclock.compare``.
"""
from __future__ import annotations

import numpy as np


def mean_ion_mass(x_c: float = 0.5, x_o: float = 0.5, x_ne: float = 0.0) -> float:
    """mu_I for a C/O(/Ne) mixture given mass fractions."""
    return 1.0 / (x_c / 12 + x_o / 16 + x_ne / 20)


def mean_molecular_weight(mu_i: float, ye: float = 0.5) -> float:
    """Ions plus electrons: 1/mu = 1/mu_I + Y_e."""
    return 1.0 / (1.0 / mu_i + ye)


def cooling_age(lum, mass=0.6, lum0=1.0, x_c=0.5, x_o=0.5, ye=0.5):
    """Mestel cooling age in years for luminosity ``lum`` [Lsun] of a ``mass`` [Msun] white dwarf."""
    lum = np.asarray(lum, float)
    mu_i = mean_ion_mass(x_c, x_o)
    mu = mean_molecular_weight(mu_i, ye)
    return 6.0e6 * (12.0 / mu_i) * mass ** (5.0 / 7.0) * (mu / 2.0) ** (-2.0 / 7.0) * (lum ** (-5.0 / 7.0) - lum0 ** (-5.0 / 7.0))
