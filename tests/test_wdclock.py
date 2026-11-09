"""Checks against stars with independent, published cooling ages and against Agol (2011).

Sirius B: Teff = 25,369 K, M = 1.018 Msun, cooling age ~126 Myr (Bond et al. 2017, ApJ 840, 70,
arXiv:1703.10625; dynamical mass from HST astrometry). 40 Eridani B: Teff ~17,200 K, M = 0.573 Msun
(Mason et al. 2017), cooling age ~122 Myr (Bond et al. 2017b, arXiv:1709.00478, who favour a thin
hydrogen layer). Agol (2011, ApJL 731, L31): continuously habitable zone for >= 3 Gyr at 0.005-0.02 AU
for 0.4-0.9 Msun, maximum habitable duration ~8 Gyr near 0.01 AU with the Bergeron et al. (2001) tracks.
"""
import numpy as np
import pytest

from wdclock import CoolingGrid, HabitableZone, mestel_age


@pytest.fixture(scope="module")
def thick():
    return CoolingGrid("thick")


@pytest.fixture(scope="module")
def thin():
    return CoolingGrid("thin")


def test_sirius_b_age(thick):
    age = thick.age(25369, 1.018)[0] / 1e6
    assert abs(age - 126) / 126 < 0.10


def test_40_eri_b_age(thick, thin):
    ages = [g.age(17200, 0.573)[0] / 1e6 for g in (thick, thin)]
    assert min(abs(a - 122) / 122 for a in ages) < 0.15


def test_age_monotonic_in_teff(thick):
    teff = np.array([20000, 10000, 6000, 4500])
    ages = thick.age(teff, 0.6)
    assert np.all(np.diff(ages) > 0)


def test_logg_inversion_roundtrip(thick):
    m = thick.mass_from_logg(8000, thick.logg(8000, 0.7)[0])[0]
    assert abs(m - 0.7) < 0.01


def test_outside_grid_is_nan(thick):
    assert np.isnan(thick.age(6000, 1.5)[0])
    assert np.isnan(thick.age(1000, 0.6)[0])


def test_chz_matches_agol(thick):
    r = HabitableZone(thick).continuously_habitable(0.6, 3e9)
    assert 0.004 < r["a_min"] < 0.010
    assert 0.014 < r["a_max"] < 0.025
    assert 4e9 < r["t_best"] < 9e9
    assert 0.008 < r["a_best"] < 0.014


def test_hz_period_order_of_hours(thick):
    hz = HabitableZone(thick).now(6000, 0.6)
    assert 0.3 < hz["p_in"][0] < 1.0 and 1.0 < hz["p_out"][0] < 3.0


def test_mestel_within_factor_two_before_crystallisation(thick):
    tr = thick.track(0.6)
    i = np.nanargmin(np.abs(np.log10(tr["lum"]) - (-3)))
    ratio = tr["age"][i] / mestel_age(1e-3, 0.6)
    assert 0.5 < ratio < 2.0
