"""Habitable zone of a cooling white dwarf.

Physics
-------
A planet's equilibrium temperature depends on the stellar flux it receives,
S = L / (4 pi a^2). Climate models give the flux limits between which liquid water
can persist on an Earth-like planet: the inner edge where a runaway greenhouse
boils the oceans, the outer edge where even a maximal CO2 greenhouse cannot stop
global glaciation. In units of the solar constant those limits are S_in = 1.0512
and S_out = 0.3438 for a Sun-like spectrum (Kopparapu et al. 2013, ApJ 765, 131,
"runaway greenhouse" and "maximum greenhouse"; the optimistic "recent Venus" and
"early Mars" limits are 1.7753 and 0.3179). Kopparapu's temperature-dependent
corrections are only valid for 2600-7200 K and we do not extrapolate them to white
dwarfs; Agol (2011, ApJL 731, L31) likewise scaled Solar-System limits. The habitable
distance is then a = sqrt(L / S) AU with L in Lsun: for a typical cool white dwarf,
L ~ 1e-4 Lsun, that is a ~ 0.01 AU, orbital periods of 4-32 hours.

Because a white dwarf only cools, its habitable zone marches inward. A planet at a
fixed orbit enters the zone when L has dropped to S_in a^2 and leaves it when L
reaches S_out a^2; the time between is its habitable duration. Agol's
"continuously habitable zone" is the range of a for which that duration exceeds a
minimum, 3 Gyr in his paper (0.005-0.02 AU for 0.4-0.9 Msun). Inside the Roche
limit, a_R = 0.0054 AU (rho_p/rho_earth)^(-1/3) (M/0.6 Msun)^(1/3) (Agol 2011), a
rocky planet is torn apart; for Teff above ~4500 K the zone lies outside it.

Caveats the numbers do not carry: a young white dwarf (Teff > 1e4 K) bathes the
planet in UV that can photolyse water (Barnes & Heller 2013; Kozakis et al. 2018);
the planet must have arrived at 0.01 AU after the giant phases, which requires
migration; tidal locking is immediate. The zone is necessary, not sufficient.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .montreal import CoolingGrid

AU_CM = 1.495978707e13
G_CGS = 6.674e-8
MSUN_G = 1.98847e33
RSUN_CM = 6.957e10
REARTH_CM = 6.371e8
YEAR_S = 3.15576e7
DAY_S = 86400.0

S_LIMITS = {                       # Kopparapu et al. 2013, Table 3, solar-spectrum values
    "runaway_greenhouse": 1.0512,
    "maximum_greenhouse": 0.3438,
    "recent_venus": 1.7753,
    "early_mars": 0.3179,
}


def hz_distance_au(lum_lsun, s_eff: float) -> np.ndarray:
    """Orbital distance [AU] where an L [Lsun] star delivers s_eff solar constants."""
    return np.sqrt(np.asarray(lum_lsun, float) / s_eff)


def period_days(a_au, mass_msun) -> np.ndarray:
    """Keplerian period [days] for a planet of negligible mass."""
    a = np.asarray(a_au, float) * AU_CM
    return 2 * np.pi * np.sqrt(a ** 3 / (G_CGS * np.asarray(mass_msun, float) * MSUN_G)) / DAY_S


def roche_limit_au(mass_msun, rho_ratio: float = 1.0) -> np.ndarray:
    """Agol (2011) eq.: a_R = 0.0054 AU (rho_p/rho_earth)^(-1/3) (M/0.6)^(1/3)."""
    return 0.0054 * rho_ratio ** (-1.0 / 3.0) * (np.asarray(mass_msun, float) / 0.6) ** (1.0 / 3.0)


def transit_probability(a_au, r_wd_rsun, r_p_rearth: float = 1.0) -> np.ndarray:
    """Geometric probability (R_wd + R_p) / a for a circular orbit."""
    return (np.asarray(r_wd_rsun, float) * RSUN_CM + r_p_rearth * REARTH_CM) / (np.asarray(a_au, float) * AU_CM)


def transit_duration_min(a_au, mass_msun, r_wd_rsun, r_p_rearth: float = 1.0) -> np.ndarray:
    """Central-transit duration [minutes] ~ (P/pi) (R_wd + R_p) / a."""
    p = period_days(a_au, mass_msun) * DAY_S
    return p / np.pi * transit_probability(a_au, r_wd_rsun, r_p_rearth) / 60.0


def transit_depth(r_wd_rsun, r_p_rearth: float = 1.0) -> np.ndarray:
    """(R_p / R_wd)^2, capped at 1: an Earth in front of an Earth-sized star is a total eclipse."""
    return np.minimum(1.0, (r_p_rearth * REARTH_CM / (np.asarray(r_wd_rsun, float) * RSUN_CM)) ** 2)


@dataclass
class HabitableZone:
    """Habitable-zone quantities for a white dwarf of given mass along its cooling track."""

    grid: CoolingGrid
    s_in: float = S_LIMITS["runaway_greenhouse"]
    s_out: float = S_LIMITS["maximum_greenhouse"]

    def now(self, teff, mass):
        """Inner/outer edge [AU] and the orbital periods [d] at those edges, today."""
        lum = self.grid.luminosity(teff, mass)
        a_in, a_out = hz_distance_au(lum, self.s_in), hz_distance_au(lum, self.s_out)
        return {"lum": lum, "a_in": a_in, "a_out": a_out,
                "p_in": period_days(a_in, mass), "p_out": period_days(a_out, mass),
                "roche": roche_limit_au(mass)}

    def duration_at(self, a_au, mass: float) -> dict:
        """How long a planet at ``a_au`` around a ``mass`` white dwarf sits in the zone.

        Returns the entry and exit cooling ages [yr] and their difference. NaN exit
        means the track ends (age > ~10-15 Gyr) before the planet leaves the zone.
        """
        tr = self.grid.track(mass)
        age, lum = tr["age"], tr["lum"]
        ok = np.isfinite(age) & np.isfinite(lum)
        age, lum = age[ok], lum[ok]
        a2 = np.asarray(a_au, float) ** 2
        # luminosity falls with age; find when it crosses S_in a^2 (enter) and S_out a^2 (exit)
        enter = _crossing(age, lum, self.s_in * a2)
        exit_ = _crossing(age, lum, self.s_out * a2)
        return {"enter": enter, "exit": exit_, "duration": exit_ - enter}

    def continuously_habitable(self, mass: float, t_min_yr: float = 3e9, n: int = 400):
        """Range of a [AU] with habitable duration >= t_min (Agol's CHZ), and the a of maximum duration."""
        a = np.logspace(-3, -0.7, n)
        d = self.duration_at(a, mass)["duration"]
        good = np.isfinite(d) & (d >= t_min_yr)
        if not good.any():
            return {"a_min": np.nan, "a_max": np.nan, "a_best": a[np.nanargmax(d)], "t_best": np.nanmax(d), "a": a, "duration": d}
        return {"a_min": a[good].min(), "a_max": a[good].max(), "a_best": a[np.nanargmax(d)],
                "t_best": np.nanmax(d), "a": a, "duration": d}


def _crossing(age, lum, target):
    """Age at which the decreasing ``lum`` first drops below ``target`` (vectorised over target)."""
    target = np.atleast_1d(np.asarray(target, float))
    out = np.full(target.shape, np.nan)
    order = np.argsort(lum)                      # increasing lum <-> decreasing age
    l_s, a_s = lum[order], age[order]
    inside = (target > l_s[0]) & (target < l_s[-1])
    out[inside] = np.interp(np.log10(target[inside]), np.log10(l_s), np.log10(np.clip(a_s, 1.0, None)))
    out[inside] = 10 ** out[inside]
    out[target >= l_s[-1]] = 0.0                 # brighter than the track ever was: already inside at t=0
    return out if out.size > 1 else out.reshape(np.shape(np.asarray(target).squeeze()))
