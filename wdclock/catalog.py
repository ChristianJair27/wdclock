"""Run the clock over a catalogue of white dwarfs.

Input: any table with an effective temperature column and either a mass or a
log g column (defaults match the Gaia EDR3 catalogue of Gentile Fusillo et al.
2021 as distributed on VizieR, J/MNRAS/508/3877: ``TeffH``, ``loggH``, ``MassH``;
pure-hydrogen solution). Output: the same rows plus

    age_thick, age_thin   cooling age [Gyr] for q_H = 1e-4 and 1e-10
    age_err_layer         |age_thick - age_thin| / 2, the hydrogen-layer systematic
    lum                   luminosity [Lsun] (thick)
    radius                radius [Rsun] (thick)
    crystallised          crystallised core mass fraction (thick)
    hz_in_au, hz_out_au   habitable-zone edges today [AU] (runaway / maximum greenhouse)
    hz_p_in_d, hz_p_out_d orbital periods at those edges [days]
    hz_roche_au           Roche limit for an Earth-density planet [AU]
    hz_dur_gyr            habitable duration [Gyr] of a planet now at the zone centre
    hz_transit_prob       geometric transit probability for an Earth at the zone centre
    hz_transit_min        central transit duration [min] for that planet
    hz_depth              transit depth (R_p/R_wd)^2 for an Earth, capped at 1

Rows outside the model grid (mass < 0.20 or > 1.30 Msun, Teff outside the tracks,
missing parameters) get NaN. Gentile Fusillo et al. warn that parameters above
40,000 K, below 0.2 or above 1.3 Msun, and for unresolved binaries are unreliable;
the He-atmosphere and mixed solutions are not used here (use ``teff_col``/``mass_col``
to point at them if you want to).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .hz import HabitableZone, hz_distance_au, period_days, roche_limit_au, transit_depth, transit_duration_min, transit_probability
from .montreal import CoolingGrid

GYR = 1e9


def run(df: pd.DataFrame, teff_col: str = "TeffH", mass_col: str | None = "MassH", logg_col: str | None = "loggH",
        thick: CoolingGrid | None = None, thin: CoolingGrid | None = None, hz_duration: bool = True) -> pd.DataFrame:
    thick = thick or CoolingGrid("thick")
    thin = thin or CoolingGrid("thin")
    out = df.copy()
    teff = out[teff_col].to_numpy(float)
    if mass_col and mass_col in out:
        mass = out[mass_col].to_numpy(float)
    else:
        mass = thick.mass_from_logg(teff, out[logg_col].to_numpy(float))
        out["mass_from_logg"] = mass
    out["age_thick"] = thick.age(teff, mass) / GYR
    out["age_thin"] = thin.age(teff, mass) / GYR
    out["age_err_layer"] = np.abs(out["age_thick"] - out["age_thin"]) / 2
    out["lum"] = thick.luminosity(teff, mass)
    out["radius"] = thick.radius(teff, mass)
    out["crystallised"] = thick.crystallised_fraction(teff, mass)
    hz = HabitableZone(thick)
    out["hz_in_au"] = hz_distance_au(out["lum"], hz.s_in)
    out["hz_out_au"] = hz_distance_au(out["lum"], hz.s_out)
    out["hz_p_in_d"] = period_days(out["hz_in_au"], mass)
    out["hz_p_out_d"] = period_days(out["hz_out_au"], mass)
    out["hz_roche_au"] = roche_limit_au(mass)
    a_mid = np.sqrt(out["hz_in_au"] * out["hz_out_au"])
    out["hz_transit_prob"] = transit_probability(a_mid, out["radius"])
    out["hz_transit_min"] = transit_duration_min(a_mid, mass, out["radius"])
    out["hz_depth"] = transit_depth(out["radius"])
    if hz_duration:
        out["hz_dur_gyr"] = _durations(hz, a_mid.to_numpy(), mass) / GYR
    return out


def _durations(hz: HabitableZone, a_mid: np.ndarray, mass: np.ndarray) -> np.ndarray:
    """Habitable duration at a_mid; the track is evaluated once per 0.01-Msun mass bin for speed."""
    out = np.full(a_mid.shape, np.nan)
    ok = np.isfinite(a_mid) & np.isfinite(mass)
    bins = np.round(mass[ok], 2)
    idx_ok = np.flatnonzero(ok)
    for m in np.unique(bins):
        sel = idx_ok[bins == m]
        out[sel] = hz.duration_at(a_mid[sel], float(m))["duration"]
    return out


def summary(res: pd.DataFrame) -> dict:
    """Headline numbers for a report: medians, oldest decile, fraction crystallised, HZ periods."""
    a = res["age_thick"].dropna()
    return {
        "n_input": int(len(res)),
        "n_with_age": int(a.size),
        "median_age_gyr": float(a.median()),
        "p16_age_gyr": float(a.quantile(0.16)),
        "p84_age_gyr": float(a.quantile(0.84)),
        "p90_age_gyr": float(a.quantile(0.90)),
        "max_age_gyr": float(a.max()),
        "median_layer_err_gyr": float(res["age_err_layer"].dropna().median()),
        "frac_crystallised_gt_half": float((res["crystallised"].dropna() > 0.5).mean()),
        "median_hz_p_in_hours": float((res["hz_p_in_d"].dropna() * 24).median()),
        "median_hz_p_out_hours": float((res["hz_p_out_d"].dropna() * 24).median()),
        "median_hz_duration_gyr": float(res["hz_dur_gyr"].dropna().median()) if "hz_dur_gyr" in res else None,
        "median_transit_prob": float(res["hz_transit_prob"].dropna().median()),
    }
