"""Command line: ``wdclock star``, ``wdclock catalog``, ``wdclock chz``.

    wdclock star --teff 25369 --mass 1.018          # Sirius B
    wdclock star --teff 6000 --logg 8.0              # mass from log g, then the clock
    wdclock catalog in.csv out.csv [--teff-col TeffH --mass-col MassH]
    wdclock chz --mass 0.6 --tmin 3                  # Agol's continuously habitable zone
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np

from .catalog import run, summary
from .hz import HabitableZone
from .montreal import CoolingGrid


def main(argv=None):
    ap = argparse.ArgumentParser(prog="wdclock", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("star"); s.add_argument("--teff", type=float, required=True)
    s.add_argument("--mass", type=float); s.add_argument("--logg", type=float)
    c = sub.add_parser("catalog"); c.add_argument("inp"); c.add_argument("out")
    c.add_argument("--teff-col", default="TeffH"); c.add_argument("--mass-col", default="MassH"); c.add_argument("--logg-col", default="loggH")
    c.add_argument("--no-duration", action="store_true", help="skip the habitable-duration column (faster)")
    z = sub.add_parser("chz"); z.add_argument("--mass", type=float, default=0.6); z.add_argument("--tmin", type=float, default=3.0, help="Gyr")
    a = ap.parse_args(argv)

    if a.cmd == "star":
        thick, thin = CoolingGrid("thick"), CoolingGrid("thin")
        mass = a.mass if a.mass is not None else float(thick.mass_from_logg(a.teff, a.logg)[0])
        hz = HabitableZone(thick).now(a.teff, mass)
        res = {"teff_K": a.teff, "mass_msun": round(mass, 3),
               "age_thick_gyr": round(float(thick.age(a.teff, mass)[0]) / 1e9, 4),
               "age_thin_gyr": round(float(thin.age(a.teff, mass)[0]) / 1e9, 4),
               "lum_lsun": float(hz["lum"][0]), "radius_rsun": float(thick.radius(a.teff, mass)[0]),
               "hz_in_au": float(hz["a_in"][0]), "hz_out_au": float(hz["a_out"][0]),
               "hz_period_in_hours": float(hz["p_in"][0] * 24), "hz_period_out_hours": float(hz["p_out"][0] * 24),
               "roche_au": float(hz["roche"])}
        print(json.dumps(res, indent=1))
    elif a.cmd == "catalog":
        import pandas as pd
        df = pd.read_csv(a.inp)
        res = run(df, a.teff_col, a.mass_col if a.mass_col in df else None, a.logg_col, hz_duration=not a.no_duration)
        res.to_csv(a.out, index=False)
        print(json.dumps(summary(res), indent=1))
    elif a.cmd == "chz":
        hz = HabitableZone(CoolingGrid("thick"))
        r = hz.continuously_habitable(a.mass, a.tmin * 1e9)
        print(json.dumps({"mass": a.mass, "t_min_gyr": a.tmin, "a_min_au": r["a_min"], "a_max_au": r["a_max"],
                          "a_longest_au": r["a_best"], "longest_duration_gyr": r["t_best"] / 1e9}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
