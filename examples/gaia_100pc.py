"""Worked example: the Gaia EDR3 white dwarfs within 100 pc.

    python examples/gaia_100pc.py path/to/wd_100pc.csv [--lang en|es] [--out figures]

The input is the Gentile Fusillo et al. (2021) catalogue cut at parallax > 10 mas and
P_WD > 0.75 (16,676 stars), with the VizieR column names (TeffH, loggH, MassH, Gmag,
GMAG = absolute G, BP_RP, and Cartesian x, y, z in pc). Writes the age table next to the
input and the six figures of ``wdclock.plots``. Prints the headline numbers as JSON.
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wdclock.catalog import run, summary  # noqa: E402
from wdclock.plots import Figures  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("catalog")
ap.add_argument("--lang", default="en")
ap.add_argument("--out", default="figures")
ap.add_argument("--ages", default=None, help="reuse an existing output table instead of recomputing")
a = ap.parse_args()

df = pd.read_csv(a.ages) if a.ages else run(pd.read_csv(a.catalog))
if not a.ages:
    out = Path(a.catalog).with_name("wd_ages.csv")
    df.to_csv(out, index=False)
    print("ages written to", out)
s = summary(df)
near = df[df["dist_pc"] < 40]["age_thick"].dropna()
s["median_age_gyr_within_40pc"] = float(near.median())
s["n_within_40pc"] = int(near.size)
print(json.dumps(s, indent=1))

F = Figures(lang=a.lang, outdir=a.out)
for fn in (F.hr_isochrones, F.age_histogram, F.hz_periods, F.map3d_age):
    print(fn(df))
print(F.mestel_vs_montreal())
print(F.hz_duration_map())
