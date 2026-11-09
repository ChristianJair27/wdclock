"""Montreal white-dwarf evolutionary (cooling) sequences.

Physics in one paragraph
------------------------
A white dwarf no longer burns anything. It is a degenerate carbon-oxygen core
(about the mass of the Sun in a volume the size of the Earth) wrapped in a thin
non-degenerate envelope of helium and, usually, hydrogen. Its luminosity is the
thermal energy of the ions leaking out through that envelope, so the star cools
and dims along a nearly fixed radius: the effective temperature and luminosity
of a white dwarf of known mass are a clock. Reading the clock needs a model of
how fast heat escapes (envelope opacity, hydrogen-layer thickness), how much heat
the core stores (ionic heat capacity, Debye cooling) and the latent heat and
phase separation released when the core crystallises after a few Gyr. The
Montreal sequences of Bedard, Bergeron, Brassard & Fontaine (2020, ApJ 901, 93)
solve all of that for C/O cores with a He mantle (q_He = 1e-2) and either a
"thick" (q_H = 1e-4) or "thin" (q_H = 1e-10) hydrogen layer, for 23 masses from
0.20 to 1.30 Msun. They are the de facto standard for white-dwarf ages; the Gaia
EDR3 catalogue of Gentile Fusillo et al. (2021) used them for its masses, so our
input parameters and our clock share the same physics.

Data source
-----------
https://www.astro.umontreal.ca/~bergeron/CoolingModels/ (files
``seq_MMM_{thick,thin}.txt``, bundle ``AllSequences.tar.gz``). The authors ask
users to acknowledge the web site and cite Bedard et al. (2020). The files are
not redistributed in this repository; ``wdclock.montreal.download()`` fetches
them into a cache directory.

File format: a 4-line header, then each model occupies three lines holding 16
numbers: Mod, Teff [K], log g [cgs], R [cm], Age [yr], L [erg/s], log Tc,
log Pc, log rho_c, Mx/M (crystallised mass fraction), log qx, L_nu [erg/s],
log(H/*), log(He/*), log(C/*), log(O/*).
"""
from __future__ import annotations

import os
import tarfile
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

LSUN = 3.828e33          # erg/s (IAU 2015 nominal)
RSUN = 6.957e10          # cm (IAU 2015 nominal)
MASSES = np.round(np.arange(0.20, 1.3001, 0.05), 2)
BASE_URL = "https://www.astro.umontreal.ca/~bergeron/CoolingModels/CoolingModels/"
COLUMNS = ["mod", "teff", "logg", "radius", "age", "lum", "log_tc", "log_pc", "log_rhoc",
           "mx", "log_qx", "lnu", "log_h", "log_he", "log_c", "log_o"]


def default_cache() -> Path:
    env = os.environ.get("WDCLOCK_MODELS")
    if env:
        return Path(env)
    return Path.home() / ".wdclock" / "montreal"


def download(dest: Path | None = None) -> Path:
    """Fetch AllSequences.tar.gz from the Montreal site and unpack it (about 0.6 MB)."""
    dest = Path(dest or default_cache())
    dest.mkdir(parents=True, exist_ok=True)
    tgz = dest / "AllSequences.tar.gz"
    if not tgz.exists():
        urllib.request.urlretrieve(BASE_URL + "AllSequences.tar.gz", tgz)
    with tarfile.open(tgz) as tf:
        tf.extractall(dest, filter="data")
    return dest


def read_sequence(path: Path) -> np.ndarray:
    """Return a structured array (one row per model) from one seq_*.txt file."""
    nums: list[float] = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("=") or s.startswith("#") or "Teff" in s or "Log(" in s or "Lnu" in s:
                continue
            nums.extend(float(x) for x in s.split())
    arr = np.array(nums).reshape(-1, 16)
    return np.rec.fromarrays(arr.T, names=COLUMNS)


@dataclass
class CoolingGrid:
    """All 23 mass sequences for one hydrogen-layer thickness.

    ``age(teff, mass)`` interpolates log(age) in log(Teff) along each mass track,
    then linearly in mass between the two bracketing tracks. The same scheme gives
    luminosity, radius and log g. Teff is the natural clock variable because it
    decreases monotonically along a cooling track after the first ~1e5 yr.
    """

    layer: str = "thick"
    path: Path | None = None
    tracks: dict[float, np.ndarray] = field(default_factory=dict, repr=False)

    def __post_init__(self):
        root = Path(self.path or default_cache())
        if not any(root.glob("seq_*_%s.txt" % self.layer)):
            download(root)
        for m in MASSES:
            f = root / f"seq_{int(round(m * 100)):03d}_{self.layer}.txt"
            t = read_sequence(f)
            # keep the monotonic cooling part: drop the first rows while Teff still rises
            i0 = int(np.argmax(t["teff"]))
            t = t[i0:]
            order = np.argsort(-t["teff"])          # Teff decreasing
            self.tracks[m] = t[order]

    # ---- single-track interpolation ------------------------------------------------
    def _on_track(self, m: float, teff: np.ndarray, col: str, log: bool = True) -> np.ndarray:
        t = self.tracks[m]
        x = np.log10(t["teff"])[::-1]            # increasing for np.interp
        y = t[col][::-1]
        if log:
            y = np.log10(np.clip(y, 1e-30, None))
        out = np.interp(np.log10(teff), x, y, left=np.nan, right=np.nan)
        return 10 ** out if log else out

    def _interp(self, teff, mass, col, log=True):
        teff = np.atleast_1d(np.asarray(teff, float))
        mass = np.atleast_1d(np.asarray(mass, float))
        teff, mass = np.broadcast_arrays(teff, mass)
        out = np.full(teff.shape, np.nan)
        ok = np.isfinite(teff) & np.isfinite(mass) & (mass >= MASSES[0]) & (mass <= MASSES[-1])
        if not ok.any():
            return out
        idx = np.clip(np.searchsorted(MASSES, mass[ok], side="right") - 1, 0, len(MASSES) - 2)
        m0, m1 = MASSES[idx], MASSES[idx + 1]
        w = (mass[ok] - m0) / (m1 - m0)
        y0 = np.empty(ok.sum()); y1 = np.empty(ok.sum())
        for k, m in enumerate(MASSES):
            s0 = idx == k
            if s0.any():
                y0[s0] = self._on_track(m, teff[ok][s0], col, log)
            s1 = idx + 1 == k
            if s1.any():
                y1[s1] = self._on_track(m, teff[ok][s1], col, log)
        if log:
            val = 10 ** ((1 - w) * np.log10(y0) + w * np.log10(y1))
        else:
            val = (1 - w) * y0 + w * y1
        out[ok] = val
        return out

    def age(self, teff, mass):
        """Cooling age in years for effective temperature(s) [K] and mass(es) [Msun]."""
        return self._interp(teff, mass, "age")

    def luminosity(self, teff, mass):
        """Luminosity in Lsun."""
        return self._interp(teff, mass, "lum") / LSUN

    def radius(self, teff, mass):
        """Radius in Rsun."""
        return self._interp(teff, mass, "radius") / RSUN

    def logg(self, teff, mass):
        return self._interp(teff, mass, "logg", log=False)

    def crystallised_fraction(self, teff, mass):
        return self._interp(teff, mass, "mx", log=False)

    def mass_from_logg(self, teff, logg):
        """Invert log g(Teff, M) for the mass. Vectorised; NaN outside the grid."""
        teff = np.atleast_1d(np.asarray(teff, float))
        logg = np.atleast_1d(np.asarray(logg, float))
        teff, logg = np.broadcast_arrays(teff, logg)
        gg = np.array([self._on_track(m, teff, "logg", log=False) for m in MASSES])  # (23, n)
        out = np.full(teff.shape, np.nan)
        for j in range(teff.size):
            col = gg[:, j]
            good = np.isfinite(col)
            if good.sum() < 2 or not (np.nanmin(col) <= logg.flat[j] <= np.nanmax(col)):
                continue
            out.flat[j] = np.interp(logg.flat[j], col[good], MASSES[good])
        return out

    def track(self, mass: float) -> dict[str, np.ndarray]:
        """Interpolated cooling track at an arbitrary mass: age [yr], teff [K], lum [Lsun], radius [Rsun]."""
        base = self.tracks[MASSES[np.argmin(np.abs(MASSES - mass))]]
        teff = base["teff"]
        return {"teff": teff, "age": self.age(teff, mass), "lum": self.luminosity(teff, mass),
                "radius": self.radius(teff, mass), "mx": self.crystallised_fraction(teff, mass)}
