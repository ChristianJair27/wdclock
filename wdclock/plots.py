"""Figures for the clock: dark background, provenance line on every plot, English or Spanish labels.

    from wdclock.plots import Figures
    F = Figures(lang="en", outdir="figures")
    F.hr_isochrones(df)        # Gaia colour-magnitude diagram with cooling isochrones from Table_DA
    F.age_histogram(df)        # thick vs thin ages of the sample
    F.mestel_vs_montreal()     # the analytic law against the full sequences
    F.hz_periods(df)           # where a habitable planet would transit (period window per star)
    F.hz_duration_map()        # habitable duration vs orbital distance for several masses
    F.map3d_age(df)            # the solar neighbourhood coloured by age
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, LogNorm

from .hz import HabitableZone
from .mestel import cooling_age as mestel_age
from .montreal import CoolingGrid, LSUN, default_cache

BG = "#070b14"
AGE_CMAP = LinearSegmentedColormap.from_list("age", ["#9fd3ff", "#ffffff", "#ffd27a", "#ff7a3d", "#7a1f1f"])
STAR_CMAP = LinearSegmentedColormap.from_list("star", ["#ff6a3d", "#ffd27a", "#ffffff", "#9fd3ff", "#6f8dff"])

T = {
    "en": dict(prov="Data: Gaia EDR3 white dwarfs (Gentile Fusillo+ 2021) · models: Montreal cooling sequences (Bédard+ 2020) · wdclock",
               hr_title="The white-dwarf cooling sequence is a clock", hr_x="Gaia BP − RP (colour)", hr_y="Absolute G magnitude",
               age="Cooling age (Gyr)", iso="{:g} Gyr", hist_title="Cooling ages of {n:,} white dwarfs within 100 pc",
               thick="thick H layer (q_H = 1e-4)", thin="thin H layer (q_H = 1e-10)", count="Number of white dwarfs",
               mestel_title="Mestel's 1952 law against the full models (0.6 Msun)", lum="Luminosity (Lsun)", age_y="Cooling age (Gyr)",
               mestel="Mestel analytic", montreal="Montreal sequence", cryst="core crystallising",
               hz_title="Where a habitable planet would transit: orbital period of the habitable zone today",
               hz_x="Orbital period (hours)", hz_y="Effective temperature (K)", inner="inner edge (runaway greenhouse)", outer="outer edge (maximum greenhouse)",
               dur_title="How long a planet stays in the habitable zone as the star cools", dur_x="Orbital distance (AU)", dur_y="Habitable duration (Gyr)",
               roche="Roche limit (Earth density)", chz="≥ 3 Gyr", map_title="The solar neighbourhood's dead stars, coloured by cooling age (100 pc)",
               pc="pc"),
    "es": dict(prov="Datos: enanas blancas de Gaia EDR3 (Gentile Fusillo+ 2021) · modelos: secuencias de enfriamiento de Montreal (Bédard+ 2020) · wdclock",
               hr_title="La secuencia de enfriamiento de las enanas blancas es un reloj", hr_x="Gaia BP − RP (color)", hr_y="Magnitud absoluta G",
               age="Edad de enfriamiento (Gyr)", iso="{:g} Gyr", hist_title="Edades de enfriamiento de {n:,} enanas blancas a menos de 100 pc",
               thick="capa de H gruesa (q_H = 1e-4)", thin="capa de H fina (q_H = 1e-10)", count="Número de enanas blancas",
               mestel_title="La ley de Mestel (1952) frente a los modelos completos (0,6 Msol)", lum="Luminosidad (Lsol)", age_y="Edad de enfriamiento (Gyr)",
               mestel="Mestel analítico", montreal="Secuencia de Montreal", cryst="núcleo cristalizando",
               hz_title="Dónde transitaría un planeta habitable: periodo orbital de la zona habitable hoy",
               hz_x="Periodo orbital (horas)", hz_y="Temperatura efectiva (K)", inner="borde interior (invernadero desbocado)", outer="borde exterior (invernadero máximo)",
               dur_title="Cuánto tiempo sigue un planeta en la zona habitable mientras la estrella se enfría", dur_x="Distancia orbital (UA)", dur_y="Duración habitable (Gyr)",
               roche="límite de Roche (densidad terrestre)", chz="≥ 3 Gyr", map_title="Las estrellas muertas del vecindario solar, coloreadas por edad (100 pc)",
               pc="pc"),
}


def style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "axes.edgecolor": "#3a4256",
        "axes.labelcolor": "#dfe6f3", "xtick.color": "#aab4c8", "ytick.color": "#aab4c8", "text.color": "#eef2fa",
        "grid.color": "#1e2638", "axes.grid": True, "grid.linewidth": 0.5, "font.family": "DejaVu Sans",
        "axes.titleweight": "semibold", "axes.titlesize": 13, "legend.frameon": False,
    })


def read_table_da(path: Path | None = None) -> pd.DataFrame:
    """Montreal synthetic photometry for pure-H atmospheres (Table_DA): Teff, log g, mass, Gaia mags, age."""
    path = Path(path or default_cache()) / "Table_DA"
    rows = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) == 44 and s[0].replace(".", "").isdigit():
                rows.append([float(x) for x in s])
    cols = ["teff", "logg", "mass", "mbol", "bc", "U", "B", "V", "R", "I", "J", "H", "Ks", "Y", "J2", "H2", "K", "W1", "W2", "W3", "W4",
            "S36", "S45", "S58", "S80", "u", "g", "r", "i", "z", "gP", "rP", "iP", "zP", "yP", "G2", "G2_BP", "G2_RP", "G3", "G3_BP", "G3_RP", "FUV", "NUV", "age"]
    return pd.DataFrame(rows, columns=cols)


class Figures:
    def __init__(self, lang: str = "en", outdir: str | Path = "figures", grid: CoolingGrid | None = None):
        self.t = T[lang]
        self.lang = lang
        self.out = Path(outdir); self.out.mkdir(parents=True, exist_ok=True)
        self.grid = grid or CoolingGrid("thick")
        style()

    def _finish(self, fig, name, extra=""):
        fig.text(0.01, 0.005, self.t["prov"] + (" · " + extra if extra else ""), fontsize=7.5, color="#8a94a8", ha="left", va="bottom")
        p = self.out / f"{name}.png"
        fig.savefig(p, dpi=170, bbox_inches="tight")
        plt.close(fig)
        return p

    # 1 -------------------------------------------------------------------------------
    def hr_isochrones(self, df: pd.DataFrame, isochrones=(0.1, 0.5, 1, 2, 4, 6, 8, 10), name="hr_isochrones"):
        tab = read_table_da()
        fig, ax = plt.subplots(figsize=(9.5, 8))
        ok = df["age_thick"].notna()
        sc = ax.scatter(df.loc[ok, "BP_RP"], df.loc[ok, "GMAG"], c=df.loc[ok, "age_thick"], cmap=AGE_CMAP, s=2.5, linewidths=0,
                        norm=LogNorm(vmin=0.03, vmax=10), alpha=0.9, rasterized=True)
        ax.scatter(df.loc[~ok, "BP_RP"], df.loc[~ok, "GMAG"], c="#3a4256", s=1.5, linewidths=0)
        # isochrones: for each age, take masses 0.45-1.05 Msun, read Teff and log g at that age from the cooling
        # sequences, then look up Gaia colour and magnitude in Table_DA (regular Teff x log g grid)
        from scipy.interpolate import RegularGridInterpolator
        tt = np.sort(tab["teff"].unique()); gg = np.sort(tab["logg"].unique())
        piv = {c: tab.pivot_table(index="teff", columns="logg", values=c).reindex(index=tt, columns=gg).to_numpy() for c in ("G3", "G3_BP", "G3_RP")}
        f_g = RegularGridInterpolator((np.log10(tt), gg), piv["G3"], bounds_error=False)
        f_c = RegularGridInterpolator((np.log10(tt), gg), piv["G3_BP"] - piv["G3_RP"], bounds_error=False)
        for age in isochrones:
            pts = []
            for m in np.arange(0.45, 1.051, 0.05):
                tr = self.grid.track(float(m)); ok = np.isfinite(tr["age"]) & (tr["age"] > 0)
                a_, t_ = tr["age"][ok], tr["teff"][ok]
                if not (a_.min() < age * 1e9 < a_.max()):
                    continue
                teff = np.interp(age * 1e9, a_, t_)
                lg = float(self.grid.logg(teff, m)[0])
                pts.append((f_c([[np.log10(teff), lg]])[0], f_g([[np.log10(teff), lg]])[0]))
            pts = np.array([p for p in pts if np.all(np.isfinite(p))])
            if len(pts) > 1:
                ax.plot(pts[:, 0], pts[:, 1], color="#ffffff", lw=0.9, alpha=0.6)
                ax.annotate(self.t["iso"].format(age), (pts[0, 0], pts[0, 1]), fontsize=8, color="#ffffff", alpha=0.85,
                            xytext=(-30, -6), textcoords="offset points")
        ax.invert_yaxis(); ax.set_xlim(-0.7, 2.1); ax.set_ylim(17.5, 8.0)
        ax.set_xlabel(self.t["hr_x"]); ax.set_ylabel(self.t["hr_y"]); ax.set_title(self.t["hr_title"])
        cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02); cb.set_label(self.t["age"]); cb.outline.set_visible(False)
        cb.set_ticks([0.03, 0.1, 0.3, 1, 3, 10]); cb.set_ticklabels(["0.03", "0.1", "0.3", "1", "3", "10"])
        return self._finish(fig, name, "isochrones: 0.45-1.05 Msun, colours from Montreal Table_DA")

    # 2 -------------------------------------------------------------------------------
    def age_histogram(self, df: pd.DataFrame, name="age_histogram"):
        fig, ax = plt.subplots(figsize=(10, 5.5))
        bins = np.arange(0, 13.5, 0.25)
        a, b = df["age_thick"].dropna(), df["age_thin"].dropna()
        ax.hist(a, bins=bins, color="#9fd3ff", alpha=0.85, label=self.t["thick"] + f" (median {a.median():.2f} Gyr)")
        ax.hist(b, bins=bins, histtype="step", color="#ffd27a", lw=1.4, label=self.t["thin"] + f" (median {b.median():.2f} Gyr)")
        ax.set_xlabel(self.t["age"]); ax.set_ylabel(self.t["count"]); ax.set_title(self.t["hist_title"].format(n=int(a.size)))
        ax.legend(loc="upper right")
        return self._finish(fig, name)

    # 3 -------------------------------------------------------------------------------
    def mestel_vs_montreal(self, mass=0.6, name="mestel_vs_montreal"):
        tr = self.grid.track(mass)
        ok = np.isfinite(tr["lum"]) & np.isfinite(tr["age"]) & (tr["age"] > 1e6)
        L, age, mx = tr["lum"][ok], tr["age"][ok] / 1e9, tr["mx"][ok]
        fig, ax = plt.subplots(figsize=(9.5, 6))
        ax.plot(L, age, color="#9fd3ff", lw=2.2, label=self.t["montreal"])
        Lm = np.logspace(-5, -1, 200)
        ax.plot(Lm, mestel_age(Lm, mass) / 1e9, color="#ffd27a", lw=1.6, ls="--", label=self.t["mestel"])
        cr = mx > 0.01
        if cr.any():
            ax.fill_betweenx([0.01, 20], L[cr].min(), L[cr].max(), color="#ff7a3d", alpha=0.12, label=self.t["cryst"])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_xaxis(); ax.set_ylim(0.01, 20)
        ax.set_xlabel(self.t["lum"]); ax.set_ylabel(self.t["age_y"]); ax.set_title(self.t["mestel_title"]); ax.legend(loc="upper left")
        return self._finish(fig, name, "Mestel form as in M. Zingale's stellar-physics notes")

    # 4 -------------------------------------------------------------------------------
    def hz_periods(self, df: pd.DataFrame, gmax=17.0, name="hz_periods"):
        d = df[(df["Gmag"] < gmax) & df["hz_p_in_d"].notna()]
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hlines(d["TeffH"], d["hz_p_in_d"] * 24, d["hz_p_out_d"] * 24, colors="#9fd3ff", alpha=0.08, lw=0.8)
        ax.scatter(d["hz_p_in_d"] * 24, d["TeffH"], s=3, color="#ff7a3d", linewidths=0, label=self.t["inner"], alpha=0.8)
        ax.scatter(d["hz_p_out_d"] * 24, d["TeffH"], s=3, color="#6f8dff", linewidths=0, label=self.t["outer"], alpha=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(self.t["hz_x"]); ax.set_ylabel(self.t["hz_y"]); ax.set_title(self.t["hz_title"])
        ax.set_yticks([4000, 6000, 10000, 20000, 40000, 80000]); ax.set_yticklabels(["4000", "6000", "10 000", "20 000", "40 000", "80 000"]); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xticks([1, 3, 10, 30, 100, 300, 1000, 3000]); ax.set_xticklabels(["1", "3", "10", "30", "100", "300", "1000", "3000"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.legend(loc="lower right", markerscale=4)
        return self._finish(fig, name, f"{len(d):,} stars with G < {gmax:g} (TESS reach)")

    # 5 -------------------------------------------------------------------------------
    def hz_duration_map(self, masses=(0.4, 0.6, 0.8, 1.0), name="hz_duration"):
        hz = HabitableZone(self.grid)
        fig, ax = plt.subplots(figsize=(9.5, 6))
        colors = ["#9fd3ff", "#ffffff", "#ffd27a", "#ff7a3d"]
        for m, c in zip(masses, colors):
            r = hz.continuously_habitable(m, 3e9)
            ax.plot(r["a"], r["duration"] / 1e9, color=c, lw=2, label=f"{m:.1f} Msun")
            ax.axvline(0.0054 * (m / 0.6) ** (1 / 3), color=c, lw=0.8, ls=":", alpha=0.6)
        ax.axhline(3, color="#aab4c8", lw=0.8, ls="--"); ax.text(0.0011, 3.15, self.t["chz"], fontsize=8, color="#aab4c8")
        ax.set_xscale("log"); ax.set_xlim(0.001, 0.2); ax.set_ylim(0, 12)
        ax.set_xlabel(self.t["dur_x"]); ax.set_ylabel(self.t["dur_y"]); ax.set_title(self.t["dur_title"])
        ax.legend(loc="upper right", title=self.t["roche"] + " ···")
        return self._finish(fig, name, "limits: runaway / maximum greenhouse (Kopparapu+ 2013), CHZ after Agol 2011")

    # 6 -------------------------------------------------------------------------------
    def map3d_age(self, df: pd.DataFrame, name="map3d_age"):
        d = df[df["age_thick"].notna()]
        fig = plt.figure(figsize=(10, 9))
        ax = fig.add_subplot(projection="3d", facecolor=BG)
        sc = ax.scatter(d["x"], d["y"], d["z"], c=d["age_thick"], cmap=AGE_CMAP, vmin=0, vmax=8, s=1.6, linewidths=0, alpha=0.9)
        ax.scatter([0], [0], [0], color="#ffd27a", s=40, marker="*")
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False; axis.pane.set_edgecolor("#1e2638"); axis.label.set_color("#aab4c8")
        ax.set_xlabel(f"x ({self.t['pc']})"); ax.set_ylabel(f"y ({self.t['pc']})"); ax.set_zlabel(f"z ({self.t['pc']})")
        ax.set_title(self.t["map_title"]); ax.view_init(elev=22, azim=35)
        cb = fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.02); cb.set_label(self.t["age"]); cb.outline.set_visible(False)
        return self._finish(fig, name)
