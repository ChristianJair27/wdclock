# What we learned building wdclock (2026-11-09)

Notes by a self-taught team: a hobbyist (Christian) and an AI agent. Written so that the next person who
wants to date a white dwarf from a catalogue does not repeat our detours.

## The physics, in the order it made sense to us

1. **Why a white dwarf is a clock.** Nothing burns. The star is a degenerate carbon-oxygen core whose electrons
   supply the pressure but almost no heat capacity; the heat is in the ions (3/2 kT each). It leaks out through a
   thin, non-degenerate envelope. Known mass plus measured temperature therefore gives elapsed time, provided
   you know (a) how fast the envelope leaks, (b) how much heat the ions hold, (c) what else releases heat.
2. **Mestel's law is (a)+(b) with Kramers opacity and ideal-gas ions.** It gives age ∝ L^(-5/7). Zingale's
   notebook writes it as t = 6 Myr (A/12)^-1 (M/Msun)^(5/7) (μ/2)^(-2/7) [L^(-5/7) − L0^(-5/7)]. Against the
   Montreal 0.6 Msun sequence it is within a factor 1.5–2 from 0.1 to 1e-4 Lsun, too fast at the hot end (it
   ignores neutrino cooling, which dominates above ~25,000 K) and wrong at the cool end.
3. **Crystallisation is (c).** Below about 1e-4 Lsun (Teff ~ 5000 K at 0.6 Msun) the core freezes, releasing
   latent heat and the gravitational energy of oxygen sinking in the C/O mixture. The full sequence shows the
   age jumping from ~3 Gyr to ~8 Gyr over a small range in L: the star "stalls". One third of the 100 pc sample
   has a core more than half crystallised. Gaia saw this stall as a pile-up in the colour-magnitude diagram
   (Tremblay et al. 2019, Nature 565, 202); our isochrones bunch up exactly there.
4. **The hydrogen layer is the leak valve.** A thick layer (10^-4 of the mass) insulates better than a thin
   one (10^-10): for most stars the age difference is a few per cent (median 0.09 Gyr here), but it is a
   systematic you cannot remove star by star without a spectrum. It also moves the stall: in the thin 0.6 Msun
   sequence the core is 19 % crystallised at 6000 K and the clock runs at 2.8 Gyr per 1000 K there, against 0 %
   and 1.2 Gyr per 1000 K for the thick one, which only reaches that pace near 5200 K. The most populated
   temperature range of the sample, 5000–6300 K, therefore lands at 5–6 Gyr in the thin-layer histogram and
   makes the bump visible in `age_histogram.png`; it is a model feature, not a burst of star formation.
5. **Gaia gives Teff and log g from three photometric bands plus a parallax.** Gentile Fusillo et al. (2021)
   fitted G, BP, RP with pure-H, pure-He and mixed atmospheres and used these same Montreal sequences for the
   mass–radius relation, so feeding their Teff and mass back into the sequences is self-consistent. The
   failure modes are instructive: Sirius B comes out at 10,900 K and 0.44 Msun because Sirius A contaminates the
   photometry (the real values are 25,400 K and 1.02 Msun). Always look at the star before trusting its age.
6. **The habitable zone of a cooling star marches inward.** With L ~ 1e-4 Lsun the zone sits at ~0.01 AU, with
   orbital periods of 4–32 hours (Agol 2011), inside the orbit of Mercury by a factor of 40 and outside the
   Roche limit (0.0054 AU for Earth density) only for Teff above ~4500 K. A planet at a fixed distance enters the
   zone when the star has dimmed enough and leaves it when the star is too dim: we reproduce Agol's continuously
   habitable zone (≥ 3 Gyr) at 0.0078–0.0174 AU for 0.6 Msun, with the longest stay, 6.2 Gyr, at 0.0117 AU
   (Agol found up to 8 Gyr with older cooling tracks; the crystallisation stall is what makes the stay long).
7. **Why it matters for a transit search.** For the TESS-bright stars (G < 17) the habitable zone today has
   periods of 100–230 hours (median), because the bright stars are the young, hot ones whose zones are still far
   out. For the cool, old stars the habitable period is 10–40 hours. A BLS that stops at a few days covers the
   habitable zone of 81 % of the G < 17 sample. If every one of those 3,339 stars had an Earth in its habitable
   zone, geometry alone (median transit probability 0.56 %) would let us see 7.5 of them transit, each for about
   two minutes, with depths near 100 %.

## Methods we copied from professionals, and why

- Interpolate in log(age) along a single-mass track as a function of log(Teff), then linearly in mass between
  the two bracketing tracks. Teff is monotonic along a cooling track after the first 1e5 yr; age spans six
  orders of magnitude, so log–log is the only interpolation that behaves.
- Validate against stars with dynamical masses and independent ages (Sirius B, 40 Eri B). The literature uses
  exactly these two because HST astrometry gives their masses to 1–3 %.
- Draw isochrones from the mass sequences, then convert to Gaia colours with the Montreal synthetic photometry
  table on its regular (Teff, log g) grid. Our first attempt interpolated the photometry table along constant
  log g and produced zig-zags at 6–10 Gyr, where very old, cool white dwarfs turn blue again (collision-induced
  absorption by H2); constant log g is not constant mass and the table is sparse there.
- State what is a measurement, what is a model, and what is a definition. The age is a model reading of a
  measurement; the habitable zone is a definition on top of a model.

## From Michael Zingale's notes

The structure of his stellar-physics course (zingale.github.io/stars) is: equation of state first, then
energy transport, then the simplest analytic solution (Mestel), then the full model and the difference between
them. We followed the same order and the comparison figure (`mestel_vs_montreal.png`) is the result. His "WD
cooling" notebook compares Mestel with a MESA track; we compare with the Montreal sequences, which are what the
Gaia catalogue uses, and we add crystallisation shading so the reader sees where and why the analytic law fails.

## What we would do next

- Add the pure-He and mixed-atmosphere solutions (Table_DB, `TeffHe`, `MassHe`) with a per-star choice from
  spectral type where known (DA vs DB/DC/DZ lists).
- Propagate the Gaia parameter uncertainties (the catalogue gives `e_TeffH`, `e_loggH`) by Monte Carlo.
- Add progenitor lifetimes with an initial-final mass relation (Cummings et al. 2018) to turn cooling ages into
  total ages, and compare the result with the published local star-formation histories properly, with the
  Gaia selection function.
- Feed `hz_p_in_d`–`hz_p_out_d` into the transit search as the period window to weight.
