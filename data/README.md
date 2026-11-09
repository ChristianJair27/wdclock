# Input data

`examples/gaia_100pc.py` expects the Gaia EDR3 white dwarf catalogue of Gentile Fusillo et al. (2021),
VizieR J/MNRAS/508/3877, cut at parallax > 10 mas and P_WD > 0.75 (16,676 rows), with the VizieR column
names: `WDJname, GaiaEDR3, RA_ICRS, DE_ICRS, Plx, Pwd, Gmag, BP_RP, TeffH, loggH, MassH` plus the derived
`GMAG` (absolute G from the parallax), `dist_pc` and Cartesian `x, y, z` in pc. The catalogue is not stored in
this repository; fetch it from VizieR (CDS or the Harvard mirror) and add the derived columns, or point the
script at any table with Teff and mass columns (`--teff-col`, `--mass-col` in the CLI).

The Montreal cooling sequences are downloaded automatically to `~/.wdclock/montreal`.
