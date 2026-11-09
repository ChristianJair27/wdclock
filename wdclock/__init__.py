"""wdclock: white-dwarf cooling ages and habitable-zone windows from Montreal evolutionary sequences.

    from wdclock import CoolingGrid, HabitableZone
    g = CoolingGrid("thick")
    g.age(25369, 1.018) / 1e6          # Sirius B: ~126 Myr (Bond et al. 2017)
    HabitableZone(g).now(6000, 0.6)    # inner/outer edge in AU and the periods in days
"""
from .montreal import CoolingGrid, MASSES, download
from .hz import HabitableZone, S_LIMITS, period_days, roche_limit_au, transit_probability, transit_duration_min, transit_depth
from .mestel import cooling_age as mestel_age

__version__ = "0.1.0"
__all__ = ["CoolingGrid", "HabitableZone", "MASSES", "S_LIMITS", "download", "mestel_age", "period_days",
           "roche_limit_au", "transit_probability", "transit_duration_min", "transit_depth"]
