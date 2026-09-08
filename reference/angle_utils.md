# Angle helpers for ggNetView layout parameters

Several layout parameters in ggNetView accept an angle. Internally
everything is stored in radians, but for convenience the public API
allows the user to pass either radians (`pi/2`) or degrees (`45`)
directly. The unit is auto-detected from the magnitude.
