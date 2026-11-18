"""Define numeric identifiers for ZSoil analysis-driver types.

This module is part of ``zsoil-output``, a Python SDK developed for use with ZSoil FEM to
read model data and output-result files.

The :class:`DriversType` constants identify initial-state, stability, time-dependent,
arc-length, dynamic, pushover, and eigenmode analysis records.
"""



class DriversType ():
    INITIAL_STATE = 1
    STABILITY = 2
    TIME_DEPENDENT = 3
    ARC_LENGTH = 4
    DYNAMICS = 5
    PUSHOVER = 6
    EIGENMODES = 7
