"""compcat - reproduction code for the HEA/DMR machine-learning screening study.

The package is deliberately small. Four modules cover everything the analysis
needs: reading the archive, applying the physical window, forming the
offset-invariant descriptors, and the handful of statistics used to put error
bars on the results.
"""
from .io import (ARCHIVE, SYSTEMS, load_energy_table, load_configurations,
                 label_composition, parsed_composition, system_of)
from .window import WINDOW, in_window, retention
from .descriptors import DESCRIPTORS, SIGNS, configuration_descriptors, composite_score
from .stats import pearson, bootstrap_ci, bootstrap_null, split_half

__version__ = "1.0.0"
__all__ = [
    "ARCHIVE", "SYSTEMS", "load_energy_table", "load_configurations",
    "label_composition", "parsed_composition", "system_of",
    "WINDOW", "in_window", "retention",
    "DESCRIPTORS", "SIGNS", "configuration_descriptors", "composite_score",
    "pearson", "bootstrap_ci", "bootstrap_null", "split_half",
]
