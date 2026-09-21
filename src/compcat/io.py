"""Reading the archive.

The archive stores one energy column per configuration and no reference state.
Two mutually inconsistent records of the slab composition exist; both are
exposed here so that any analysis can be repeated against either. See
`scripts/08_composition_record_test.py` for the test that distinguishes them.
"""
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

ARCHIVE = Path(__file__).resolve().parents[2] / "data" / "raw"

SYSTEMS = {
    "AgAuCuPdPt": ("Ag Au Cu Pd Pt".split(), "fcc"),
    "AlCeCuZnZr": ("Al Ce Cu Zn Zr".split(), "bcc"),
    "AlCoCuFeNi": ("Al Co Cu Fe Ni".split(), "bcc"),
    "AlCoMoNiZn": ("Al Co Mo Ni Zn".split(), "bcc"),
    "CaCeCoLaNi": ("Ca Ce Co La Ni".split(), "fcc"),
    "CoCrNiVZn":  ("Co Cr Ni V Zn".split(),  "bcc"),
    "CoCuFeMnNi": ("Co Cu Fe Mn Ni".split(), "fcc"),
    "CoCuNiSnZn": ("Co Cu Ni Sn Zn".split(), "fcc"),
    "CoFeMgMnNi": ("Co Fe Mg Mn Ni".split(), "bcc"),
    "CoNiPdRhRu": ("Co Ni Pd Rh Ru".split(), "fcc"),
    "GaNbTaTiZr": ("Ga Nb Ta Ti Zr".split(), "bcc"),
    "HfNbTaTiZr": ("Hf Nb Ta Ti Zr".split(), "bcc"),
}

# species the archived composition parser invents; they are not in any slab
SPURIOUS = {"C", "F", "H", "O", "N"}

_DIR = re.compile(r"outputs/([A-Za-z0-9]+)-(FCC|BCC)/")
_PAIR = re.compile(r"([A-Z][a-z]?)(\d+)")


def _open(name):
    p = ARCHIVE / name
    if p.suffix == ".gz" or (ARCHIVE / (name + ".gz")).exists():
        p = p if p.suffix == ".gz" else ARCHIVE / (name + ".gz")
        return gzip.open(p, "rt", newline="")
    return open(p, newline="")


def system_of(traj_file):
    """Alloy system for a trajectory path.

    Rows belonging to CoNiPdRhRu carry a path with no system prefix, a
    consequence of that per-system file having been lost and the rows recovered
    from the combined table.
    """
    head = traj_file.split("/")[0]
    if head.endswith(("+fcc", "+bcc")):
        return head.split("+")[0]
    return "CoNiPdRhRu"


def label_composition(directory_label):
    """Composition read from the trajectory directory name, e.g.
    Co20Ni17Pd32Rh20Ru19. Subscripts sum to 108."""
    pairs = _PAIR.findall(directory_label)
    total = sum(int(n) for _, n in pairs)
    return {e: int(n) / total for e, n in pairs}


def parsed_composition(composition_json):
    """Metal composition recovered from the archived composition field by
    discarding the species the parser invents and renormalising."""
    d = json.loads(composition_json)
    metals = {k: v for k, v in d.items() if k not in SPURIOUS}
    total = sum(metals.values())
    return {k: v / total for k, v in metals.items()}


def load_energy_table(name="all_energy_results.csv.gz", lattice=None, systems=None):
    """Yield the archived rows, optionally restricted by lattice or system.

    Each row is augmented with `system`, `composition_label` and `lattice`.
    """
    with _open(name) as fh:
        for row in csv.DictReader(fh):
            m = _DIR.search(row["traj_file"])
            if not m:
                continue
            lab, lat = m.group(1), m.group(2).lower()
            if lattice and lat != lattice:
                continue
            sysname = system_of(row["traj_file"])
            if systems and sysname not in systems:
                continue
            try:
                row["energy"] = float(row["energy"])
            except (TypeError, ValueError):
                continue
            row["system"] = sysname
            row["composition_label"] = lab
            row["lattice"] = lat
            yield row


def load_configurations(name="all_energy_results.csv.gz", lattice="fcc",
                        systems=None, window=None):
    """Group the archive by configuration.

    Returns {system: {(composition_label, position): {adsorbate: E_arch}}}.
    A configuration is one slab at one surface position; the adsorbates on it
    share a reference state, which is what makes within-configuration
    differences meaningful.
    """
    from .window import in_window
    out = defaultdict(lambda: defaultdict(dict))
    for row in load_energy_table(name, lattice=lattice, systems=systems):
        if window is not None and not in_window(row["energy"], window):
            continue
        key = (row["composition_label"], row["position"])
        out[row["system"]][key][row["adsorbate"]] = row["energy"]
    return {k: dict(v) for k, v in out.items()}


def load_parsed_compositions(system="CoNiPdRhRu",
                             name="df_ads_alignn_results.csv.gz"):
    """The second composition record, keyed by directory label."""
    out = {}
    with _open(name) as fh:
        for row in csv.DictReader(fh):
            m = _DIR.search(row["traj_file"])
            if not m or system_of(row["traj_file"]) != system:
                continue
            out[m.group(1)] = parsed_composition(row["composition"])
    return out
