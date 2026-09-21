"""Reference benchmark on clean Ni(111).

Establishes, on a surface whose energetics are known independently, what the
archived quantity can and cannot be. Ni(111) 3x3x4, 12 A vacuum, bottom two
layers fixed, relaxed to fmax = 0.05 eV/A with MACE-MP-0.

On the CHx series the result tracks published DFT with r = 0.95 and a mean
absolute deviation of 0.77 eV, which is the warrant for the diagnostic. The
agreement does not extend to the oxygen-bearing species, so inferences drawn
from this benchmark are restricted to the CHx series.

Produces: data/derived/ni111_reference.json
Reproduces: Table 5, Figure 7 and Section 2.9.
Runtime: tens of minutes on one core. Requires ase and mace-torch.
"""
import json
import numpy as np
from ase import Atoms
from ase.build import fcc111, add_adsorbate, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from mace.calculators import mace_mp

CALC = mace_mp(model="medium", default_dtype="float64", device="cpu")

FMAX, STEPS = 0.05, 200
VAC = 12.0


def relax(atoms, label, fix_bottom=None):
    atoms.calc = CALC
    if fix_bottom is not None:
        atoms.set_constraint(FixAtoms(mask=[a.position[2] < fix_bottom for a in atoms]))
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=FMAX, steps=STEPS)
    e = atoms.get_potential_energy()
    print(f"  {label:22s} E = {e:12.4f} eV   ({len(atoms)} atoms)", flush=True)
    return atoms, e


def gas(name, spec):
    a = spec() if callable(spec) else spec.copy()
    a.set_cell([15, 15, 15]); a.center(); a.pbc = False
    return relax(a, f"gas {name}")[1]


# ---------------------------------------------------------------- gas phase
GAS_SPECS = {
    "CH4": lambda: molecule("CH4"),
    "CH3": lambda: molecule("CH3"),
    "CH2": lambda: molecule("CH2_s1A1d"),
    "CH":  lambda: Atoms("CH", positions=[[0, 0, 0], [0, 0, 1.12]]),
    "C":   lambda: Atoms("C",  positions=[[0, 0, 0]]),
    "H":   lambda: Atoms("H",  positions=[[0, 0, 0]]),
    "O":   lambda: Atoms("O",  positions=[[0, 0, 0]]),
    "CO":  lambda: molecule("CO"),
    "CO2": lambda: molecule("CO2"),
    "H2":  lambda: molecule("H2"),
    "H2O": lambda: molecule("H2O"),
    "OH":  lambda: molecule("OH"),
}

print("gas-phase references (MACE-MP-0):", flush=True)
EGAS = {k: gas(k, v) for k, v in GAS_SPECS.items()}

# ---------------------------------------------------------------- clean slab
print("\nclean Ni(111) 3x3x4, 12 A vacuum, bottom 2 layers fixed:", flush=True)
clean = fcc111("Ni", size=(3, 3, 4), vacuum=VAC / 2)
zs = np.unique(np.round(clean.positions[:, 2], 3))
fixz = zs[1] + 0.1                       # fix the bottom two layers
clean, E_SLAB = relax(clean, "clean Ni(111)", fix_bottom=fixz)

# ---------------------------------------------------------------- adsorption
SITES = {"C": "fcc", "H": "fcc", "O": "fcc", "CH": "fcc", "CH2": "fcc",
         "CH3": "fcc", "CO": "fcc", "CH4": "ontop", "CO2": "ontop",
         "H2O": "ontop", "OH": "fcc", "H2": "ontop"}
HEIGHT = {"CH4": 3.2, "CO2": 3.2, "H2": 2.6, "H2O": 2.2}

print("\nadsorption on Ni(111):", flush=True)
res = {}
for name in GAS_SPECS:
    slab = clean.copy()
    ads = GAS_SPECS[name]()
    h = HEIGHT.get(name, 1.6)
    add_adsorbate(slab, ads, height=h, position=SITES[name])
    slab.set_constraint(FixAtoms(mask=[a.position[2] < fixz for a in slab]))
    try:
        _, e = relax(slab, f"{name}/Ni(111)", fix_bottom=fixz)
        res[name] = {"E_slab_ads": e, "Eads": e - E_SLAB - EGAS[name]}
    except Exception as exc:
        print(f"  !! {name}: {exc}", flush=True)

out = {"E_slab_clean": E_SLAB, "E_gas": EGAS, "adsorption": res,
       "model": "MACE-MP-0 medium", "slab": "Ni(111) 3x3x4, 12 A vacuum"}
json.dump(out, open("../data/derived/ni111_reference.json", "w"), indent=2)

print("\n" + "=" * 62)
print(f"{'species':8s} {'Eads (eV)':>12s}   {'literature DFT (eV)':>22s}")
LIT = {"C": -6.9, "H": -2.8, "O": -5.6, "CH": -6.4, "CH2": -4.0, "CH3": -1.9,
       "CO": -1.9, "CH4": -0.15, "OH": -3.3, "H2O": -0.35, "CO2": -0.1, "H2": -0.1}
for k, v in res.items():
    print(f"{k:8s} {v['Eads']:12.2f}   {LIT.get(k, float('nan')):22.2f}")
