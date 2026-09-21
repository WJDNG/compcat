"""Facet control for the bcc branch.

Rebuilds one affected composition on both (111) and (110) with the same builder
and the same size, and recomputes adsorption with MACE-MP-0 against the
closed-shell references of script 02.

The facet choice is a real error - bcc(111) gives four atomic levels at 0.83 A
and a mean surface coordination of 7.0, against 2.03 A and 10.0 for (110), which
lies 0.70 eV per atom lower. It is not the cause of the archived magnitudes:
both facets return adsorption energies inside the physical range.

Produces: data/derived/bcc_facet_control.json
Reproduces: the bcc paragraph of Section 3.7 and the control row of Table S6.
Runtime: hours on one core. Requires ase and mace-torch.
"""
import json, itertools
import numpy as np
from ase import Atoms
from ase.build import bcc110, bcc111, molecule
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.neighborlist import natural_cutoffs, NeighborList
from mace.calculators import mace_mp

RNG = np.random.default_rng(0)
CALC = mace_mp(model="medium", default_dtype="float64", device="cpu")
FMAX, STEPS, VAC = 0.05, 300, 12.0

# one bcc system from the archived screen
ELEMENTS = ["Al", "Co", "Cu", "Fe", "Ni"]          # AlCoCuFeNi (bcc)
A0 = 2.87                                          # bcc lattice constant, A


def decorate(slab, seed=0):
    rng = np.random.default_rng(seed)
    n = len(slab)
    per = n // len(ELEMENTS)
    syms = list(itertools.chain.from_iterable([[e] * per for e in ELEMENTS]))
    syms += [ELEMENTS[0]] * (n - len(syms))
    rng.shuffle(syms)
    slab.set_chemical_symbols(syms)
    return slab


def layer_stats(slab):
    z = np.sort(np.unique(np.round(slab.positions[:, 2], 2)))
    d = np.diff(z)
    return len(z), (float(d.min()) if len(d) else 0.0)


def mean_cn(slab, surface_only=True):
    cut = natural_cutoffs(slab, mult=1.05)
    nl = NeighborList(cut, self_interaction=False, bothways=True)
    nl.update(slab)
    zmax = slab.positions[:, 2].max()
    cns = []
    for i in range(len(slab)):
        if surface_only and slab.positions[i, 2] < zmax - 0.5:
            continue
        cns.append(len(nl.get_neighbors(i)[0]))
    return float(np.mean(cns))


def relax(atoms, label, fixz):
    atoms.calc = CALC
    atoms.set_constraint(FixAtoms(mask=[a.position[2] < fixz for a in atoms]))
    BFGS(atoms, logfile=None).run(fmax=FMAX, steps=STEPS)
    e = atoms.get_potential_energy()
    print(f"  {label:26s} E = {e:12.3f} eV  ({len(atoms)} atoms)", flush=True)
    return atoms, e


def gas(spec):
    a = spec()
    a.set_cell([15, 15, 15]); a.center(); a.pbc = False
    a.calc = CALC
    if len(a) > 1:
        BFGS(a, logfile=None).run(fmax=FMAX, steps=200)
    return a.get_potential_energy()


GAS = {"CH4": lambda: molecule("CH4"), "H2": lambda: molecule("H2"),
       "H2O": lambda: molecule("H2O"), "CO": lambda: molecule("CO"),
       "CH3": lambda: molecule("CH3"),
       "C": lambda: Atoms("C", positions=[[0, 0, 0]]),
       "H": lambda: Atoms("H", positions=[[0, 0, 0]]),
       "O": lambda: Atoms("O", positions=[[0, 0, 0]])}
print("gas phase (MACE-MP-0):", flush=True)
EG = {k: gas(v) for k, v in GAS.items()}
for k, v in EG.items():
    print(f"  {k:4s} {v:10.3f}", flush=True)

# closed-shell chemical potentials, same convention as the Ni(111) benchmark
MU = {"C": EG["CH4"] - 2 * EG["H2"], "H": .5 * EG["H2"],
      "O": EG["H2O"] - EG["H2"], "CH3": EG["CH4"] - .5 * EG["H2"],
      "CO": EG["CO"]}
ADS = ["C", "H", "O", "CH3", "CO"]

out = {}
for facet, builder in (("111", bcc111), ("110", bcc110)):
    print(f"\n=== bcc({facet}) ===", flush=True)
    slab = builder("Fe", size=(3, 3, 4), a=A0, vacuum=VAC / 2)
    slab = decorate(slab, seed=1)
    nlev, dmin = layer_stats(slab)
    zs = np.unique(np.round(slab.positions[:, 2], 3))
    fixz = zs[len(zs) // 2] + 0.05
    print(f"  {len(slab)} atoms, {nlev} atomic levels, min spacing {dmin:.2f} A", flush=True)
    clean, E0 = relax(slab.copy(), f"clean bcc({facet})", fixz)
    cn = mean_cn(clean)
    print(f"  mean surface coordination {cn:.2f}", flush=True)

    res = {}
    for name in ADS:
        s = clean.copy()
        a = GAS[name]()
        # place above the highest atom, roughly at a hollow between top atoms
        top = s.positions[s.positions[:, 2] > s.positions[:, 2].max() - 0.4]
        xy = top[:, :2].mean(axis=0) if len(top) > 2 else top[0, :2]
        a.translate([xy[0] - a.positions[:, 0].mean(),
                     xy[1] - a.positions[:, 1].mean(),
                     s.positions[:, 2].max() + 1.7 - a.positions[:, 2].min()])
        s += a
        try:
            _, e = relax(s, f"{name}/bcc({facet})", fixz)
            res[name] = {"E": e, "Eform": e - E0 - MU[name]}
            print(f"      -> Eform({name}) = {res[name]['Eform']:+.2f} eV", flush=True)
        except Exception as exc:
            print(f"  !! {name}: {exc}", flush=True)
    out[facet] = {"E_clean": E0, "n_atoms": len(clean), "n_levels": nlev,
                  "min_spacing": dmin, "mean_surface_cn": cn, "ads": res}

json.dump(out, open("../data/derived/bcc_facet_control.json", "w"), indent=2)
print("\n" + "=" * 64)
print(f"{'species':8s} {'bcc(111)':>12s} {'bcc(110)':>12s}")
for k in ADS:
    a = out["111"]["ads"].get(k, {}).get("Eform", float("nan"))
    b = out["110"]["ads"].get(k, {}).get("Eform", float("nan"))
    print(f"{k:8s} {a:12.2f} {b:12.2f}")
print("\ngeometry:")
for f in ("111", "110"):
    d = out[f]
    print(f"  bcc({f}): {d['n_atoms']} atoms, {d['n_levels']} levels, "
          f"dmin {d['min_spacing']:.2f} A, surface CN {d['mean_surface_cn']:.2f}")
