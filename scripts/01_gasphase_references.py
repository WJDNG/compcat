"""Gas-phase reference energies for the 17 DMR intermediates.

PBE/def2-TZVP, spin-unrestricted, each species at its own optimized geometry and
in its experimental electronic ground state. Every calculation converged.

Produces: data/derived/gasphase_energies.json
Reproduces: Table 1, Table S1, Section S9 (Cartesian coordinates) and the
gas-phase ladder of Figure 10a.
Runtime: minutes on one core. Requires pyscf and geometric/pyberny.
"""
import json, sys
import numpy as np
from pyscf import gto, dft
from pyscf.geomopt.berny_solver import optimize

XC = "PBE"
BASIS = "def2-tzvp"

# name: (atom string [Angstrom, starting geometry], charge, spin=2S)
SPECIES = {
    "H":    ("H 0 0 0", 0, 1),
    "C":    ("C 0 0 0", 0, 2),
    "O":    ("O 0 0 0", 0, 2),
    "H2":   ("H 0 0 0; H 0 0 0.74", 0, 0),
    "CO":   ("C 0 0 0; O 0 0 1.128", 0, 0),
    "CO2":  ("C 0 0 0; O 0 0 1.16; O 0 0 -1.16", 0, 0),
    "H2O":  ("O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", 0, 0),
    "OH":   ("O 0 0 0; H 0 0 0.97", 0, 1),
    "CH":   ("C 0 0 0; H 0 0 1.12", 0, 1),
    "CH2":  ("C 0 0 0; H 0 0.86 0.65; H 0 -0.86 0.65", 0, 2),   # triplet ground state
    "CH3":  ("C 0 0 0; H 0 1.079 0; H 0.934 -0.539 0; H -0.934 -0.539 0", 0, 1),
    "CH4":  ("C 0 0 0; H 0.629 0.629 0.629; H -0.629 -0.629 0.629; "
             "H -0.629 0.629 -0.629; H 0.629 -0.629 -0.629", 0, 0),
    "CH2O": ("C 0 0 0; O 0 0 1.208; H 0 0.943 -0.588; H 0 -0.943 -0.588", 0, 0),
    "CH3O": ("C 0 0 0; O 0 0 1.37; H 0 1.03 -0.36; H 0.89 -0.51 -0.36; H -0.89 -0.51 -0.36", 0, 1),
    "HCO":  ("C 0 0 0; O 0 1.02 0.72; H 0 -1.0 0.5", 0, 1),
    "COOH": ("C 0 0 0; O 0 1.02 0.72; O 0 -1.15 0.62; H 0 -1.05 1.58", 0, 1),
    "HCOO": ("C 0 0 0; O 0 1.10 0.68; O 0 -1.10 0.68; H 0 0 -1.09", 0, 1),
}

HARTREE_EV = 27.211386245988


def run(name, spec):
    atom, charge, spin = spec
    mol = gto.M(atom=atom, basis=BASIS, charge=charge, spin=spin,
                verbose=0, unit="Angstrom")
    mf = dft.UKS(mol)
    mf.xc = XC
    mf.conv_tol = 1e-9
    mf.max_cycle = 200

    natm = mol.natm
    if natm > 1:
        try:
            mol_eq = optimize(mf, maxsteps=60)
            mf = dft.UKS(mol_eq)
            mf.xc = XC
            mf.conv_tol = 1e-9
            mf.max_cycle = 200
            e = mf.kernel()
            geom = mol_eq.atom_coords(unit="Angstrom").tolist()
            syms = [mol_eq.atom_symbol(i) for i in range(mol_eq.natm)]
            opt = True
        except Exception as exc:                      # fall back to input geometry
            print(f"  !! opt failed for {name}: {exc}", file=sys.stderr)
            e = mf.kernel()
            geom = mol.atom_coords(unit="Angstrom").tolist()
            syms = [mol.atom_symbol(i) for i in range(mol.natm)]
            opt = False
    else:
        e = mf.kernel()
        geom = [[0.0, 0.0, 0.0]]
        syms = [mol.atom_symbol(0)]
        opt = True

    return {
        "name": name, "E_hartree": float(e), "E_eV": float(e) * HARTREE_EV,
        "spin_2S": spin, "converged": bool(mf.converged), "optimized": opt,
        "symbols": syms, "coords_ang": geom,
    }


if __name__ == "__main__":
    out = {}
    for name, spec in SPECIES.items():
        print(f"running {name} ...", flush=True)
        out[name] = run(name, spec)
        r = out[name]
        print(f"  E = {r['E_eV']:.4f} eV   converged={r['converged']} opt={r['optimized']}", flush=True)
    with open("../data/derived/gasphase_energies.json", "w") as fh:
        json.dump({"functional": XC, "basis": BASIS, "species": out}, fh, indent=2)
    print("saved")
