"""Screening test of C28 hyperdiamond's dynamical stability.

Method: MACE-MP-0 machine-learned interatomic potential (trained on PBE DFT data),
full relaxation of cell and atoms, then finite-displacement phonons with phonopy.
This is a screening-level calculation, not DFT. A result here is a reason to run DFT, not a replacement for it.

Usage:
    python stability.py diamond     # sanity check against known diamond values
    python stability.py hyper       # relax hyperdiamond and compute its phonons
Optional: --model small|medium|large   --supercell N   --fmax 0.005
"""
import argparse, json, sys, time
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
import spglib
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

THZ_TO_CM = 33.35641


def get_calc(model):
    from mace.calculators import mace_mp
    return mace_mp(model=model, default_dtype="float64", device="cpu")


def spacegroup(atoms, symprec=1e-3):
    return spglib.get_spacegroup((atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers), symprec=symprec)


def primitive(atoms, symprec=1e-3):
    cell, pos, num = spglib.standardize_cell((atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers),
                                             to_primitive=True, no_idealize=False, symprec=symprec)
    return Atoms(numbers=num, cell=cell, scaled_positions=pos, pbc=True)


def relax(atoms, calc, fmax, log):
    atoms.calc = calc
    opt = BFGS(FrechetCellFilter(atoms), logfile=log)
    opt.run(fmax=fmax, steps=600)
    return atoms


def bond_lengths(atoms, cutoff=1.75):
    d = atoms.get_all_distances(mic=True)
    vals = np.round(d[(d > 0.5) & (d < cutoff)], 3)
    uniq, counts = np.unique(vals, return_counts=True)
    return {float(u): int(c // 2) for u, c in zip(uniq, counts)}


def phonons(atoms, calc, n, tag):
    cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:], scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(cell, supercell_matrix=np.diag(n), primitive_matrix=np.eye(3))
    ph.generate_displacements(distance=0.01)
    scs = ph.supercells_with_displacements
    print(f"  {len(scs)} displaced supercells of {len(scs[0])} atoms", flush=True)
    forces = []
    for i, sc in enumerate(scs):
        t = time.time()
        a = Atoms(symbols=sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
        print(f"  displacement {i + 1}/{len(scs)}  {time.time() - t:.1f} s", flush=True)
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_mesh([8, 8, 8] if len(atoms) > 10 else [20, 20, 20])
    mesh = ph.get_mesh_dict()["frequencies"]
    ph.run_qpoints([[0, 0, 0]])
    gamma = ph.get_qpoints_dict()["frequencies"][0]
    ph.auto_band_structure(npoints=51, plot=False)
    band = ph.get_band_structure_dict()
    bmin = min(float(np.min(f)) for f in band["frequencies"])
    try:
        import matplotlib
        matplotlib.use("Agg")
        ph.auto_band_structure(npoints=51, plot=True).savefig(f"{tag}_phonon_bands.png", dpi=140)
    except Exception as e:  # plotting is optional
        print("  (band plot skipped:", e, ")")
    ph.save(f"{tag}_phonopy.yaml", settings={"force_constants": False})
    return {"n_displacements": len(scs), "supercell_atoms": len(scs[0]),
            "min_freq_mesh_THz": float(mesh.min()), "min_freq_bandpath_THz": bmin,
            "max_freq_THz": float(mesh.max()), "max_freq_cm-1": float(mesh.max() * THZ_TO_CM),
            "gamma_THz": [round(float(x), 3) for x in gamma]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["diamond", "hyper", "custom"])
    ap.add_argument("--model", default="medium")
    ap.add_argument("--supercell", default=None, help="N or N,N,N")
    ap.add_argument("--cif", default=None, help="structure file for the custom mode")
    ap.add_argument("--fmax", type=float, default=0.005)
    ap.add_argument("--tag", default=None, help="label for output files, so runs do not overwrite each other")
    args = ap.parse_args()
    calc = get_calc(args.model)
    sc = [int(x) for x in args.supercell.split(",")] if args.supercell else None
    if sc and len(sc) == 1:
        sc = sc * 3
    tag = args.tag or args.what
    out = {"model": f"MACE-MP-0 {args.model}", "what": args.what}

    if args.what == "diamond":
        at = bulk("C", "diamond", a=3.567)
        n = sc or [4, 4, 4]
    else:
        at = primitive(read(args.cif if args.what == "custom" else "hyperdiamond.cif"))
        n = sc or [2, 2, 2]
        out["start_bonds_A"] = bond_lengths(at)
    print(f"{args.what}: {len(at)} atoms in primitive cell, start space group {spacegroup(at)}", flush=True)

    t = time.time()
    relax(at, calc, args.fmax, f"{tag}_relax.log")
    a_conv = float((4 * at.get_volume()) ** (1 / 3))  # conventional cubic edge of an fcc lattice
    out.update({"relax_seconds": round(time.time() - t, 1), "space_group_after": spacegroup(at),
                "a_conventional_A": round(a_conv, 4), "cell_lengths_angles": [round(float(x), 4) for x in at.cell.cellpar()],
                "volume_A3_per_atom": round(float(at.get_volume() / len(at)), 4), "energy_eV_per_atom": float(at.get_potential_energy() / len(at)),
                "max_force_eV_A": float(np.abs(at.get_forces()).max()), "bonds_A": bond_lengths(at),
                "density_g_cm3": round(len(at) * 12.011 * 1.66053907 / at.get_volume(), 4)})
    write(f"{tag}_relaxed.cif", at)
    print(json.dumps(out, indent=1), flush=True)

    print("phonons...", flush=True)
    t = time.time()
    out["phonons"] = phonons(at, calc, n, tag)
    out["phonon_seconds"] = round(time.time() - t, 1)
    lo = min(out["phonons"]["min_freq_mesh_THz"], out["phonons"]["min_freq_bandpath_THz"])
    out["verdict"] = ("no imaginary modes found (dynamically stable at this level of theory)" if lo > -0.3
                      else f"imaginary mode found, lowest {lo:.2f} THz (dynamically unstable at this level of theory)")
    json.dump(out, open(f"{tag}_result.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    sys.exit(main())
