"""C28X solids with X = Si, Ge, Sn: the zincblende arrangement of C28C with the single carbon replaced.

Zeger et al. (1995) suggested these "would be equally stable, and with similar electronic properties" but did
not calculate them. This script builds the starting structures and the element reference energies
(MACE-MP-0 medium). The phonon test itself is run with stability.py in custom mode.

    python group4.py build      # writes c28si.cif, c28ge.cif, c28sn.cif
    python group4.py refs       # relaxes diamond-structure C, Si, Ge, Sn and writes group4_refs.json
    python group4.py summary    # formation energies, densities and phonon verdicts
"""
import json, sys
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.data import atomic_masses, atomic_numbers

X = {"Si": 1.89, "Ge": 1.95, "Sn": 2.14}          # starting C-X bond lengths (typical single bonds), relaxed afterwards
A0 = {"C": 3.567, "Si": 5.431, "Ge": 5.658, "Sn": 6.489}


def build():
    cage = read("../c28-cage.xyz").positions
    r_apex = float(np.linalg.norm(cage[0]))
    for el, b in X.items():
        a = 4 * (r_apex + b) / np.sqrt(3)
        cell = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * a
        at = Atoms("C28" + el, positions=np.vstack([cage, [[a / 4] * 3]]), cell=cell, pbc=True)
        d = at.get_all_distances(mic=True)
        np.fill_diagonal(d, 9)
        print(f"C28{el}: a = {a:.3f} A, C-{el} bonds: {np.sum(np.abs(d[28] - b) < 0.01)}, shortest other contact {d[d > 0.5].min():.3f} A")
        write(f"c28{el.lower()}.cif", at)


def refs(model="medium"):
    from mace.calculators import mace_mp
    from ase.optimize import BFGS
    from ase.filters import FrechetCellFilter
    calc = mace_mp(model=model, default_dtype="float64", device="cpu")
    out = {}
    for el, a in A0.items():
        at = bulk(el, "diamond", a=a)
        at.calc = calc
        BFGS(FrechetCellFilter(at), logfile=None).run(fmax=0.002, steps=300)
        out[el] = {"energy_eV_per_atom": float(at.get_potential_energy() / 2), "a_A": round(float((4 * at.get_volume()) ** (1 / 3)), 4), "a_measured_A": a}
        print(el, out[el], flush=True)
    json.dump(out, open("group4_refs.json", "w"), indent=1)


def summary():
    ref = json.load(open("group4_refs.json"))
    rows = []
    for el in ["C", "Si", "Ge", "Sn"]:
        tag = "c28" + el.lower()
        try:
            r = json.load(open(f"{tag}_result.json"))
        except FileNotFoundError:
            continue
        at = read(f"{tag}_relaxed.cif")
        n = len(at)
        e_form = (r["energy_eV_per_atom"] * n - 28 * ref["C"]["energy_eV_per_atom"] - ref[el]["energy_eV_per_atom"]) / n
        mass = sum(atomic_masses[atomic_numbers[s]] for s in at.get_chemical_symbols())
        d = at.get_all_distances(mic=True)
        xi = n - 1
        cx = np.sort(d[xi][d[xi] > 0.5])[:4]
        p = r["phonons"]
        rows.append({"solid": "C28" + el, "space_group": r["space_group_after"], "a_A": r["a_conventional_A"],
                     "formation_energy_eV_per_atom_vs_elements": round(float(e_form), 3),
                     "density_g_cm3": round(float(mass * 1.66053907 / at.get_volume()), 3),
                     "cage_to_X_bond_A": [round(float(x), 3) for x in cx],
                     "min_freq_mesh_THz": round(p["min_freq_mesh_THz"], 3), "min_freq_bandpath_THz": round(p["min_freq_bandpath_THz"], 3),
                     "lowest_optical_gamma_THz": p["gamma_THz"][3:6], "verdict": r["verdict"], "bonds_A": r["bonds_A"]})
        print(json.dumps(rows[-1]))
    json.dump(rows, open("group4_summary.json", "w"), indent=1)


if __name__ == "__main__":
    {"build": build, "refs": refs, "summary": summary}[sys.argv[1]]()
