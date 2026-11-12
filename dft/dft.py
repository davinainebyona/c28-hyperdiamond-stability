"""DFT check of C28 hyperdiamond with Quantum ESPRESSO (PBE, PAW), run on GitHub Actions.

Steps (each a separate job, because one runner is limited to 6 hours):
    python dft.py relax NAME          # variable-cell relaxation        NAME = diamond | hyper
    python dft.py disp NAME INDEX     # forces for one displaced cell
    python dft.py collect NAME        # zone-centre phonons from all displaced cells

Phonons here use the primitive cell itself as the supercell, so they are exact only at the
zone centre (Gamma). That is where the soft cage-rotation mode sits.
"""
import json, os, subprocess, sys
import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
import spglib
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
PSEUDO = {"C": "C.pbe-n-kjpaw_psl.1.0.0.UPF"}
ECUTWFC, ECUTRHO = 45.0, 360.0          # Ry; SSSP-recommended for this pseudopotential
KPTS = {"diamond": (8, 8, 8), "hyper": (2, 2, 2)}
DISPLACEMENT = 0.03                      # Angstrom; larger than the ML run so soft-mode forces sit above DFT noise
NP = os.environ.get("NP", "4")
THZ_TO_CM = 33.35641


def symmetrised_primitive(atoms):
    cell, pos, num = spglib.standardize_cell((atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers),
                                             to_primitive=True, no_idealize=False, symprec=1e-3)
    return Atoms(numbers=num, cell=cell, scaled_positions=pos, pbc=True)


def start_structure(name):
    if name == "diamond":
        return bulk("C", "diamond", a=3.57)
    return symmetrised_primitive(read(os.path.join(HERE, "..", "hyperdiamond-stability", "hyper_relaxed.cif")))


def run_pw(atoms, label, calculation, kpts, conv_thr):
    inp = {"control": {"calculation": calculation, "pseudo_dir": os.path.join(HERE, "pseudo"), "outdir": os.path.join(HERE, f"tmp_{label}"),
                       "prefix": label, "tprnfor": True, "tstress": True, "disk_io": "low", "forc_conv_thr": 1.0e-4, "etot_conv_thr": 1.0e-6},
           "system": {"ecutwfc": ECUTWFC, "ecutrho": ECUTRHO, "occupations": "fixed"},
           "electrons": {"conv_thr": conv_thr, "mixing_beta": 0.5, "electron_maxstep": 200},
           "ions": {}, "cell": {"press_conv_thr": 0.3}}
    fin, fout = os.path.join(HERE, f"{label}.in"), os.path.join(HERE, f"{label}.out")
    write(fin, atoms, format="espresso-in", input_data=inp, pseudopotentials=PSEUDO, kpts=kpts)
    # ASE also writes empty &FCP and &RISM blocks, which older pw.x versions do not recognise
    import re
    drop = "FCP|RISM|IONS|CELL" if calculation == "scf" else "FCP|RISM"
    nl = chr(10)
    txt = re.sub("&(" + drop + ")[^/]*/[^" + nl + "]*" + nl, "", open(fin).read())
    with open(fin, "w", newline=nl) as f:
        f.write(txt)
    with open(fout, "w") as f:
        subprocess.run(["mpirun", "--oversubscribe", "-np", NP, "pw.x", "-in", fin], stdout=f, stderr=subprocess.STDOUT, cwd=HERE, check=False)
    text = open(fout, errors="replace").read()
    if "JOB DONE" not in text:
        sys.exit(f"pw.x did not finish for {label}; last lines:\n" + "\n".join(text.splitlines()[-40:]))
    return read(fout, format="espresso-out", index=-1), text


def save_structure(name, atoms, extra):
    os.makedirs(RES, exist_ok=True)
    d = {"cell": atoms.cell[:].tolist(), "scaled_positions": atoms.get_scaled_positions().tolist(), "numbers": atoms.numbers.tolist()}
    d.update(extra)
    json.dump(d, open(os.path.join(RES, f"{name}_relaxed.json"), "w"), indent=1)
    write(os.path.join(RES, f"{name}_dft_relaxed.cif"), atoms)


def load_structure(name):
    d = json.load(open(os.path.join(RES, f"{name}_relaxed.json")))
    return Atoms(numbers=d["numbers"], cell=d["cell"], scaled_positions=d["scaled_positions"], pbc=True), d


def bonds(atoms, cutoff=1.75):
    d = atoms.get_all_distances(mic=True)
    v, c = np.unique(np.round(d[(d > 0.5) & (d < cutoff)], 3), return_counts=True)
    return {float(a): int(b // 2) for a, b in zip(v, c)}


def relax(name):
    a0 = start_structure(name)
    at, text = run_pw(a0, f"{name}_relax", "vc-relax", KPTS[name], 1.0e-9)
    at = Atoms(numbers=at.numbers, cell=at.cell[:], scaled_positions=at.get_scaled_positions(), pbc=True)
    sg = spglib.get_spacegroup((at.cell[:], at.get_scaled_positions(), at.numbers), symprec=1e-3)
    at_sym = symmetrised_primitive(at)
    energy = None
    for line in text.splitlines():
        if line.strip().startswith("!"):
            energy = float(line.split("=")[1].split()[0]) * 13.605693123  # Ry -> eV, last occurrence is the final structure
    extra = {"space_group": sg, "n_atoms": len(at_sym), "a_conventional_A": round(float((4 * at_sym.get_volume()) ** (1 / 3)), 4),
             "energy_eV_per_atom": energy / len(at), "density_g_cm3": round(len(at_sym) * 12.011 * 1.66053907 / at_sym.get_volume(), 4),
             "bonds_A": bonds(at_sym), "n_ionic_steps": text.count("number of scf cycles") or None,
             "settings": {"functional": "PBE", "pseudopotential": PSEUDO["C"], "ecutwfc_Ry": ECUTWFC, "ecutrho_Ry": ECUTRHO, "kpoints": list(KPTS[name])}}
    save_structure(name, at_sym, extra)
    print(json.dumps(extra, indent=1))


def phonopy_for(atoms):
    cell = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell[:], scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(cell, supercell_matrix=np.eye(3, dtype=int), primitive_matrix=np.eye(3))
    ph.generate_displacements(distance=DISPLACEMENT)
    return ph


def disp(name, index):
    atoms, _ = load_structure(name)
    ph = phonopy_for(atoms)
    scs = ph.supercells_with_displacements
    json.dump({"n_displacements": len(scs)}, open(os.path.join(RES, f"{name}_ndisp.json"), "w"))
    if index >= len(scs):
        print(f"index {index} not needed: only {len(scs)} displacements")
        return
    sc = scs[index]
    a = Atoms(symbols=sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
    out, _ = run_pw(a, f"{name}_disp{index}", "scf", KPTS[name], 1.0e-11)
    np.save(os.path.join(RES, f"{name}_forces_{index}.npy"), out.get_forces())
    print("max force eV/A:", float(np.abs(out.get_forces()).max()))


def collect(name):
    atoms, meta = load_structure(name)
    ph = phonopy_for(atoms)
    n = len(ph.supercells_with_displacements)
    missing = [i for i in range(n) if not os.path.exists(os.path.join(RES, f"{name}_forces_{i}.npy"))]
    if missing:
        sys.exit(f"missing force files for displacements {missing} of {n}")
    forces = np.array([np.load(os.path.join(RES, f"{name}_forces_{i}.npy")) for i in range(n)])
    forces -= forces.mean(axis=1, keepdims=True)  # remove net force (numerical drift)
    ph.forces = forces
    ph.produce_force_constants()
    ph.symmetrize_force_constants()
    ph.run_qpoints([[0, 0, 0]])
    g = ph.get_qpoints_dict()["frequencies"][0]
    optical = g[3:]
    out = {"what": name, "method": "Quantum ESPRESSO, PBE, PAW; finite displacements in the primitive cell (Gamma only)",
           "displacement_A": DISPLACEMENT, "n_displacements": n, "relaxed": {k: v for k, v in meta.items() if k not in ("cell", "scaled_positions", "numbers")},
           "acoustic_THz": [round(float(x), 3) for x in g[:3]], "lowest_optical_THz": [round(float(x), 3) for x in optical[:9]],
           "highest_THz": round(float(g.max()), 3), "highest_cm-1": round(float(g.max() * THZ_TO_CM), 1),
           "gamma_THz": [round(float(x), 3) for x in g]}
    out["verdict"] = ("no imaginary optical mode at Gamma" if optical.min() > -0.3 else f"imaginary mode at Gamma, lowest {optical.min():.2f} THz")
    json.dump(out, open(os.path.join(RES, f"{name}_dft_result.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "gamma_THz"}, indent=1))


if __name__ == "__main__":
    cmd, name = sys.argv[1], sys.argv[2]
    if cmd == "relax":
        relax(name)
    elif cmd == "disp":
        disp(name, int(sys.argv[3]))
    elif cmd == "collect":
        collect(name)
    else:
        sys.exit("unknown command")
