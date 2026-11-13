"""Follow-up analyses for C28 hyperdiamond, all at MACE-MP-0 (screening) level.

    python extras.py softmode            # what the lowest optical mode at Gamma physically is
    python extras.py elastic             # C11, C12, C44, bulk modulus, Born stability (diamond control included)
    python extras.py md --temp 300       # Langevin molecular dynamics; does the bond network survive?
"""
import argparse, json, time
import numpy as np
from ase import Atoms, units
from ase.build import bulk
from ase.io import read
from ase.optimize import BFGS
from ase.neighborlist import neighbor_list

EV_A3_TO_GPA = 160.21766208


def get_calc(model):
    from mace.calculators import mace_mp
    return mace_mp(model=model, default_dtype="float64", device="cpu")


def bond_graph(atoms, cutoff=1.75):
    i, j = neighbor_list("ij", atoms, cutoff)
    nb = [set() for _ in atoms]
    for a, b in zip(i, j):
        nb[a].add(int(b))
    return nb


def cages_of(atoms):
    """Split atoms into cages: remove bonds joining two 4-coordinate (apex) atoms, then take connected components."""
    nb = bond_graph(atoms)
    apex = {k for k, s in enumerate(nb) if len(s) == 4}
    seen, cages = set(), []
    for s in range(len(atoms)):
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            a = stack.pop()
            comp.append(a)
            for b in nb[a]:
                if a in apex and b in apex:
                    continue
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
        cages.append(sorted(comp))
    return cages, sorted(apex)


def softmode():
    import phonopy
    ph = phonopy.load("hyper_phonopy.yaml", log_level=0)
    ph.run_qpoints([[0, 0, 0]], with_eigenvectors=True)
    d = ph.get_qpoints_dict()
    freqs, vecs = d["frequencies"][0], d["eigenvectors"][0].real
    prim = ph.primitive
    atoms = Atoms(symbols=prim.symbols, cell=prim.cell, scaled_positions=prim.scaled_positions, pbc=True)
    cages, apex = cages_of(atoms)
    out = {"cages": [len(c) for c in cages], "n_apex": len(apex), "modes": []}
    for m in range(3, 9):
        u = vecs[:, m].reshape(-1, 3)
        rec = {"freq_THz": round(float(freqs[m]), 3), "cages": []}
        tot = float((u ** 2).sum())
        trans = []
        for c in cages:
            ref = atoms.positions[c[0]]
            r = np.array([atoms.get_distance(c[0], k, mic=True, vector=True) for k in c]) + ref
            r = r - r.mean(0)
            uc = u[c]
            t = uc.mean(0)
            w = uc - t
            inertia = sum(np.eye(3) * (ri @ ri) - np.outer(ri, ri) for ri in r)
            omega = np.linalg.solve(inertia, sum(np.cross(ri, wi) for ri, wi in zip(r, w)))
            rot = np.cross(omega, r)
            internal = w - rot
            trans.append(t)
            rec["cages"].append({"translation_frac": round(float(len(c) * (t @ t) / tot), 3),
                                 "rotation_frac": round(float((rot ** 2).sum() / tot), 3),
                                 "internal_frac": round(float((internal ** 2).sum() / tot), 3)})
        rec["cage_translations_antiparallel"] = bool(len(trans) == 2 and trans[0] @ trans[1] < 0)
        rec["apex_share_of_motion"] = round(float((u[apex] ** 2).sum() / tot), 3)
        out["modes"].append(rec)
    json.dump(out, open("softmode_result.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def elastic_of(atoms, calc, delta=0.005):
    atoms = atoms.copy()
    atoms.calc = calc
    cell0, spos = atoms.cell[:].copy(), atoms.get_scaled_positions()

    def stress_for(eps):
        a = Atoms(numbers=atoms.numbers, cell=cell0 @ (np.eye(3) + eps), scaled_positions=spos, pbc=True)
        a.calc = calc
        BFGS(a, logfile=None).run(fmax=0.003, steps=300)  # relax internal coordinates at fixed strain
        return a.get_stress(voigt=True)

    e1 = np.zeros((3, 3)); e1[0, 0] = delta
    e4 = np.zeros((3, 3)); e4[1, 2] = e4[2, 1] = delta / 2  # engineering shear gamma = delta
    s1p, s1m = stress_for(e1), stress_for(-e1)
    s4p, s4m = stress_for(e4), stress_for(-e4)
    c11 = (s1p[0] - s1m[0]) / (2 * delta) * EV_A3_TO_GPA
    c12 = ((s1p[1] - s1m[1]) + (s1p[2] - s1m[2])) / 2 / (2 * delta) * EV_A3_TO_GPA
    c44 = (s4p[3] - s4m[3]) / (2 * delta) * EV_A3_TO_GPA
    return {"C11_GPa": round(float(c11), 1), "C12_GPa": round(float(c12), 1), "C44_GPa": round(float(c44), 1),
            "bulk_modulus_GPa": round(float((c11 + 2 * c12) / 3), 1),
            "born_stable": bool(c11 - c12 > 0 and c11 + 2 * c12 > 0 and c44 > 0)}


def elastic(model):
    calc = get_calc(model)
    out = {"model": f"MACE-MP-0 {model}"}
    for name, f in (("diamond", "diamond_relaxed.cif"), ("hyperdiamond", "hyper_relaxed.cif")):
        t = time.time()
        out[name] = elastic_of(read(f), calc)
        out[name]["seconds"] = round(time.time() - t, 1)
        print(name, out[name], flush=True)
    out["diamond_experiment"] = {"C11_GPa": 1079, "C12_GPa": 124, "C44_GPa": 578, "bulk_modulus_GPa": 443}
    json.dump(out, open("elastic_result.json", "w"), indent=1)


def md(model, temp, steps, tag):
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    atoms = read("hyper_relaxed.cif")
    atoms.calc = get_calc(model)
    nb0 = bond_graph(atoms)
    rng = np.random.default_rng(28)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp, rng=rng)
    Stationary(atoms)
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=temp, friction=0.01 / units.fs, rng=rng)
    log, t0 = [], time.time()
    mean_pos, n_avg = np.zeros((len(atoms), 3)), 0
    ref = atoms.positions.copy()
    for s in range(0, steps, 50):
        dyn.run(50)
        nb = bond_graph(atoms, cutoff=1.9)
        broken = sum(len(a - b) for a, b in zip(nb0, nb)) // 2
        new = sum(len(b - a) for a, b in zip(nb0, nb)) // 2
        d = atoms.get_all_distances(mic=True)
        maxbond = max(d[a, b] for a in range(len(atoms)) for b in nb0[a])
        log.append({"step": s + 50, "T_K": round(float(atoms.get_temperature()), 1),
                    "E_pot_eV_per_atom": round(float(atoms.get_potential_energy() / len(atoms)), 4),
                    "bonds_broken": int(broken), "bonds_new": int(new), "longest_original_bond_A": round(float(maxbond), 3)})
        if s + 50 > steps // 2:  # average the second half
            disp = atoms.positions - ref
            disp -= np.round(disp @ np.linalg.inv(atoms.cell[:])) @ atoms.cell[:]
            mean_pos += disp; n_avg += 1
        if (s // 50) % 10 == 9:
            print(log[-1], f"{time.time() - t0:.0f}s", flush=True)
    mean_disp = np.linalg.norm(mean_pos / max(n_avg, 1), axis=1)
    out = {"model": f"MACE-MP-0 {model}", "temperature_K": temp, "steps": steps, "timestep_fs": 1.0, "atoms": len(atoms),
           "max_bonds_broken_at_any_sample": max(x["bonds_broken"] for x in log),
           "max_new_bonds_at_any_sample": max(x["bonds_new"] for x in log),
           "longest_original_bond_seen_A": max(x["longest_original_bond_A"] for x in log),
           "mean_T_second_half_K": round(float(np.mean([x["T_K"] for x in log[len(log) // 2:]])), 1),
           "max_time_averaged_displacement_A": round(float(mean_disp.max()), 3),
           "network_intact": bool(all(x["bonds_broken"] == 0 and x["bonds_new"] == 0 for x in log)), "trace": log}
    json.dump(out, open(f"{tag}_result.json", "w"), indent=1)
    print({k: v for k, v in out.items() if k != "trace"})


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["softmode", "elastic", "md"])
    ap.add_argument("--model", default="medium")
    ap.add_argument("--temp", type=float, default=300)
    ap.add_argument("--steps", type=int, default=3000)
    a = ap.parse_args()
    if a.what == "softmode":
        softmode()
    elif a.what == "elastic":
        elastic(a.model)
    else:
        md(a.model, a.temp, a.steps, f"md_{int(a.temp)}K")
