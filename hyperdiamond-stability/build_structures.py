"""Build starting structures for two relatives of C28 hyperdiamond.

    c28c.cif              C28 cage + one carbon atom in the zincblende arrangement (Zeger et al. 1995, a = 9.685 A)
    hyperlonsdaleite.cif  C28 cages on the lonsdaleite (hexagonal diamond) lattice (Seifert et al. 2005 proposed this phase)

The cage geometry is the one used for hyperdiamond. Both structures are starting points for relaxation,
not published coordinates.
"""
import itertools
import numpy as np
from ase import Atoms
from ase.io import read, write

cage = read("../c28-cage.xyz").positions          # 28 atoms, centred on the origin; atoms 0-3 are the apexes
V = cage[:4] / np.linalg.norm(cage[:4], axis=1)[:, None]
R_APEX = float(np.linalg.norm(cage[0]))
D_CC = 2 * R_APEX + 1.539                          # cage-centre separation in hyperdiamond


def report(name, atoms):
    d = atoms.get_all_distances(mic=True)
    np.fill_diagonal(d, 9)
    coord = (d < 1.75).sum(1)
    vals, counts = np.unique(np.round(d[d < 1.75], 3), return_counts=True)
    print(f"{name}: {len(atoms)} atoms, coordination counts {dict(zip(*np.unique(coord, return_counts=True)))}, "
          f"bond lengths {dict(zip(vals.tolist(), (counts // 2).tolist()))}, shortest non-bonded {d[d >= 1.75].min():.3f} A")


# --- C28C: fcc cell, cage at the origin, single carbon at (1/4,1/4,1/4)
a = 9.685
cell = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]) * a
pos = np.vstack([cage, [[a / 4] * 3]])
c28c = Atoms("C29", positions=pos, cell=cell, pbc=True)
report("C28C", c28c)
write("c28c.cif", c28c)

# --- hyperlonsdaleite: ideal wurtzite-type site lattice with cage-centre separation D_CC
ah = D_CC * np.sqrt(8 / 3)
ch = ah * np.sqrt(8 / 3)
hcell = np.array([[ah, 0, 0], [-ah / 2, ah * np.sqrt(3) / 2, 0], [0, 0, ch]])
u = 3 / 8
sites = np.array([[1 / 3, 2 / 3, 0], [2 / 3, 1 / 3, .5], [1 / 3, 2 / 3, u], [2 / 3, 1 / 3, .5 + u]]) @ hcell
images = [np.array(s) @ hcell for s in itertools.product((-1, 0, 1), repeat=3)]
positions = []
for s in sites:
    nb = [t + im - s for t in sites for im in images if abs(np.linalg.norm(t + im - s) - D_CC) < 0.05]
    assert len(nb) == 4, len(nb)
    N = np.array(nb) / D_CC
    best = None
    for perm in itertools.permutations(range(4)):       # the cage is achiral, so any orthogonal map is allowed
        U, _, Wt = np.linalg.svd(N[list(perm)].T @ V)
        Q = U @ Wt
        err = np.abs(V @ Q.T - N[list(perm)]).max()
        if best is None or err < best[0]:
            best = (err, Q)
    assert best[0] < 1e-6, best[0]
    positions.append(cage @ best[1].T + s)
hl = Atoms("C112", positions=np.vstack(positions), cell=hcell, pbc=True)
hl.wrap()
report("hyperlonsdaleite", hl)
write("hyperlonsdaleite.cif", hl)
print(f"hex cell a = {ah:.3f} A, c = {ch:.3f} A")
