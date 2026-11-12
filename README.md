# C28 hyperdiamond: a dynamical-stability screening

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24154048.svg)](https://doi.org/10.5281/zenodo.24154048)

C28 hyperdiamond is a hypothetical carbon crystal proposed in 1992: diamond's lattice with a 28-atom fullerene cage on every site. It has been studied with density-functional theory (DFT) but never synthesised.

This repository asks one question: **is the structure dynamically stable, meaning free of vibrational modes that would tear it apart?** It answers at screening level with a machine-learned interatomic potential, and states plainly what that does and does not show.

## Why I did this

This project started with a physics talk at my school about carbon structures: C60, nanotubes and the wider fullerene family. C28 was the one that caught my interest. It is the smallest fullerene seen in quantity, and its four open bonds point in the same directions as the bonds of a single carbon atom in diamond.

The speaker also talked about potential for producing a hyperdiamond, a crystal built from cages in place of atoms. A material that light, fully bonded, and with a band gap closer to a semiconductor's than diamond's sounded like something that could matter a great deal if it could ever be made. I wanted to look into it properly, not just read about it.

When I did, I found that the idea is over thirty years old and that the main calculations date from the 1990s. I could not find a recent study asking a basic question: would the structure hold together, or would its own vibrations pull it apart? That question can be tested on a laptop, so I tested it.

I am an A-level student, not a researcher. This is a screening result with clear limits, which are set out above. My aim is to add one small, checkable piece of evidence and to make everything open, so that anyone with better tools can confirm it, correct it, or take it further.

## Result

With the MACE-MP-0 potential at three model sizes, C28 hyperdiamond keeps its Fd-3m symmetry after full relaxation and shows **no imaginary phonon modes**.

| Run | Cubic cell edge (Å) | Energy above diamond (eV/atom) | Softest optical mode at Γ (THz) | Next mode (THz) | Imaginary modes |
|---|---|---|---|---|---|
| small, 2×2×2 supercell | 16.028 | 0.617 | 1.48 | 4.71 | none |
| medium, 2×2×2 supercell | 15.945 | 0.580 | 0.87 | 4.19 | none |
| large, 2×2×2 supercell | 15.991 | 0.651 | 1.32 | 4.97 | none |
| medium, 3×3×3 supercell | 15.945 | 0.580 | 0.87 | 4.19 | none |
| Published DFT-LDA (Zeger et al. 1995) | 15.777 | 0.744 | n/a | n/a | n/a |

All three models find a triply degenerate low-frequency mode well separated from the rest of the spectrum. Its frequency is the least certain number here, and it is where an instability would appear first.

![Phonon bands of C28 hyperdiamond, MACE-MP-0 medium](hyperdiamond-stability/hyper_phonon_bands.png)

## What this is not

- **It is not DFT.** MACE-MP-0 is a general-purpose potential trained on Materials Project data. It was never trained on C28 cages.
- **The models are imperfect on carbon.** For diamond they give:

  | Model | Cell edge (measured 3.567 Å) | Highest phonon (measured 1332 cm⁻¹) |
  |---|---|---|
  | small | 3.582 | 997 |
  | medium | 3.570 | 1172 |
  | large | 3.532 | 1331 |

  No single model gets both right. They also disagree on hyperdiamond's cage-to-cage bond (1.45 to 1.50 Å, against 1.539 Å in published DFT).
- **Three sizes of one model family are not three independent methods.** Their agreement is weaker evidence than it looks.
- **Harmonic and zero temperature only.** No anharmonic or finite-temperature effects.
- **Prior work not fully reviewed.** I have not read Seifert, Enyashin and Heine (2005) in full, so I cannot rule out that it addressed stability with an earlier method.

The claim this supports: *a general-purpose machine-learned potential, at three sizes, finds no dynamical instability in C28 hyperdiamond.* It is a lead, not a confirmed result.

## Status and next step

A DFT calculation (Quantum ESPRESSO, PBE) of the relaxed 56-atom cell and its zone-centre phonons is in progress, to test the soft mode directly. Results will be added to this repository.

## Method

1. **Starting structure.** `hyperdiamond.cif` was built with exact Td cage symmetry from the relaxed bond lengths in Zeger et al. (1995): 1.407, 1.479 and 1.510 Å inside the cage, 1.539 Å between cages, a = 15.777 Å. It was checked for 12 pentagons, 4 hexagons, 42 bonds per cage, and 4-coordinate apex atoms.
2. **Relaxation.** Cell and atoms relaxed together (ASE BFGS with a Frechet cell filter) to forces below 0.005 eV/Å.
3. **Phonons.** Finite displacements of 0.01 Å with phonopy: 10 symmetry-distinct displacements of a 448-atom (2×2×2) or 1512-atom (3×3×3) supercell of the 56-atom primitive cell.
4. **Control.** The same procedure on diamond, to measure each model's error on a known material.

## Reproduce

Requires Python 3.12. Runs on a CPU; the 2×2×2 calculation takes about a minute.

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install -r requirements.txt
cd hyperdiamond-stability
python stability.py diamond
python stability.py hyper
python stability.py hyper --model large --tag hyper_large
python stability.py hyper --supercell 3 --tag hyper_medium_sc3
```

## Files

| Path | Contents |
|---|---|
| `hyperdiamond.cif` | Starting structure, 224-atom conventional cell |
| `c28-cage.xyz` | A single C28 cage |
| `c28-hyperdiamond.html` | Interactive 3D explainer page (a page fragment; it loads three.js from a CDN) |
| `hyperdiamond-stability/stability.py` | The relaxation and phonon script |
| `hyperdiamond-stability/*_result.json` | Numbers for each run |
| `hyperdiamond-stability/*_relaxed.cif` | Relaxed structures |
| `hyperdiamond-stability/*_phonopy.yaml` | Displacements and forces, enough to rebuild the force constants |
| `hyperdiamond-stability/*_phonon_bands.png` | Phonon band plots |

Runs without a suffix used the medium model and a 2×2×2 supercell.

## References

- Guo et al., "Uranium stabilization of C28: a tetravalent fullerene", Science 257, 1661 (1992).
- Bylander and Kleinman, Phys. Rev. B 47, 10967 (1993).
- Kaxiras, Zeger, Antonelli and Juan, Phys. Rev. B 49, 8446 (1994). https://doi.org/10.1103/PhysRevB.49.8446
- Zeger, Juan, Kaxiras and Antonelli, Phys. Rev. B 52, 2125 (1995). https://arxiv.org/abs/cond-mat/9502051
- Seifert, Enyashin and Heine, Phys. Rev. B 72, 012102 (2005). https://doi.org/10.1103/PhysRevB.72.012102
- Batatia et al., "A foundation model for atomistic materials chemistry" (MACE-MP-0). https://arxiv.org/abs/2401.00096
- Togo, "First-principles phonon calculations with phonopy and phono3py", J. Phys. Soc. Jpn. 92, 012001 (2023).

## How this was made

The calculations, scripts and write-up were produced with the help of an AI coding assistant (Claude). The structure was derived from the published bond lengths and verified numerically; the literature values quoted come from papers read directly.

## Citing and licence

If you use this work, please cite it using `CITATION.cff` (GitHub shows a "Cite this repository" button).

Archived on Zenodo: https://doi.org/10.5281/zenodo.24154048

- Code: MIT licence, see `LICENSE`.
- Text, data and figures: Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/. Reuse is welcome with credit.
