# NOTES: MATE70020 Problem Set 2, AI-generated reference solution

> AI-generated (Claude Code, Claude Opus 5.5, 10 November 2026). This is for comparison only and is not a student submission.

## Deliverables

| File | What it is |
|---|---|
| `submission.pdf` / `submission.tex` | The 3-page answer (a)–(d), with AI notice box and footer on every page. The LaTeX follows the user's PSet 1 template: siunitx units, `\maketitle`, a Method section with a footnote link to the runs, and labelled, cross-referenced tables and figure. Kept on top of that template: the AI notice, the footer and the Name/CID placeholder |
| `runs/BN.in`, `runs/*.UPF` | Byte-identical copies of the course files (`cmp` checked). The originals in the repo root are untouched |
| `runs/BN.bands.in` | Main band-structure input (Γ–K–M, 91 k-points) |
| `runs/BN.bands_dense.in`, `runs/BN.bands_nearK.in` | Gap-location checks (271 k-points on the full path; 201 k-points within 0.17 Å⁻¹ of K) |
| `runs/BN.kz.in` | Part (a) check: eigenvalues at k_z = 0 vs k_z = π/c |
| `runs/bands_pp.in` | `bands.x` input |
| `runs/01_scf_bands/` … `04_kz/` | Each job: SCF output, bands output, XML with all eigenvalues, `bands.x` files |
| `runs/99_sanity_ecut120/` | 120 Ry numerical sanity check. **Not** used for any reported number |
| `runs/scripts/run_all.sh` | Reproduces every run, then the analysis and the figure |
| `runs/scripts/analyse.py` | Every number in the PDF, plus all self-checks → `runs/results.txt`, `runs/results.json` |
| `runs/scripts/plot_bands.py` | Figure → `runs/bands.pdf` (used in LaTeX) and `runs/bands.png` |
| `runs/scripts/sanity_ecut.sh` | Runs the 120 Ry check |
| `runs/scripts/qe_env_packages.txt` | Exact conda-forge package versions used |

## Environment and commands

- Quantum ESPRESSO **7.6** (conda-forge build `mpi_openmpi_hdb26816_104`, OpenMPI 5.0.11, OpenBLAS 0.3.34), run with 4 MPI processes.
- Python 3 with numpy 2.5 and matplotlib 3.11. LaTeX is TeX Live (pdflatex).

```bash
# QE via micromamba (micromamba itself fetched as a conda-forge package, because micro.mamba.pm was blocked)
curl -O https://conda.anaconda.org/conda-forge/linux-64/micromamba-2.3.0-0.tar.bz2
tar xjf micromamba-2.3.0-0.tar.bz2 bin/micromamba
MAMBA_ROOT_PREFIX=$PWD/mamba ./bin/micromamba create -y -p $PWD/qe-env -c conda-forge qe
```

```bash
cd runs && QE_BIN=/path/to/qe-env/bin NP=4 bash scripts/run_all.sh   # about 2 min
```

```bash
cd runs && QE_BIN=/path/to/qe-env/bin bash scripts/sanity_ecut.sh    # optional, about 1.5 min
```

```bash
pdflatex submission.tex && pdflatex submission.tex
```

Inside `run_all.sh`, each job is `pw.x -in BN.in` (SCF), then `pw.x -in BN.<x>.in` (`calculation='bands'`), then `bands.x -in bands_pp.in`. Each job runs in its own directory.

## Changes to the supplied input (all of them)

- **SCF:** none. `BN.in` was run exactly as supplied: LDA-PZ, `ecutwfc = 60`, `ibrav = 4`, `A = 2.504`, `C = 20.0`, same positions and pseudopotentials, `6 6 1 0 0 0`, `nbnd = 8`, default `conv_thr = 1e-6` Ry.
- **Bands runs:** only two changes.
  - Added `calculation = 'bands'` to `&CONTROL`.
  - Replaced the `K_POINTS` card with `crystal_b` Γ (0,0,0) 60 → K (1/3,1/3,0) 30 → M (1/2,0,0) 1. The weights are chosen so the spacing is uniform (0.028 Å⁻¹), because |ΓK| = 2|KM|.
- `nbnd = 8` was kept: 4 occupied and 4 empty bands, which is enough to see the CBM.
- **Checks only:**
  - `K_POINTS crystal` with 4 points (k_z test).
  - Denser `crystal_b` paths.
  - `ecutwfc = 120` in `99_sanity_ecut120` only.

## Problems hit and how they were fixed

1. **Ubuntu `apt` QE 6.7 crashes on start-up.** Every `pw.x` run aborted with `*** buffer overflow detected ***` in `__snprintf_chk`. The binary is stripped, so this is a packaging/fortify bug and not an input problem. **Fix:** use QE 7.6 from conda-forge instead, via micromamba as above.
2. **Network restrictions.**
   - `www.quantum-espresso.org`, `micro.mamba.pm`, `api.anaconda.org`, nature.com and Crossref were not reachable from the container.
   - **INPUT_PW:** read from the QE source file `PW/Doc/INPUT_PW.def` on GitLab (develop branch). It gives `v1 = a(1,0,0), v2 = a(-1/2,sqrt(3)/2,0), v3 = a(0,0,c/a)` for ibrav=4, and says `A`, `B`, `C` are in Å. This was cross-checked against the crystal and reciprocal axes that `pw.x` prints. The PDF cites the official web page, which is the page the problem sheet links to.
   - **Nature Communications reference:** checked against catalogue records found by web search (PMC, RePEc/IDEAS, the Nottingham repository): title, the 12 authors, volume 10, article 2639, year 2019, DOI 10.1038/s41467-019-10610-5. The abstract states a 6.1 eV direct gap measured by deep-UV photoluminescence and reflectance. The full text was not read.
3. **matplotlib was not installed.** Fixed with `pip install matplotlib`.
4. **`pw.x` stops printing eigenvalues when there are more than 100 k-points.** The dense paths therefore have no eigenvalues in the `.out` file. **Fix:** keep `pwscf.save/data-file-schema.xml` from each bands run as `BN.<x>.xml` and parse that. The `bands.x` `.gnu` data agree with the XML to 5×10⁻⁵ eV.
5. **The fine path near K first showed the VBM one step (0.0017 Å⁻¹) past K, 0.24 meV higher than at K.**
   - Band 4 jitters at the 0.1–0.3 meV level, and the jumps follow changes in the plane-wave count (5758→5762→5766). Band 5 is smooth.
   - This is basis-set discretisation at the fixed 60 Ry cutoff (sharp |k+G|² < E_cut sphere).
   - At 120 Ry the band-4 maximum is exactly at K and symmetric around it.
   - Likewise, the occupied-band "k_z dispersion" of ≤ 8.6 meV at 60 Ry drops to ≤ 0.013 meV at 120 Ry.
   - **Conclusion:** the VBM is at K. The PDF reports the 60 Ry results and says "same location".

## Assumptions

- `A` and `C` are lengths in Å. `C = 20.0` is the cell height c, not the ratio c/a; QE prints a(3) = 7.987 alat, which confirms this.
- K = (1/3,1/3,0) and M = (1/2,0,0) in crystal coordinates.
  - In Cartesian coordinates, K is at 60° and M at 30° from x.
  - M is the midpoint of the side from K to K′ = (2/3,−1/3,0), which matches the problem-sheet sketch (K corner, M on the adjacent side, K′ the next corner).
  - Any symmetry-equivalent choice gives the same bands.
- The calculation is spin-unpolarised, as in BN.in, so there are 2 electrons per band.
- Valence charges were taken from the UPF headers: B 3, N 5. That gives 8 electrons and 4 occupied bands.
- The 6.1 eV experimental value is taken as given by the question. The experiment measures an optical (PL/reflectance) gap. The PDF does not discuss excitonic or quasiparticle effects, because they are beyond Lectures 2–4 and the markers advise against going beyond what is asked.
- Model name in the AI notice: Claude Opus 5.5. The session metadata reports the configured and last-served model as `claude-opus-5-5`; the marketing name is used rather than the raw ID string.
- Date: 10 November 2026, the session date and the PSet 2 due date in the Lecture 3 notes.

## Self-verification (all done before writing; full output in `runs/results.txt`)

| Check | Result |
|---|---|
| a_i · b_j = 2π δ_ij | max deviation 8.9×10⁻¹⁶ Å·Å⁻¹ → PASS |
| b_i vs `pw.x` "reciprocal axes" (2π/a units) | max deviation 4.6×10⁻⁷ (print precision) → PASS |
| K, M Cartesian vs hexagon geometry | \|ΓK\| = 4π/3a = 1.6728 Å⁻¹ (corner = side length); \|ΓM\| = 2π/(√3 a) = 1.4487 Å⁻¹ = (√3/2)\|ΓK\| = \|b₁\|/2; \|KM\| = \|ΓK\|/2. All to 10⁻¹⁶ → PASS |
| Wigner–Seitz test | K equidistant from 3 nearest G (corner); M from 2 (side midpoint) |
| SCF convergence | all 4 SCF runs (and the 120 Ry one) print "convergence has been achieved in 7 iterations"; E = −25.77488163 Ry; estimated accuracy 3.2×10⁻⁷ Ry < 10⁻⁶ |
| Bands runs | no "eigenvalues not converged" warnings (`ethr` 1.25×10⁻⁸) |
| Gap location, dense path | 271 k: VBM K, CBM K, gap 4.5504 eV (same as 91 k) |
| Gap location, fine path near K | 201 k: CBM at K; band-4 maximum within 0.24 meV of its K value (basis noise, see above); gap at K 4.5504 eV |
| SCF grid (contains K) | HOMO −4.3395, LUMO 0.2051 → 4.545 eV. The 6 meV difference from the bands run is due to the default `conv_thr` |
| Part (a) k_z test | occupied bands shift ≤ 8.6 meV between k_z = 0 and π/c (≤ 0.013 meV at 120 Ry), vs 3.6 eV in-plane dispersion of band 4 |
| Full reproducibility | `run_all.sh` rerun from scratch → `results.txt` byte-identical |
| PSet re-read | (a) explained; (b) vectors from the docs plus reciprocal vectors, symbolic and numeric; (c) K and M (and Γ) as combinations of b₁, b₂; (d) plot Γ→K→M, direct/indirect, value, and "how does it compare" each answered explicitly |

Other numbers that are not in the PDF:
- Direct gap at Γ: 6.19 eV (this one is in the PDF). At M: 5.58 eV.
- Occupied bandwidth: 17.56 eV.
- Gap at 120 Ry: 4.578 eV, i.e. +0.03 eV from the cutoff. The 60 Ry value of 4.55 eV is therefore converged to about 0.03 eV.

## Things deliberately left out of the PDF (possible marker questions)

- Empty "nearly-free-electron" vacuum states at Γ (bands 6–8) do disperse along k_z, by up to 0.39 eV. This is why part (a) says the **occupied** KS energies are k_z-independent; only the occupied states build n(r).
- The 36 grid points reduce to 7 irreducible k-points by symmetry. This is not asked.
- Bulk hBN is indirect (5.95 eV, per the same paper's abstract). This is not asked.
- Excitons and GW: see Assumptions.

## Self-assessment against the marks

| Part | Marks | Assessment |
|---|---|---|
| (a) | 1 | Full: links the 20 Å vacuum/periodic images, k_z-independence of occupied states, the density sum over k (Lecture 4), and the thin BZ (\|b₃\| = 0.314 vs 2.897 Å⁻¹), with a numerical check. Risk: low. |
| (b) | 2 | Full: doc vectors with values and units, cross products shown, symbolic and numeric b_i in Å⁻¹ (and 2π/a units), check a_i·b_j = 2πδ_ij. Risk: low. |
| (c) | 3 | Full: derivation from the BZ definition (perpendicular bisectors; K equidistant from Γ, b₁, b₂ → x = y = 1/3; M = b₁/2), crystal and Cartesian coordinates, geometric check, consistent with the sketch. Risk: low. A marker may expect a different but equivalent K (e.g. (2/3,1/3,0) is **not** equivalent; (−1/3,2/3,0) is). The answer gives the one adjacent to M as in the sketch. |
| (d) | 4 | Full: plot with labels, VBM = 0 and occupied/empty distinguished; electron count → 4 occupied bands; direct gap at K, 4.55 eV; compared with 6.1 eV (−1.55 eV, −25%); explanation from Lecture 3 (LDA/GGA gap underestimation, hybrids/PBE0 open the gap) plus Lecture 2 (KS eigenvalues not physical) and a one-line self-interaction justification. Risk: low to medium only if a marker wants the exciton/optical-gap caveat, which was deliberately omitted (beyond the course). |

Expected: 10/10 for a generous marker; worst plausible 9/10.
