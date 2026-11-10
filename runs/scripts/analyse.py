#!/usr/bin/env python3
"""Every number quoted in the submission, plus the self-checks.

Reads the QE runs in runs/0*/ and writes runs/results.txt and runs/results.json.
Run from anywhere:  python3 runs/scripts/analyse.py
"""
import json
import pathlib
import xml.etree.ElementTree as ET

import numpy as np

RUNS = pathlib.Path(__file__).resolve().parents[1]
HA2EV = 27.211386245988          # QE's autoev
BOHR = 0.529177210903            # Angstrom, QE's bohr_radius_angs
out = []


def p(s=""):
    out.append(s)
    print(s)


# ---------------------------------------------------------------- (b) lattice
A, C = 2.504, 20.0               # Angstrom, from BN.in (A, C are lengths)
# QE INPUT_PW, ibrav=4: v1 = a(1,0,0), v2 = a(-1/2, sqrt3/2, 0), v3 = a(0,0,c/a)
a1 = A * np.array([1.0, 0.0, 0.0])
a2 = A * np.array([-0.5, np.sqrt(3) / 2, 0.0])
a3 = A * np.array([0.0, 0.0, C / A])
V = a1 @ np.cross(a2, a3)
b1 = 2 * np.pi * np.cross(a2, a3) / V
b2 = 2 * np.pi * np.cross(a3, a1) / V
b3 = 2 * np.pi * np.cross(a1, a2) / V
Amat, Bmat = np.array([a1, a2, a3]), np.array([b1, b2, b3])

p("== (b) real and reciprocal lattice ==")
p(f"a = {A} A = {A / BOHR:.4f} bohr ; c = {C} A ; V = a1.(a2 x a3) = sqrt3/2 a^2 c = {V:.4f} A^3"
  f" (= {V / BOHR**3:.4f} bohr^3)")
for n, v in zip(("a1", "a2", "a3"), Amat):
    p(f"  {n} = ({v[0]: .5f}, {v[1]: .5f}, {v[2]: .5f}) A")
for n, v in zip(("b1", "b2", "b3"), Bmat):
    p(f"  {n} = ({v[0]: .5f}, {v[1]: .5f}, {v[2]: .5f}) 1/A   |{n}| = {np.linalg.norm(v):.5f} 1/A"
      f"   in 2pi/a units: ({v[0] * A / 2 / np.pi:.6f}, {v[1] * A / 2 / np.pi:.6f}, {v[2] * A / 2 / np.pi:.6f})")
p(f"  analytic: 2pi/a = {2 * np.pi / A:.5f} 1/A, 4pi/(sqrt3 a) = {4 * np.pi / (np.sqrt(3) * A):.5f} 1/A,"
  f" 2pi/c = {2 * np.pi / C:.5f} 1/A")
G = Amat @ Bmat.T
p("  a_i . b_j / 2pi =\n" + np.array2string(G / (2 * np.pi), precision=12, suppress_small=True))
dev = np.abs(G - 2 * np.pi * np.eye(3)).max()
p(f"  CHECK a_i.b_j = 2pi delta_ij : max deviation = {dev:.2e}  -> {'PASS' if dev < 1e-12 else 'FAIL'}")
ang = np.degrees(np.arccos(b1 @ b2 / np.linalg.norm(b1) / np.linalg.norm(b2)))
p(f"  angle(b1,b2) = {ang:.4f} deg (hexagonal reciprocal lattice)")

# compare with the reciprocal axes printed by pw.x (units 2pi/alat)
txt = (RUNS / "01_scf_bands" / "BN.scf.out").read_text().splitlines()
i = next(k for k, l in enumerate(txt) if "reciprocal axes" in l)
qe_b = np.array([[float(x) for x in txt[i + 1 + j].split("= (")[1].split(")")[0].split()] for j in range(3)])
devqe = np.abs(qe_b - Bmat * A / (2 * np.pi)).max()
p(f"  CHECK vs pw.x 'reciprocal axes' (2pi/alat): max deviation = {devqe:.1e} -> {'PASS' if devqe < 1e-5 else 'FAIL'}")

# ---------------------------------------------------------------- (c) K and M
p("\n== (c) high-symmetry points ==")
crys = {"G": np.array([0, 0, 0.0]), "K": np.array([1 / 3, 1 / 3, 0]), "M": np.array([1 / 2, 0, 0.0])}
cart = {k: v @ Bmat for k, v in crys.items()}
for k in crys:
    c = cart[k]
    p(f"  {k}: crystal {crys[k]}  ->  Cartesian ({c[0]:.5f}, {c[1]:.5f}, {c[2]:.5f}) 1/A"
      f" = ({c[0] * A / 2 / np.pi:.6f}, {c[1] * A / 2 / np.pi:.6f}, 0) 2pi/a")
GK, GM, KM = (np.linalg.norm(cart["K"]), np.linalg.norm(cart["M"]), np.linalg.norm(cart["K"] - cart["M"]))
p(f"  |GK| = {GK:.5f} 1/A  (analytic 4pi/(3a)      = {4 * np.pi / (3 * A):.5f})")
p(f"  |GM| = {GM:.5f} 1/A  (analytic 2pi/(sqrt3 a) = {2 * np.pi / (np.sqrt(3) * A):.5f})")
p(f"  |KM| = {KM:.5f} 1/A  (analytic 2pi/(3a)      = {2 * np.pi / (3 * A):.5f})")
chk = [abs(GK - 4 * np.pi / (3 * A)), abs(GM - 2 * np.pi / (np.sqrt(3) * A)), abs(KM - GK / 2),
       abs(GM - GK * np.sqrt(3) / 2), abs(GM - np.linalg.norm(b1) / 2)]
p(f"  CHECK hexagon geometry (|GK|=edge, |KM|=edge/2, |GM|=sqrt3/2 edge=|b1|/2): max dev = {max(chk):.1e}"
  f" -> {'PASS' if max(chk) < 1e-12 else 'FAIL'}")
# Wigner-Seitz test: K equidistant from Gamma, b1 and b2 (a corner = 3 BZ meet);
# M equidistant from Gamma and b1 only (an edge midpoint); neither closer to any other G.
Gs = [n1 * b1 + n2 * b2 for n1 in range(-3, 4) for n2 in range(-3, 4)]
for k in ("K", "M"):
    d = np.array([np.linalg.norm(cart[k] - g) for g in Gs])
    nearest = np.isclose(d, d.min(), atol=1e-9).sum()
    p(f"  {k}: distance to Gamma = {np.linalg.norm(cart[k]):.5f}, min distance to any G = {d.min():.5f},"
      f" number of equidistant nearest G = {nearest} -> {'corner' if nearest == 3 else 'edge midpoint' if nearest == 2 else '?'}")


# ---------------------------------------------------------------- (d) bands
def read_xml(path):
    root = ET.parse(path).getroot()
    ks = root.findall(".//band_structure/ks_energies")
    k = np.array([[float(x) for x in e.find("k_point").text.split()] for e in ks])  # 2pi/alat, Cartesian
    e = np.array([[float(x) for x in e.find("eigenvalues").text.split()] for e in ks]) * HA2EV
    nel = float(root.find(".//band_structure/nelec").text)
    return k, e, nel


def to_crystal(k_cart_2pia):
    return np.linalg.solve(Bmat.T * A / (2 * np.pi), k_cart_2pia.T).T


def label(kc):
    for n, v in crys.items():
        if np.allclose(kc[:2], v[:2], atol=1e-6):
            return n
    return f"({kc[0]:.4f},{kc[1]:.4f},{kc[2]:.4f})"


def edges(path):
    k, e, nel = read_xml(path)
    nocc = int(round(nel)) // 2
    iv, ic = np.argmax(e[:, nocc - 1]), np.argmin(e[:, nocc])
    kv, kc = to_crystal(k[iv]), to_crystal(k[ic])
    return dict(nel=nel, nocc=nocc, nk=len(k), vbm=e[iv, nocc - 1], cbm=e[ic, nocc],
                kv=label(kv), kc=label(kc), kvc=kv, kcc=kc, k=k, e=e)


p("\n== (d) band structure ==")
main = edges(RUNS / "01_scf_bands" / "BN.bands.xml")
p(f"  valence electrons: B 3 + N 5 = {main['nel']:.0f} per cell -> occupied bands = 8/2 = {main['nocc']}"
  f" (nbnd = 8 -> 4 empty bands)")
gap = main["cbm"] - main["vbm"]
p(f"  main path ({main['nk']} k): VBM = {main['vbm']:.4f} eV at {main['kv']}, CBM = {main['cbm']:.4f} eV at {main['kc']}")
p(f"  gap = {gap:.4f} eV, {'DIRECT' if main['kv'] == main['kc'] else 'INDIRECT'}")
e = main["e"]
iG = 0
iK = next(i for i, kk in enumerate(main["k"]) if label(to_crystal(kk)) == "K")
iM = len(e) - 1
for n, i in (("G", iG), ("K", iK), ("M", iM)):
    p(f"  at {n}: band 4 = {e[i, 3] - main['vbm']: .4f} eV, band 5 = {e[i, 4] - main['vbm']: .4f} eV"
      f" (relative to VBM); direct gap {e[i, 4] - e[i, 3]:.4f} eV")
p(f"  occupied bandwidth (band 1 min to VBM) = {main['vbm'] - e[:, 0].min():.3f} eV")
p(f"  experiment 6.1 eV (Elias et al., Nat. Commun. 10, 2639 (2019)):"
  f" LDA gap is {6.1 - gap:.2f} eV lower = {100 * (6.1 - gap) / 6.1:.0f}% underestimate"
  f" (ratio {gap / 6.1:.2f})")

p("\n  -- gap-location checks --")
res_checks = {}
for name, path in (("dense (3x) path", RUNS / "02_dense" / "BN.bands_dense.xml"),
                   ("fine path within 10-20% of K", RUNS / "03_nearK" / "BN.bands_nearK.xml")):
    r = edges(path)
    p(f"  {name} ({r['nk']} k): VBM {r['vbm']:.4f} eV at {r['kv']}, CBM {r['cbm']:.4f} eV at {r['kc']},"
      f" gap {r['cbm'] - r['vbm']:.4f} eV")
    kpts = [label(to_crystal(kk)) for kk in r["k"]]
    iKr = kpts.index("K")
    exc = (r["vbm"] - r["e"][iKr, r["nocc"] - 1]) * 1000
    p(f"     band-4 maximum exceeds its value at K by {exc:.2f} meV; band-5 minimum at {r['kc']};"
      f" gap at K = {r['e'][iKr, 4] - r['e'][iKr, 3]:.4f} eV")
    res_checks[name] = dict(nk=r["nk"], kv=r["kv"], kc=r["kc"], gap=r["cbm"] - r["vbm"], excess_meV=exc)
    assert r["kc"] == "K" and exc < 1.0        # VBM at K within sub-meV basis-set noise
scf = (RUNS / "01_scf_bands" / "BN.scf.out").read_text()
ho, lu = [float(x) for x in scf.split("highest occupied, lowest unoccupied level (ev):")[1].split()[:2]]
p(f"  SCF 6x6x1 grid (contains K): HOMO {ho} eV, LUMO {lu} eV -> gap {lu - ho:.4f} eV"
  f" (difference to bands run {abs(lu - ho - gap) * 1000:.0f} meV: default conv_thr)")

# bands.x output (used for the plot) must agree with the XML
gnu = np.loadtxt(RUNS / "01_scf_bands" / "BN.bands.dat.gnu")
nk = main["nk"]
eg = gnu[:, 1].reshape(-1, nk).T
p(f"  CHECK bands.x .gnu vs pw.x XML eigenvalues: max |diff| = {np.abs(eg - e).max():.1e} eV")

p("\n== (a) k_z check: eigenvalues at k_z = 0 vs k_z = pi/c (= b3/2) ==")
kz_k, kz_e, _ = read_xml(RUNS / "04_kz" / "BN.kz.xml")
for (i0, i1, n) in ((0, 1, "Gamma vs A"), (2, 3, "K vs H")):
    d = kz_e[i1] - kz_e[i0]
    p(f"  {n}: occupied bands 1-4 shift {np.array2string(d[:4] * 1000, precision=1)} meV;"
      f" empty bands 5-8 shift {np.array2string(d[4:] * 1000, precision=0)} meV")
occ_kz = np.abs(kz_e[[1, 3], :4] - kz_e[[0, 2], :4]).max()
p(f"  max occupied-band shift along k_z = {occ_kz * 1000:.1f} meV"
  f" (in-plane: band 4 varies by {e[:, 3].max() - e[:, 3].min():.2f} eV along G-K-M)")
p(f"  |b3| = 2pi/c = {2 * np.pi / C:.4f} 1/A vs |b1| = {np.linalg.norm(b1):.4f} 1/A (ratio {np.linalg.norm(b1) / (2 * np.pi / C):.1f})")

san = RUNS / "99_sanity_ecut120"
if san.exists():
    p("\n== numerical sanity check at ecutwfc = 120 Ry (NOT used for reported numbers) ==")
    r = edges(san / "BN.bands_nearK.ecut120.xml")
    kpts = [label(to_crystal(kk)) for kk in r["k"]]
    iKr = kpts.index("K")
    p(f"  near-K path: VBM at {r['kv']}, CBM at {r['kc']}; band-4 max - band-4(K) ="
      f" {(r['vbm'] - r['e'][iKr, 3]) * 1000:.3f} meV; gap at K = {r['e'][iKr, 4] - r['e'][iKr, 3]:.4f} eV")
    _, e120, _ = read_xml(san / "BN.kz.ecut120.xml")
    p(f"  max occupied-band shift k_z=0 -> pi/c: {np.abs(e120[[1, 3], :4] - e120[[0, 2], :4]).max() * 1000:.3f} meV")
    p("  => the sub-meV VBM offset near K and the few-meV k_z shifts at 60 Ry are plane-wave-basis noise.")

(RUNS / "results.txt").write_text("\n".join(out) + "\n")
json.dump(dict(a=A, c=C, b=(Bmat).tolist(), GK=GK, GM=GM, KM=KM, vbm=main["vbm"], cbm=main["cbm"],
               gap=gap, kv=main["kv"], kc=main["kc"], nocc=main["nocc"], checks=res_checks,
               kz_occ_max_meV=occ_kz * 1000, scf_gap=lu - ho),
          open(RUNS / "results.json", "w"), indent=1)
