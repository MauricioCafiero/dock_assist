#!/usr/bin/env python3
"""Validation suite for vina_dock's metal-site docking (the AD4Zn path).

This makes the validation record in METAL_DOCKING.md §5 reproducible from a
clean clone. That record was produced ad hoc against a receptor at
`/tmp/ad4zntest/proteinH.pdb`, a scratch path that no longer exists, so its
numbers could not be re-checked by anyone. Everything needed is built here
from the PDB entry instead.

Benchmark: **1S63**, human farnesyltransferase, the AutoDock4Zn tutorial's own
case. The catalytic Zn sits at (18.142, 132.126, 5.224) coordinated by
Asp B297 / Cys B299 / His B362, and the co-crystallized inhibitor `778` is a
nitrile-bearing chelator that should redock with a nitrogen on the zinc.

  METAL_DOCKING.md §5 expects: engine `ad4zn`, best ≈ -12 kcal/mol,
  "COORDINATES ZN at 1.91 A" (crystal Zn-N ≈ 2.0-2.2 A).

Receptor construction: protein `ATOM` records plus the ZN ion, dropping
everything else -- the inhibitor being redocked (so it cannot bias the
result), the FPP cofactor, sucrose and waters. The AD4Zn path hands the
receptor to meeko, which templates Zn but not FPP/sucrose, so this is also
what that path can actually consume.

Checks, fast ones first:

  detect    detect_metal_sites finds the Zn at the recorded coordinates with
            its recorded coordinating residues. No docking, seconds.
  error     a box centre with no ion nearby raises a DockError that points at
            detect_metal_sites_agent. No docking, seconds.
  redock    the full AD4Zn pipeline: meeko receptor prep -> zinc_pseudo TZ
            atoms -> per-ligand AutoGrid -> Vina --scoring ad4. Reports score,
            donor-to-zinc distance, and pose agreement with the crystal
            inhibitor measured in the crystal frame with no superposition.
            Minutes, and needs the adfr-suite and vina conda envs.

Usage:
  python3 test_metal.py              # all checks
  python3 test_metal.py --fast       # detect + error only, no docking
  python3 test_metal.py --mutate-mg  # also redock with ZN->MG (tier engine)
"""
import argparse
import math
import os
import sys
import urllib.request

import vina_dock as vd

HERE = os.path.dirname(os.path.abspath(__file__))
TESTDIR = os.path.join(HERE, "..", "pdb_files", "metal_tests")

PDB_ID = "1S63"
ZN_XYZ = (18.142, 132.126, 5.224)
ZN_CONTACTS = {"ASP B297", "CYS B299", "HIS B362"}
LIGAND_HET = "778"
# inhibitor SMILES as recorded in METAL_DOCKING.md §5
LIGAND_SMILES = "N#Cc1ccc(Cn2cncc2CN2CCN(c3cccc(Cl)c3)C(=O)C2)cc1"
EXPECT_SCORE = -12.0      # kcal/mol, ad4zn engine
EXPECT_COORD = 1.91       # Angstroms, donor N to Zn


def fetch_receptor(mutate_mg=False):
    """Download 1S63 and build the docking receptor: protein + the metal ion.

    Everything else is dropped: the inhibitor `778` (it is what we redock, so
    leaving it in would both block the site and bias the result), the FPP
    cofactor, sucrose and waters. With `mutate_mg` the ZN record is rewritten
    as MG, which is how METAL_DOCKING.md exercised a non-zinc engine on a
    receptor whose geometry is otherwise identical.
    """
    os.makedirs(TESTDIR, exist_ok=True)
    raw = os.path.join(TESTDIR, f"{PDB_ID}.pdb")
    if not os.path.exists(raw):
        print(f"  fetching {PDB_ID} from RCSB ...")
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{PDB_ID}.pdb", raw)
    tag = "mg" if mutate_mg else "zn"
    out = os.path.join(TESTDIR, f"{PDB_ID.lower()}_{tag}.pdb")
    n_prot = n_metal = 0
    with open(raw) as src, open(out, "w") as dst:
        for line in src:
            if line.startswith("ATOM"):
                dst.write(line)
                n_prot += 1
            elif line.startswith("HETATM") and line[17:20].strip() == "ZN":
                if mutate_mg:
                    line = (line[:17] + " MG" + line[20:76] + "MG" + line[78:])
                dst.write(line)
                n_metal += 1
        dst.write("END\n")
    print(f"  {os.path.basename(out)}: {n_prot} protein atoms + {n_metal} "
          f"{'MG' if mutate_mg else 'ZN'} ion")
    return out, raw


def crystal_ligand(raw_pdb):
    """Heavy atoms of the co-crystallized inhibitor, for pose comparison."""
    out = []
    with open(raw_pdb) as fh:
        for line in fh:
            if line.startswith("HETATM") and line[17:20].strip() == LIGAND_HET:
                el = (line[76:78].strip() or line[12:16].strip()[0]).upper()
                if el != "H":
                    out.append((el, (float(line[30:38]), float(line[38:46]),
                                     float(line[46:54]))))
    return out


def check_detect(receptor):
    """detect_metal_sites finds the Zn where METAL_DOCKING.md says it is."""
    print("\n--- detect: metal-site detection -------------------------------")
    sites = vd.detect_metal_sites(receptor)
    zn = [s for s in sites if s["element"] == "ZN"]
    if not zn:
        print(f"  FAIL: no ZN found (sites: {[s['element'] for s in sites]})")
        return False
    s = zn[0]
    d = math.dist(s["xyz"], ZN_XYZ)
    got = set(s["contacts"])
    ok_pos = d <= 0.01
    ok_con = ZN_CONTACTS <= got
    print(f"  ZN at ({s['xyz'][0]:.3f}, {s['xyz'][1]:.3f}, {s['xyz'][2]:.3f})"
          f"  [{d:.3f} A from the recorded position]  {'OK' if ok_pos else 'FAIL'}")
    print(f"  contacts: {', '.join(sorted(got)) or 'none'}")
    missing = ZN_CONTACTS - got
    print(f"  expected {', '.join(sorted(ZN_CONTACTS))}: "
          f"{'OK' if ok_con else 'MISSING ' + ', '.join(sorted(missing))}")
    return ok_pos and ok_con


def check_error_path(receptor):
    """A centre with no ion nearby fails cleanly, pointing at the right tool."""
    print("\n--- error: box centre with no metal ion -------------------------")
    far = (ZN_XYZ[0] + 60.0, ZN_XYZ[1], ZN_XYZ[2])
    try:
        vd.dock_metal_site(receptor, [LIGAND_SMILES], center=far)
    except vd.DockError as e:
        msg = str(e)
        ok = "detect_metal_sites_agent" in msg
        print(f"  DockError raised: {msg.splitlines()[0]}")
        print(f"  points at detect_metal_sites_agent: {'OK' if ok else 'FAIL'}")
        return ok
    except Exception as e:
        print(f"  FAIL: raised {type(e).__name__} instead of DockError: {e}")
        return False
    print("  FAIL: no error raised for a centre with no ion")
    return False


def check_redock(receptor, raw_pdb, center, expect_engine, label):
    """Full pipeline redock, compared with the crystal pose and the record."""
    print(f"\n--- redock: {label} ---------------------------------------")
    report = vd.dock_metal_site(receptor, [LIGAND_SMILES], center=center,
                                exhaustiveness=16, seed=1)
    print(report)

    engine = None
    for line in report.splitlines():
        if line.strip().startswith("engine:"):
            engine = line.split(":", 1)[1].split()[0]
    score = coord = None
    for line in report.splitlines():
        t = line.strip()
        if t.startswith("score:"):
            try:
                score = float(t.split(":")[1].split()[0])
            except (IndexError, ValueError):
                pass
        if "COORDINATES" in t:
            try:
                coord = float(t.split("at")[1].split()[0])
            except (IndexError, ValueError):
                pass

    print("\n  checks:")
    ok_engine = engine == expect_engine
    print(f"    engine: {engine} (expected {expect_engine}) "
          f"{'OK' if ok_engine else 'FAIL'}")
    if score is None:
        print("    score: not parsed -- FAIL")
        return False
    print(f"    score: {score:+.2f} kcal/mol "
          f"(record: {EXPECT_SCORE:+.1f} for ad4zn)")
    if coord is None:
        print("    coordination: pose does not coordinate the ion")
    else:
        print(f"    coordination: {coord:.2f} A "
              f"(record: {EXPECT_COORD:.2f} A; crystal Zn-N 2.0-2.2 A)")

    # pose agreement with the crystal inhibitor, in the crystal frame.
    # Deliberately graph-free (centroid + overlap) rather than an RMSD: see
    # test_covalent.pose_metrics for why a PDBQT round-trip makes atom-level
    # matching unreliable.
    sdf = os.path.splitext(os.path.abspath(receptor))[0] + "_m0.sdf"
    ref = crystal_ligand(raw_pdb)
    if os.path.exists(sdf) and ref:
        from rdkit import Chem
        mol = next(iter(Chem.SDMolSupplier(sdf, removeHs=True)), None)
        if mol is not None:
            conf = mol.GetConformer()
            pose = [tuple(conf.GetAtomPosition(i))
                    for i in range(mol.GetNumAtoms())]
            cd = [sum(p[k] for p in pose) / len(pose) for k in range(3)]
            cx = [sum(r[1][k] for r in ref) / len(ref) for k in range(3)]
            near = sum(1 for _el, x in ref
                       if any(math.dist(x, p) <= 2.0 for p in pose))
            print(f"    pose vs crystal {LIGAND_HET}: centroid "
                  f"{math.dist(cd, cx):.2f} A, "
                  f"heavy-atom overlap {100 * near / len(ref):.0f}% "
                  f"({len(pose)} docked vs {len(ref)} crystal atoms)")
    else:
        print(f"    pose SDF not found at {sdf}")
    return ok_engine and score is not None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fast", action="store_true",
                    help="detection and error checks only, no docking")
    ap.add_argument("--mutate-mg", action="store_true",
                    help="also redock a ZN->MG mutant (exercises a tier engine)")
    args = ap.parse_args()

    print("Metal-site docking validation (1S63, AutoDock4Zn benchmark)")
    receptor, raw = fetch_receptor()

    results = {"detect": check_detect(receptor),
               "error": check_error_path(receptor)}
    if not args.fast:
        results["redock (ad4zn)"] = check_redock(
            receptor, raw, ZN_XYZ, "ad4zn", "Zn site, AD4Zn engine")
        if args.mutate_mg:
            mg_rec, _ = fetch_receptor(mutate_mg=True)
            results["redock (ad4map/Mg)"] = check_redock(
                mg_rec, raw, ZN_XYZ, "ad4map", "ZN->MG mutant, tier-0 engine")

    print(f"\n{'=' * 64}\nSUMMARY\n{'=' * 64}")
    for name, ok in results.items():
        print(f"  {name:22} {'PASS' if ok else 'FAIL'}")
    return 0 if all(results.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
