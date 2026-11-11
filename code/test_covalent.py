#!/usr/bin/env python3
"""Redocking validation for vina_dock's covalent (tethered) docking engine.

Four classic covalent complexes, each with a crystallographic LINK record
giving the real protein-ligand bond, are redocked from the FREE drug SMILES
and the predicted pose is compared with the crystal adduct:

  6OIM  KRAS G12C      sotorasib    Cys A12   acrylamide    SG-C25 1.81 A
  7RFS  SARS-CoV-2 Mpro nirmatrelvir Cys A145 nitrile       SG-C3  1.80 A
  4G5J  EGFR           afatinib     Cys A797  crotonamide   SG-C30 1.84 A
  5P9J  BTK            ibrutinib    Cys A481  acrylamide    SG-CAA 1.85 A

Each case tests a different part of the engine: the Michael-acceptor and
nitrile warhead reactions, and four unrelated binding sites. The receptors are
the protein-only chain A of each entry (all HETATMs dropped, so the crystal
ligand itself cannot bias the result).

Metrics reported per case, all computed in the crystal frame with NO
superposition -- the tether puts the pose in the crystal's own coordinates, so
any agreement is real and not fitted:

  adduct score   Vina score of the refined adduct (ranks a series only)
  drift          anchor-atom RMSD after Vina's local refinement; the check
                 that the covalent attachment survived
  warhead err    distance from the docked reacting carbon to the crystal's
                 LINK atom. The tether fixes this atom's distance from SG but
                 not its direction, so this measures the predicted geometry.
  centroid err   docked vs crystal ligand centroid
  overlap        fraction of crystal heavy atoms with a docked heavy atom
                 within 2 A -- a graph-free shape agreement measure
  rmsd           symmetry-corrected in-place RMSD when the docked and crystal
                 molecular graphs can be matched, else "n/a"

Usage:
  python3 test_covalent.py                 # all cases, default conformers
  python3 test_covalent.py --case 6OIM     # one case
  python3 test_covalent.py --n-confs 50    # faster, less thorough
Receptors/structures are expected in ../pdb_files/covalent_tests (fetch_pdbs()
downloads any that are missing).
"""
import argparse
import math
import os
import sys
import urllib.request

import vina_dock as vd

HERE = os.path.dirname(os.path.abspath(__file__))
TESTDIR = os.path.join(HERE, "..", "pdb_files", "covalent_tests")

# (pdb id, chain, resseq, drug name, free-drug SMILES, HET code, LINK atom,
#  crystal bond length)
CASES = [
    ("6OIM", "A", 12, "sotorasib",
     "C[C@H]1CN(CCN1C2=NC(=O)N(C3=NC(=C(C=C32)F)C4=C(C=CC=C4F)O)"
     "C5=C(C=CN=C5C(C)C)C)C(=O)C=C", "MOV", "C25", 1.81),
    ("7RFS", "A", 145, "nirmatrelvir",
     "CC1([C@@H]2[C@H]1[C@H](N(C2)C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F)"
     "C(=O)N[C@@H](C[C@@H]3CCNC3=O)C#N)C", "4WI", "C3", 1.80),
    ("4G5J", "A", 797, "afatinib",
     "CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)"
     "O[C@H]4CCOC4", "0WN", "C30", 1.84),
    ("5P9J", "A", 481, "ibrutinib",
     "C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N",
     "8E8", "CAA", 1.85),
]


def fetch_pdbs(cases):
    """Download any missing test structures and write protein-only receptors."""
    os.makedirs(TESTDIR, exist_ok=True)
    for pid, chain, *_ in cases:
        raw = os.path.join(TESTDIR, f"{pid}.pdb")
        if not os.path.exists(raw):
            url = f"https://files.rcsb.org/download/{pid}.pdb"
            print(f"  fetching {pid} from RCSB ...")
            urllib.request.urlretrieve(url, raw)
        prot = os.path.join(TESTDIR, f"{pid}_prot.pdb")
        if not os.path.exists(prot):
            n = 0
            with open(raw) as src, open(prot, "w") as dst:
                for line in src:
                    if line.startswith("ATOM") and line[21] == chain:
                        dst.write(line)
                        n += 1
                dst.write("END\n")
            print(f"  wrote {os.path.basename(prot)} ({n} protein atoms)")


def crystal_reference(pid, chain, resseq, het, link_atom):
    """Crystal geometry for one case: anchor, ligand atoms and the LINK atom."""
    raw = os.path.join(TESTDIR, f"{pid}.pdb")
    anchor, lig, link = {}, [], None
    with open(raw) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            name = line[12:16].strip()
            if line.startswith("ATOM") and line[22:26].strip() == str(resseq):
                if name in ("CA", "CB", "SG"):
                    anchor[name] = xyz
            elif line.startswith("HETATM") and line[17:20].strip() == het:
                el = (line[76:78].strip() or name[0]).upper()
                if el != "H":
                    lig.append((el, xyz))
                if name == link_atom:
                    link = xyz
    return anchor, lig, link


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def pose_metrics(pose_sdf, resname, lig_ref, link_ref):
    """Compare a docked covalent pose with the crystal adduct, in place.

    Returns a dict of the metrics documented in the module docstring. The
    docked pose still carries the three side-chain mimic atoms, which are
    stripped before comparison because the crystal HET record does not
    include them (the S belongs to the cysteine, not the ligand).
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    mol = next(iter(Chem.SDMolSupplier(pose_sdf, removeHs=True)), None)
    if mol is None:
        return {"error": f"could not read {pose_sdf}"}
    # Locate the side-chain mimic. Poses exported by meeko's mk_export.py
    # carry proper bond orders and hydrogens, so the engine's own
    # H-counting pattern works; a pose that came back through a bare PDBQT
    # round-trip has lost its nonpolar hydrogens, and only the heavy-atom
    # DEGREE pattern can match there. Try both.
    hits = ()
    for mimic_smarts in (vd.COVALENT_RESIDUES[resname][2],
                         {"CYS": "[CD1][CD2][SD2]",
                          "SER": "[CD1][CD2][OD2]"}[resname]):
        hits = mol.GetSubstructMatches(Chem.MolFromSmarts(mimic_smarts))
        if len(hits) == 1:
            break
    if len(hits) != 1:
        return {"error": f"mimic {mimic_smarts} matched {len(hits)}x in pose"}
    ca, cb, nu = hits[0]
    conf = mol.GetConformer()

    # the reacting carbon is the mimic nucleophile's other neighbour: that is
    # the atom the crystal LINK record points at
    warhead = None
    for nb in mol.GetAtomWithIdx(nu).GetNeighbors():
        if nb.GetIdx() != cb:
            warhead = nb.GetIdx()
    out = {}
    if warhead is not None and link_ref is not None:
        p = conf.GetAtomPosition(warhead)
        out["warhead_err"] = dist((p.x, p.y, p.z), link_ref)

    # ligand-only heavy atoms of the docked pose (mimic removed)
    keep = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetIdx() not in (ca, cb, nu) and a.GetAtomicNum() > 1]
    dpos = []
    for i in keep:
        p = conf.GetAtomPosition(i)
        dpos.append((mol.GetAtomWithIdx(i).GetSymbol().upper(), (p.x, p.y, p.z)))
    out["n_docked"] = len(dpos)
    out["n_crystal"] = len(lig_ref)
    if dpos and lig_ref:
        cd = [sum(q[1][k] for q in dpos) / len(dpos) for k in range(3)]
        cx = [sum(q[1][k] for q in lig_ref) / len(lig_ref) for k in range(3)]
        out["centroid_err"] = dist(cd, cx)
        near = 0
        for el, x in lig_ref:
            if any(dist(x, y) <= 2.0 for _, y in dpos):
                near += 1
        out["overlap"] = near / len(lig_ref)

    # symmetry-corrected in-place RMSD, when the graphs can be matched
    out["rmsd"] = None
    try:
        edit = Chem.RWMol(mol)
        for i in sorted((ca, cb, nu), reverse=True):
            edit.RemoveAtom(i)
        stripped = edit.GetMol()
        Chem.SanitizeMol(stripped)
        ref_block = ["HETATM%5d %-4s LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00"
                     "          %2s" % (i + 1, f"{el}{i + 1}", x[0], x[1], x[2], el)
                     for i, (el, x) in enumerate(lig_ref)]
        ref = Chem.MolFromPDBBlock("\n".join(ref_block) + "\nEND\n",
                                   removeHs=True, proximityBonding=True)
        if ref is not None and ref.GetNumAtoms() == stripped.GetNumAtoms():
            from rdkit.Chem import AllChem
            ref = AllChem.AssignBondOrdersFromTemplate(stripped, ref)
            out["rmsd"] = rdMolAlign.CalcRMS(stripped, ref)
    except Exception as e:
        out["rmsd_note"] = f"{type(e).__name__}"
    return out


def run_case(case, ga_runs):
    pid, chain, resseq, drug, smiles, het, link_atom, bond = case
    # Dock against the RAW PDB, not the protein-only copy. The engine does its
    # own HETATM handling (keep cofactors and ions, drop waters and the ligand
    # bound at the covalent site), so feeding it a pre-stripped receptor would
    # bypass exactly the code that needs testing -- and would not match what
    # the agent passes in, which is whatever get_protein_from_pdb fetched.
    # Each of the four cases carries a LINK record to its target residue, so
    # this also checks the LINK-based removal of the bound adduct.
    prot = os.path.join(TESTDIR, f"{pid}.pdb")
    print(f"\n{'=' * 72}\n{pid}  {drug}  ->  Cys {chain}{resseq}   "
          f"(crystal SG-{link_atom} = {bond:.2f} A)\n{'=' * 72}")
    anchor, lig_ref, link_ref = crystal_reference(pid, chain, resseq, het,
                                                  link_atom)
    # the engine names the pose <stem>_cov0.sdf; delete any pose left by an
    # earlier run FIRST. Without this a failed docking silently validates the
    # previous run's file and reports its metrics as though they were new --
    # which is exactly what happened on 7RFS.
    pose = os.path.splitext(os.path.abspath(prot))[0] + "_cov0.sdf"
    if os.path.exists(pose):
        os.remove(pose)
    try:
        report = vd.dock_covalent(prot, [smiles], chain=chain, resseq=resseq,
                                  ga_runs=ga_runs)
    except Exception as e:
        print(f"  !! docking raised {type(e).__name__}: {e}")
        return {"pdb": pid, "drug": drug, "ok": False,
                "note": f"{type(e).__name__}"}
    print(report)
    if not os.path.exists(pose):
        print("  !! no pose written -- docking failed for this case")
        return {"pdb": pid, "drug": drug, "ok": False}
    # take the residue type from the structure rather than assuming Cys, so a
    # serine case can be added to CASES without touching this
    target_resname, _ = vd._residue_anchor(prot, chain, resseq)
    m = pose_metrics(pose, target_resname, lig_ref, link_ref)
    m.update({"pdb": pid, "drug": drug, "ok": "error" not in m})
    print("\n  validation vs crystal:")
    if "error" in m:
        print(f"    error: {m['error']}")
    else:
        print(f"    warhead atom error : "
              f"{m.get('warhead_err', float('nan')):.2f} A")
        print(f"    centroid error     : "
              f"{m.get('centroid_err', float('nan')):.2f} A")
        print(f"    heavy-atom overlap : {m.get('overlap', 0) * 100:.0f}%  "
              f"({m.get('n_docked')} docked vs {m.get('n_crystal')} crystal "
              f"atoms)")
        print(f"    symmetry RMSD      : "
              f"{('%.2f A' % m['rmsd']) if m.get('rmsd') else 'n/a'}")
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", action="append",
                    help="PDB id to run (repeatable); default all")
    ap.add_argument("--ga-runs", type=int, default=None,
                    help=f"AutoDock4 GA runs per ligand (engine default "
                         f"{vd.COVALENT_GA_RUNS}); lower is faster")
    args = ap.parse_args()

    cases = CASES
    if args.case:
        want = {c.upper() for c in args.case}
        cases = [c for c in CASES if c[0] in want]
        if not cases:
            sys.exit(f"no such case; known: {', '.join(c[0] for c in CASES)}")

    print("Covalent docking validation")
    fetch_pdbs(cases)
    rows = [run_case(c, args.ga_runs) for c in cases]

    print(f"\n{'=' * 72}\nSUMMARY\n{'=' * 72}")
    print(f"{'PDB':6} {'drug':14} {'warhead':>9} {'centroid':>9} "
          f"{'overlap':>8} {'rmsd':>7}")
    for r in rows:
        if not r.get("ok"):
            print(f"{r['pdb']:6} {r['drug']:14} {'FAILED':>9}")
            continue
        print(f"{r['pdb']:6} {r['drug']:14} "
              f"{r.get('warhead_err', float('nan')):8.2f}A "
              f"{r.get('centroid_err', float('nan')):8.2f}A "
              f"{r.get('overlap', 0) * 100:7.0f}% "
              f"{('%.2f' % r['rmsd']) if r.get('rmsd') else 'n/a':>7}")


if __name__ == "__main__":
    main()
