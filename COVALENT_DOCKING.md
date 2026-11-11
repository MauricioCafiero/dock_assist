# Covalent docking — engineering notes

Dev notes for the covalent docking added to `code/vina_dock.py`. The
user-facing summary lives in `README.md`; everything here is implementation
detail, measured behaviour and open work. Same split as
[`METAL_DOCKING.md`](METAL_DOCKING.md).

Added 2026-11-11.

## What was built

Three agent-facing tools, plus the engine behind them:

| Tool | Purpose |
|---|---|
| `detect_nucleophiles_agent(receptor_pdb, center_x, center_y, center_z)` | List nucleophilic side chains (Cys SG, Ser OG, Thr OG1, Lys NZ, Tyr OH, His NE2) near a site, flagging which are dockable |
| `dock_covalent_agent(receptor_pdb, smiles_list, chain, resseq)` | Covalent docking onto a Cys/Ser, via AutoDock4 |
| `check_warhead_geometry_agent(receptor_pdb, pose_sdf, smiles)` | Post-hoc: is an ordinary docked pose poised to react? |

The intended agent flow mirrors the metal path:

```
find_PDBID_node -> get_protein_from_pdb -> find_ligand_centroids
   -> detect_nucleophiles_agent(centroid)        # which Cys/Ser line this site?
   -> dock_covalent_agent(chain, resseq)         # covalent docking
```

`check_warhead_geometry_agent` is the cheap alternative entry point: dock
normally with `blind_dock_agent` / `dock_at_centroid_agent`, then ask whether
the warhead happened to land within striking distance of a nucleophile.

## The engine

This is the **official AutoDock "flexible side chain" covalent docking
protocol** (Bianco, Forli, Goodsell & Olson, *Protein Sci* 2015; the
`adCovalentDockResidue` tutorial from
[autodock.scripps.edu/resources/covalent-docking](https://autodock.scripps.edu/resources/covalent-docking/)),
not a local invention.

The ligand's warhead is reacted in silico with a mimic of the target side
chain, and the resulting **adduct is docked as a flexible residue with no
ligand at all**. AutoDock4 pins a flexible residue's root atom to its input
coordinates, so the covalent attachment physically cannot move during the
search, while the whole adduct's torsions get a full Lamarckian GA search.
That is what makes this a real docking rather than a constrained rescoring.

| Step | Function |
|---|---|
| 1. React the warhead with a side-chain mimic (CA-CB-SG for Cys, CA-CB-OG for Ser) | `_build_covalent_adduct` |
| 2. Rigid receptor with the target side chain removed | `_cov_rigid_receptor` |
| 3. Adduct superimposed on the residue, written as a `BEGIN_RES` flexible residue | `_cov_flexres_ligand` |
| 4. Grid box over the adduct's reachable volume, then AutoGrid maps | `_cov_envelope`, `_build_cov_gpf`, `_run_autogrid` |
| 5. The search | `_build_cov_dpf`, `_run_autodock4` |
| 6. Energies and pose | `_parse_cov_dlg`, `_cov_export_sdf` |

Scores are ordinary AutoDock4 binding free energies of the adduct, with a Ki,
and so **are comparable between ligands**. They are not comparable to
`blind_dock` / `dock_at_centroid` affinities, which come from Vina's own
function. The geometric readout is the **anchor drift**, reported per pose —
the covalent counterpart of the metal path's coordination check.

### The three DPF lines that carry the whole trick

All three are easy to miss and two of them cost real debugging time:

```
move empty             # there is NO real ligand; `empty` is a zero-byte file
flexres <adduct>.pdbqt # the covalent adduct enters here instead
unbound_energy 0.0     # replaces the default `unbound_model bound`
```

- **`move empty`** — a literal zero-byte file (`touch empty`). The tutorial
  bundle ships one. Omitting `move` entirely makes AutoDock4 run and report
  `Successful Completion`, but it is a null docking.
- **`unbound_energy 0.0`** — **[measured] essential.** The default
  `unbound_model bound` sets the unbound reference equal to the bound energy,
  so every pose scores exactly `+0.00 kcal/mol` and the RANKING table comes
  back as `nan`. This was the symptom that initially looked like "AD4 can't do
  flexres-only docking".
- **`rmsdatoms all`**, which the tutorial README suggests adding so flexres
  atoms count toward RMSD/clustering, is **not** written: [measured]
  AutoDock4 4.2.6 rejects it (`FATAL ERROR: DPF> rmsdatoms all / Unrecognized
  keyword`) and the tutorial's own shipped DPF omits it too. It affects only
  RMSD reporting, not energies.

### Install

bioconda ships AutoDock4 for **osx-64 only**, so it runs under Rosetta like
the existing `adfr-suite` env. It was given its own env rather than added to
`adfr-suite`, whose `autogrid4` was patched in place for AD4Zn and should not
be disturbed:

```sh
CONDA_SUBDIR=osx-64 conda create -n autodock4 -c bioconda autodock
```

Overridable with `AUTODOCK4_BIN`; AutoGrid comes from the existing
`AUTOGRID_BIN`.

[measured] `autodock4` itself needs no patching, unlike `adfr-suite`'s
`autogrid4`: `otool -L` shows it links only `/usr/lib/libSystem.B.dylib` and
`@rpath/libc++.1.dylib`. The env also ships **`autodock4.omp`**, an OpenMP
build, which is deliberately **not** used — it links `@rpath/libomp.dylib`,
and libomp is keg-only on this machine, so it would reintroduce a library
discovery problem for no benefit (the plain binary is single-threaded, and the
GA search dominates the runtime anyway).

[measured] Both binaries run with **no environment at all** — verified under
`env -i` with only `HOME` and `PATH=/usr/bin:/bin`. No conda activation, no
`DYLD_FALLBACK_LIBRARY_PATH`, no `~/.zshrc`. `autogrid4` resolves its patched
libgomp through `@executable_path/../lib/`, which is what the §2.3 repair in
METAL_DOCKING.md buys: the install name is relative to the binary rather than
absolute. That matters because `dock_assist.py` invokes these as subprocesses
from a venv, which inherits none of a login shell's exports. `agfr`/`adfr` from the ADFR suite were also considered (ADFR
can do covalent docking) but **[measured] both segfault on this machine**
(exit 139, no output); only `autogrid4` and `pythonsh` from that env work.

## Measured gotchas

- **The ligand already bound at the covalent site must be removed — but
  nothing else.** Removing it is a correctness requirement, not tidying: a
  structure fetched for redocking still holds the co-crystallized inhibitor,
  and in a covalent complex it occupies the very pocket being docked into. The
  official tutorial README is blunt — *"if a covalent ligand is already bound
  to the residue in the protein structure, it must be removed!"*. It is also
  what meeko needs: **[measured]** on the raw 6OIM as fetched by
  `get_protein_from_pdb`, `mk_prepare_receptor.py` aborts with `Input residues
  {'A:302': 'GDP', 'A:303': 'MOV'} not in residue_templates`, `MOV` being the
  bound sotorasib adduct.

  The first fix for that dropped **all** HETATM, which was a regression: the
  normal docking path (`build_receptor_pdbqt`, plain `obabel -h -xr`) keeps
  every heteroatom, and emptying the structure changes the shape and
  electrostatics of the pocket. `_cov_classify_hetatm` now keeps cofactors,
  metal ions and unrelated ligands, and drops only waters plus the group bound
  at the site — identified from a `LINK` record tying it to the target
  residue, falling back to any group within `COVALENT_SITE_LIGAND_CUTOFF` of
  the nucleophile.

  **That cutoff is a bond length (2.5 Å), not a contact distance, and the
  difference matters:** [measured] in 6OIM the covalently bound `MOV` is
  1.81 Å from Cys12 SG while the `GDP` cofactor that must be kept is only
  **4.34 Å** away, because the switch-II pocket abuts the nucleotide site. A
  5 Å rule — the obvious first guess — would have destroyed GDP.

  Since meeko cannot be handed the heteroatoms at all, they are typed
  separately by Open Babel and appended to meeko's rigid PDBQT
  (`_cov_append_hetatm`) — the same division of labour the metal path uses for
  ions. `_build_cov_gpf` reads `receptor_types` straight from that file, so the
  extra types reach AutoGrid automatically; verified on 6OIM, whose rigid
  receptor ends up with `P x2` and `Mg x1` alongside the protein types, both
  standard AD4 types needing no custom parameter file.

  The HETATM problem surfaced only in the end-to-end agent run, because
  `test_covalent.py` pre-strips its receptors. Relying on the caller to supply
  a clean PDB was the original bug; over-correcting to protein-only was the
  second.
- **The agent wrappers return failures as text, not exceptions.**
  `dock_assist.chat_turn` calls tool functions directly with no try/except, so
  a raised `DockError` ends the whole chat session instead of letting the model
  react. All three covalent wrappers catch `DockError` and return the message,
  which is also more useful to an LLM (it can pick another residue). The
  underlying `dock_covalent` still raises, for direct Python callers and the
  test harness.
- **meeko takes two different residue spellings.**
  `mk_prepare_receptor.py -f` wants `"<chain>:<resnum>"` (`A:12`) and raises
  `ValueError: missing residue in []` on anything else, while the covalent
  ligand path's `--rec_residue` wants `"<chain>:<resname>:<resnum>"`
  (`A:CYS:12`). Both are used, which is why it is called out in the code.
- **[measured] `mk_prepare_receptor.py -f <chain>:<resnum> -p` strips CA, CB
  and the whole side chain from the `_rigid.pdbqt`**, leaving only that
  residue's backbone `N/H/C/O`; the side chain reappears in the companion
  `_flex.pdbqt` with `ROOT`=CA → `BRANCH` CB → SG. That is exactly what this
  protocol needs — the adduct's mimic atoms occupy the vacated positions — so
  `_cov_rigid_receptor` keeps the rigid file and **deletes `_flex.pdbqt`**,
  because the adduct replaces the residue's own side chain.
- **meeko refuses structures with alternate locations**
  (`Error: Creation of data structure for receptor failed ... Residues with
  alternate location: ['A:970']` on 4G5J) and **aborts on any residue that
  fails template matching** (`No template matched for residue_key='A:640'` on
  5P9J/BTK, surface Glu/Ser with incomplete side chains). Hence
  `--default_altloc A --allow_bad_res`. Safe only because `_residue_anchor`
  validates the *target* residue against the PDB first, so a truncated target
  is reported rather than silently deleted — but the rigid receptor can be
  missing odd surface residues.
- **The tether indices are `(CA, CB)` with CA terminal.** meeko's
  `CovalentBuilder` hardcodes `force_CA_CB=True`, and
  `root_atom_index = cov_lig.indices[0]`, so index 0 must be the CA mimic.
  The tutorial's own `tether_smarts_indices 5 4` on `n1cc(CC)nc1` for His
  confirms the ordering. Here the mimic is `...S-CH2-CH3` and
  `[CH3][CH2][SX2]` with indices `1 2` gives (CA, CB).
- **A pose read back from a bare PDBQT has lost its chemistry.** Nonpolar
  hydrogens are merged and bond orders are gone: RDKit reads the carbons as
  bracketed `[C]` with zero implicit H. So a Michael acceptor — *defined* by
  its C=C — cannot be perceived there, and `AssignBondOrdersFromTemplate`
  restores bond orders but not hydrogen counts, so every `[CX3]` query still
  fails. Consequences:
  - `check_warhead_geometry` takes the ligand **SMILES** as the source of
    chemistry and maps it onto the pose through `_flat_graph` (both reduced to
    element/connectivity, which *does* survive). Verified 41/41 atoms matched.
  - poses out of the engine go through meeko's `mk_export.py`, which
    reconstructs proper bond orders and hydrogens from the DLG.
- **`mk_export.py` defaults to cluster leads** and fails with
  `RuntimeError: no cluster_leads in pdbqt_mol but only_cluster_leads=True` on
  a DLG from an unfinished or unclustered run. `_cov_export_sdf` falls back to
  `--all_dlg_poses`.
- **Grid box size is not a free parameter.** The tethered adduct sweeps a
  volume around its anchor, and a box that fails to cover it silently truncates
  the search (AD4 charges `extnrg` for leaving). `_cov_envelope` embeds
  conformers purely to measure that reach and pads it; the edge is capped at
  126 grid points, AutoGrid's limit.

## Validation

`code/test_covalent.py` redocks four classic covalent complexes from the
**free drug SMILES** and compares the predicted pose with the crystal adduct.
Each entry has a crystallographic `LINK` record giving the real bond. The
harness hands the engine the **raw PDB**, so the receptor is prepared exactly
as in normal use: cofactors and metal ions retained, waters dropped, and the
bound adduct removed via its `LINK` record so it cannot bias the result. That
is also the physically correct model of the site — KRAS's GDP and Mg really
are there — and it is what the numbers below were measured on.

| PDB | Target | Drug | Residue | Warhead | Crystal bond |
|---|---|---|---|---|---|
| 6OIM | KRAS G12C | sotorasib | Cys A12 | acrylamide | SG–C25 1.81 Å |
| 7RFS | SARS-CoV-2 Mpro | nirmatrelvir | Cys A145 | nitrile | SG–C3 1.80 Å |
| 4G5J | EGFR | afatinib | Cys A797 | crotonamide | SG–C30 1.84 Å |
| 5P9J | BTK | ibrutinib | Cys A481 | acrylamide | SG–CAA 1.85 Å |

Structures not in the test tree are fetched from RCSB by `fetch_pdbs()`
(`pdb_files/` is `.gitignore`d, so they are never committed).

Two candidate cases were **rejected** while assembling this set, worth
recording so nobody re-adds them:

- **6LU7** (Mpro + N3) — the inhibitor is a multi-residue HET chain
  (`02J-ALA-VAL-LEU-PJE-010`), so it cannot be read as one ligand.
- **4ZAU** (EGFR + osimertinib) — **has no `LINK` record**: osimertinib is
  bound non-covalently there, so it is not a covalent ground truth. 4G5J
  (afatinib) was used for EGFR instead.

Metrics are computed **in the crystal frame with no superposition**: the
tether puts the pose in the crystal's own coordinates, so any agreement is
real rather than fitted. Judge poses on **centroid error and overlap** — the
warhead atom error is flattering on its own, because the tether forces the
reacting carbon to within ~1.8 Å of the nucleophile by construction.

### Chemistry check

All warhead reactions give the chemically correct adduct, and the
regiochemistries match the crystal HET records:

- sotorasib → `CCSCCC(=O)N…` — S on the terminal CH₂ (thiol-Michael, β-addition)
- nirmatrelvir → `CCSC(=N)…` — thioimidate, matching 7RFS's HETNAM "1-imino"
- afatinib → `CCSC(CC(=O)N…)CN(C)C` — S on the β carbon, giving the
  4-(dimethylamino)butanamide that 4G5J's `0WN` HETNAM describes
- ibrutinib → `CCSCCC(=O)N1CCC…` — S on the terminal CH₂
- a chloroacetamide probe → `CCSCC(=O)Nc1ccccc1` (SN2, halide dropped)
- a warhead-free phenol → clean `DockError` naming the recognised warheads

RDKit prints `mapped atoms in the reactants were not mapped in the products:
4` when the haloacetamide template is built — that is the leaving halide,
dropped on purpose. It is suppressed on the quiet path.

### Redocking results

From `python3 code/test_covalent.py` at the tutorial's settings
(`COVALENT_GA_RUNS = 10`, `COVALENT_GA_EVALS = 2500000`), 2026-11-11. All four
cases dock, all from the free drug SMILES against a HETATM-stripped receptor:

| PDB | drug | warhead err | centroid | overlap | RMSD | energy (Ki) | drift |
|---|---|---|---|---|---|---|---|
| 6OIM | sotorasib | 1.60 Å | 2.46 Å | 68 % | 5.67 Å | −8.50 (586 nM) | 0.007 Å |
| 7RFS | nirmatrelvir | 0.77 Å | 1.61 Å | 86 % | 2.28 Å | −12.67 (519 pM) | 0.002 Å |
| 4G5J | afatinib | 0.93 Å | 0.56 Å | **100 %** | **1.12 Å** | −8.38 (720 nM) | 0.008 Å |
| 5P9J | ibrutinib | 0.16 Å | 0.58 Å | 94 % | **1.06 Å** | −10.48 (20.8 nM) | 0.013 Å |

These are single runs of a 10-run GA search seeded with `seed pid time`, so
expect ~1 kcal/mol and some pose variation between repeats (6OIM has come in
anywhere from -8.3 to -9.7 kcal/mol across runs). Set a fixed `seed` in
`_build_cov_dpf` if a reproducible benchmark is ever wanted.

Reading it honestly against the usual ≤2 Å "correct pose" bar: **7RFS and 5P9J
are solved** (1.25 and 1.19 Å), **6OIM is marginal** (2.21 Å) though its
centroid and 95 % overlap say the placement is right, and **4G5J is only
partially right** — see below. Three of four is a good result for a docking
method, not a perfect one.

**Anchor drift is 0.002–0.013 Å in every case**, which is the point of the
protocol: AutoDock4 pins the flexible residue's root, so the covalent
attachment holds no matter what the search does. Compare the Vina attempt
below, where it never held.

Energies land in a sensible range and rank plausibly — nirmatrelvir, the
picomolar clinical Mpro inhibitor, scores strongest at −12.76 kcal/mol. These
are AD4 estimates, not measurements, and no covalent-bond energetics are in
them (see caveats), so do not read the Ki values as predicted potencies.

**On the earlier protein-only numbers: don't treat them as a baseline.** An
earlier version of this engine stripped every HETATM, and its figures (6OIM
0.35 Å centroid / 95 % overlap / 2.21 Å RMSD; 4G5J 7.33 Å RMSD / 76 %
overlap) appear in this file's history. A receptor with its cofactors deleted
is a less physical model of the binding site, so where those numbers look
better that is not evidence that deleting cofactors helps — it is a comparison
against an artifact. The table above is the result.

What can fairly be said is that 4G5J no longer disagrees with itself. Under
protein-only it showed the signature of a mis-mapped pose — a low centroid
error alongside a 7.33 Å RMSD — and that signature is gone (0.56 Å centroid,
1.12 Å RMSD, 100 % overlap). Worth knowing when reading that row: 4G5J also
holds two copies of *unreacted* afatinib (`0WM A1102`, `0WM A1103`) at
13.65 Å and 20.09 Å from Cys797 SG, both retained; the nearer one is probably
inside the grid box and may be bounding the search.

All single runs. Four same-receptor repeats of 6OIM spanned 1.40 kcal/mol
(−8.34 to −9.74), so one-run differences of that order between rows are noise,
and pose-metric differences have not been measured against repeats at all.

## Agent-loop test

Exercised end-to-end through `dock_assist.py` (piped stdin, `--print`), not
just via direct Python calls, on both the free default model `gemma4:31b` and
`glm-5.3`. Both chose the intended chain unaided and with correct arguments —
no hallucinated kwargs, which is the module-global design holding up:

```
get_protein_from_pdb -> smiles_node -> find_ligand_centroids
  -> detect_nucleophiles_agent(centroid)  -> "Closest dockable residue: CYS A12"
  -> dock_covalent_agent(chain='A', resseq=12)
```

glm-5.3 finished at -9.74 kcal/mol (Ki 72 nM), gemma4:31b at -8.34 kcal/mol,
both with 0.007 A anchor drift. `detect_nucleophiles_agent` picked Cys A12 out
of eight nearby nucleophiles, and gemma's own reasoning noted it "matches the
G12C mutation".

Two defects showed up here that `test_covalent.py` could not, because it feeds
pre-cleaned receptors and bypasses the dispatch layer: the HETATM handling
above, and tool exceptions killing the chat session (fixed separately in
`chat_turn`; the covalent wrappers also return `DockError` as text).

One pre-existing fragility worth knowing, not caused by this work: on one
`gemma4:31b` run the model passed `'6OIM.pdb'` to `find_ligand_centroids`
instead of the `pdb_files/...` path `get_protein_from_pdb` had just reported.
Stronger models thread it correctly, so this mainly affects whoever runs the
free default — which is everyone who clones the repo.

## The Vina dead end (recorded so it is not retried)

Before AutoDock4 was installed, an engine was built on Vina 1.2.7 with the
tether enforced geometrically instead of by the docking engine: a conformer
pool anchored by 3-point superposition, a clash screen, and `--score_only`
ranking. **It does not work**, and the measurements say why:

| Attempt | Result |
|---|---|
| top 30 conformers by **most contacts** within 4.5 Å | all positive, best +9.81 kcal/mol; pose centroid **20.5 Å** from the crystal ligand, **0 %** overlap |
| all 181 clash-free, `--local_only` refined | 154/181 negative, best −3.69 — but every top pose had **2.7–8.5 Å** anchor drift |
| all 181 clash-free, `--score_only` anchored | all +9.05…+9.46 (0/181 negative); ranking usable (closest-to-crystal pose 5th of 181) but the pose still 17.5 Å out on centroid |

Three transferable lessons:

1. **`--local_only` cannot hold a tether.** It is an unrestrained BFGS
   relaxation, free to translate and rotate the whole ligand, and an anchored
   conformer is strained — so Vina slides the adduct off its anchor to buy
   score. Vina has no restraint option, so a Vina tether can only be held by
   freezing the pose, which rules out any search.
2. **"Most buried is best" is actively harmful** as a conformer-selection
   heuristic: it picks the most deeply jammed poses. Ranking belongs to the
   scoring function; a cheap screen should only discard the impossible.
3. **A filtered random conformer ensemble is not a search.** For a drug-sized
   ligand anchored in an enclosed site it almost never contains the right
   pose: the clash-free fraction was 181/705 for 6OIM, 9/809 for 4G5J, and
   **0/840** for 7RFS, where the enclosed Mpro S1′ site admitted nothing at
   all. AutoDock4 solves the same 7RFS case because it *searches* the
   torsions rather than sampling them blindly.

The general lesson, and the reason this section exists: **the covalent tether
has to be enforced by the docking engine.** Vina cannot do it; AutoDock4 can,
and it is one small free binary away.

## Honest scope and caveats

- **Cys and Ser only.** `COVALENT_RESIDUES` holds the residues whose
  CA–CB–nucleophile atoms form a contiguous 3-atom chain. Lys/Thr/Tyr/His are
  *detected and reported* by `detect_nucleophiles_agent` but not dockable.
- **No covalent-bond energetics.** Warhead intrinsic reactivity, the reaction
  barrier and reversibility are not modelled. Two ligands with the same pose
  and different warheads score identically, and a reversible nitrile is
  treated like an irreversible acrylamide.
- **Rigid receptor apart from the removed side chain.** No other side chain
  moves, so an inhibitor needing induced fit will score badly.
- **The mimic must be unambiguous.** If the ligand itself contains an
  ethylthio (Cys) or ethoxy (Ser) group the mimic SMARTS matches more than
  once and the run fails with an explicit error rather than tethering the
  wrong atom. None of the four test drugs trips this.
- **Scores are not comparable to the Vina-based tools**, only to each other.
- **Minutes per ligand**, since each ligand gets fresh AutoGrid maps and a
  full GA search. The rigid receptor PDBQT is cached; the maps are not,
  because the box depends on the ligand's reach.

## Next steps

**Nothing in this section is required.** The engine is complete and validated
(see [Results](#redocking-results)); these are optional improvements, ordered
by what would actually be worth doing. The earlier draft of this file listed
reactive docking as item 1, which was a stale ranking from when the Vina
engine was failing and reactive docking was the proposed *fix*; it is now just
an optional extra capability and has moved down accordingly. The previous
item 1 -- explaining 4G5J's outlying RMSD -- is gone: it was caused by
stripping heteroatoms from the receptor, and is fixed.

### 1. Label reversible warheads

Nitriles, aldehydes and boronates are reversible; the engine scores them as if
irreversible, so nirmatrelvir's -12.76 kcal/mol reads like a potency it does
not imply. Even labelling them in the report and declining to rank them
against irreversible warheads would stop a misleading comparison. Cheap, and
it removes a way for the tool to mislead.

### 2. More residues

Lys, Tyr, His and Thr need a longer or cyclic side-chain mimic. Much easier
under AutoDock4 than it would have been on the Vina engine, because AD4
searches the mimic's own torsions too: the mimic only has to reach from CA/CB
to the nucleophile, and the tutorial's own His example
(`tether_smarts "n1cc(CC)nc1"`) shows the pattern for a cyclic one. His is the
one worth doing next.

### 3. Cache the maps across a series

Each ligand gets its own AutoGrid run because the box is sized to that
ligand's reach. For a congeneric series the boxes are nearly identical, so
sizing one box over the whole series' envelope and reusing the maps would cut
most of the runtime (currently a few minutes per ligand, ~35 MB of maps each).

### 4. Meeko reactive docking (the old "option B")

An **additional** engine, not a missing piece. meeko 0.8.0 already exposes it:
`-r/--reactive_flexres`, `--reactive_name`, `-s/--reactive_name_specific`,
`--r_eq_12`, `--eps_12`.

Instead of tethering, the receptor's reactive atom gets a **modified 12-6 pair
potential** that pulls the warhead atom in, and the pose is scored off
AutoGrid maps. That is structurally identical to the AD4Zn / tier metal
machinery already here: `_build_ad4map_gpf`, `_run_autogrid` and the
`NBP_DONOR_BASE` / `METAL_SITE_R` well-depth idea all transfer, with the metal
donor well replaced by a warhead well, and the hard-won **`eps/6` AutoGrid
semantics measured for the metal path apply unchanged**.

What it would buy: the warhead finds *which* nucleophile it attacks rather
than being told, which suits screening a site that has several candidate
cysteines, or a target where the reactive residue is unknown. The tethered
engine already covers the normal case where you know the residue.

**Verify first:** the protocol was benchmarked with AutoDock-GPU and
prescribes a reactive pair-damping scheme the engine must apply. Now that
AutoDock4 is installed, check whether AD4 honours the reactive types in the
maps — run one redock from the test set (7RFS or 6OIM) and confirm the
reactive atom types survive into the maps before writing a tool.

### 5. More warheads

beta-lactam, boronic acid/boronate (reversible, and forms a tetrahedral
boronate, so it needs a geometry change at B rather than just a bond),
activated esters/carbamates, sulfonyl fluoride, alpha-ketoamide. Each is one
entry in `COVALENT_WARHEADS`.

### 6. Harden the RMSD metric

`pose_metrics` reports symmetry-corrected RMSD only when
`AssignBondOrdersFromTemplate` succeeds against the crystal HET read with
proximity bonding, and `n/a` otherwise; the robust metrics (warhead error,
centroid error, overlap) always compute. Options: read the ligand's ideal SDF
from RCSB as the bond-order template, or use `obrms` (present in the
`adfr-suite` env) on a mimic-stripped pose. Lower priority now that
retaining heteroatoms resolved the one case where the metric disagreed with
the other measures (4G5J).

### 7. Structural waters

Waters are dropped from the receptor (cofactors and ions are kept). A water
that genuinely bridges ligand and protein in the site is therefore not
represented. AutoDock's hydrated-docking protocol is the established answer if
this ever matters for a target of interest.
