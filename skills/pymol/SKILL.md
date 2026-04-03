---
name: pymol
description: Extract a QM cluster model from a protein-ligand PDB file. Guides ligand identification, residue selection, amide capping, and charge determination. Use when preparing a cluster for QM or QM/MM calculations.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# Cluster Model Extraction

Write a Python script tailored to the user's PDB, assembled from the building blocks below. Inspect the structure first, ask questions, then compose only what is needed.

## Pre-flight Checklist

Read the PDB and fill in this checklist. For items you can resolve automatically, fill in the answer. For ambiguous items, ask the user.

**Important:** `byres (system within X of ligand)` naturally selects from all chains. Do not filter by chain — multiple chains means a dimer interface, include them all.

```
═══════════════════════════════════════════════════════════
  CLUSTER EXTRACTION — PRE-FLIGHT CHECKLIST
═══════════════════════════════════════════════════════════

  PDB file:     _______________
  PDB format:   [ ] Standard    [ ] CHARMM (needs normalization)

─── Ligand ────────────────────────────────────────────────

  Residue name: _______________
  Description:  _______________
  Charge:       ___   (reasoning: _________________________)

─── Cutoff ────────────────────────────────────────────────

  Distance:     [ ] 3.5 Å (ligand + immediate contacts)
                [ ] 4.5 Å (compact active site — recommended)
                [ ] 5.0 Å (standard active site)
                [ ] ___ Å (custom)

─── Protonation states ────────────────────────────────────

  NOTE: For CHARMM-GUI PDBs, residue names may NOT be updated
  when protonation states are changed. A protonated ASP may
  still be named "ASP" but carry an HD2 atom. Always verify
  by hydrogen counting (see Block below), not residue names.

  Histidines near ligand:
    Resi ____  Chain ____  [ ] HID/HSD (Nδ-H, 0)
                           [ ] HIE/HSE (Nε-H, 0)
                           [ ] HIP/HSP (+1)

  Non-standard protonation (if any):
    Resi ____  Chain ____  Name: ____  Charge: ___
    (ASH/ASPP = protonated ASP → 0;  GLH/GLUP = protonated GLU → 0;
     LYN = deprotonated LYS → 0)

─── Cysteines ─────────────────────────────────────────────

  (SG–SG < 2.5 Å = disulfide, neutral;
   SG–metal < 3.0 Å = CYM, −1)

  Resi ____  Chain ____  [ ] Neutral  [ ] CYM (−1)

─── Metals ────────────────────────────────────────────────

  [ ] None
  [ ] Present: ____  Resi: ____  Oxidation: ____
      Coordinating: _______________

─── Spin ──────────────────────────────────────────────────

  [ ] Singlet (1)
  [ ] Open-shell: multiplicity = ___

═══════════════════════════════════════════════════════════
```

### Key guidance

**Ligand identification** — Do NOT rely on HETATM records. CHARMM PDBs write everything as ATOM. Find residue names that are not standard amino acids (ALA–VAL), not HIS variants (HSD/HSE/HSP/HID/HIE/HIP), not caps (ACE/NME/NMA), and not solvent (HOH/WAT/TIP3/SOL). What remains is the ligand.

**CHARMM-GUI protonation** — CHARMM-GUI may add or remove protons from titratable residues without renaming them. A protonated Asp may still be called "ASP" but carry HD1 or HD2. A deprotonated Lys may still be called "LYS" with only 2 HZ atoms. Do NOT rely on residue names alone. Always run the "Verify protonation by hydrogen counting" block and use those results for all charge calculations.

**Cysteines** — Measure SG–SG distances with `cmd.get_distance`. Also check SG–metal distances for metal coordination.

**Metals** — Search for Zn, Fe, Cu, Mg, Mn, Co, Ni. A system may have multiple metal sites — only the ones that fall within the cutoff of the ligand belong in the cluster. For metals that ARE in the cluster, also force-include their coordinating residues (typically HIS, CYS, MET — check metal–N/S distances < 3.0 Å) even if those residues are outside the cutoff. Do NOT force-include metals that are far from the ligand.

---

## Runtime

**NEVER use `python3` or `python` to run these scripts.** Always run via PyMOL:

```bash
pymol -cq -r your_script.py
```

`-c` = headless (no GUI), `-q` = quiet (no banner). PyMOL handles startup and shutdown — scripts do not need `pymol.finish_launching` or `cmd.quit`.

---

## Code Building Blocks

### Block: Boilerplate

```python
from pymol import cmd
import numpy as np
```

### Block: Load structure

```python
cmd.load('FILL_PDB_PATH', 'system')
```

### Block: Normalize CHARMM PDB

Remaps segment IDs to unique chain letters. Safe to run on any PDB — no-op if not needed.

```python
chain_segs = {}
cmd.iterate('system and name CA',
            'chain_segs.setdefault(chain, set()).add(segi)',
            space={'chain_segs': chain_segs})

if any(len(segs) > 1 for segs in chain_segs.values()):
    segi_order = []
    _seen = set()
    cmd.iterate('system', 'lst.append(segi) if segi not in _seen else None; _seen.add(segi)',
                space={'lst': segi_order, '_seen': _seen})

    chain_pool = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    segi_to_chain = {s: chain_pool[i] for i, s in enumerate(segi_order)}
    cmd.alter('system', 'chain = segi_to_chain[segi]', space={'segi_to_chain': segi_to_chain})
    cmd.sort()
    for segi, ch in segi_to_chain.items():
        print(f"  {segi} → {ch}")
```

### Block: Fix element column

Fixes blank/wrong element columns in CHARMM PDBs. Handles 2-letter elements (Cl, Br, Fe, etc.) and strips lone pairs/Drude particles.

```python
TWO_LETTER_ELEMS = {'CL', 'BR', 'FE', 'ZN', 'MG', 'NA', 'CA', 'CU', 'MN',
                    'LI', 'NI', 'CO', 'SE', 'MO', 'CR'}

PROTEIN_RESNAMES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP', 'ASH', 'GLH', 'LYN',
    'ACE', 'NME', 'NMA',
}

def _infer_element(atom_name, resn):
    clean = atom_name.strip().upper()
    if resn.strip().upper() in PROTEIN_RESNAMES:
        for ch in clean:
            if ch.isalpha(): return ch
    for elem in TWO_LETTER_ELEMS:
        if clean.startswith(elem):
            return elem[0] + elem[1].lower()
    for ch in clean:
        if ch.isalpha(): return ch
    return 'X'

# Strip lone pairs, Drude particles, massless sites
lp_atoms = []
cmd.iterate('system', 'lp_atoms.append(index) if name.strip().upper().startswith(("LP","DRUDE","MW")) else None',
            space={'lp_atoms': lp_atoms})
if lp_atoms:
    cmd.remove('index ' + '+'.join(str(i) for i in lp_atoms))
    print(f"Removed {len(lp_atoms)} non-physical atoms")

# Fix elements
atom_data = []
cmd.iterate('system', 'atom_data.append((index, name, resn, elem))',
            space={'atom_data': atom_data})
fixes = {}
for idx, name, resn, current_elem in atom_data:
    correct = _infer_element(name, resn)
    if current_elem.strip().upper() != correct.strip().upper():
        fixes[idx] = correct
if fixes:
    cmd.alter('system', 'elem = fixes.get(index, elem)', space={'fixes': fixes})
    print(f"Fixed element column for {len(fixes)} atoms")
```

### Block: Select inner region

```python
# FILL
LIGAND_RESN = 'LIG'
CUTOFF = 4.5

cmd.select('ligand', f'system and resn {LIGAND_RESN}')
cmd.select('inner', f'byres (system within {CUTOFF} of ligand) or ligand')
```

### Block: Include coordinating residues for metals in the cluster

Only metals already inside the cutoff selection get their coordination shells
added. Metals far from the ligand are left out entirely.

```python
# Find metals that landed in the inner selection
METAL_ELEMS = {'ZN', 'FE', 'CU', 'MG', 'MN', 'NI', 'CO', 'MO', 'CR'}
metals_in_inner = []
cmd.iterate('inner',
            'metals_in_inner.append((index, chain, resi, resn, elem)) '
            'if elem.strip().upper() in METAL_ELEMS else None',
            space={'metals_in_inner': metals_in_inner, 'METAL_ELEMS': METAL_ELEMS})

if metals_in_inner:
    print(f"Metals in cluster: {len(metals_in_inner)}")
    for idx, chain, resi, resn, elem in metals_in_inner:
        print(f"  {elem} ({resn}) chain {chain} resi {resi}")
        # Add residues coordinating this metal (N, O, S within 3.0 Å)
        cmd.select('_coord',
                    f'byres ((elem N+O+S) and system within 3.0 of (index {idx}))')
        cmd.select('inner', 'inner or _coord')
    cmd.delete('_coord')
    print(f"After adding coordination shells: {cmd.count_atoms('inner')} atoms")
else:
    # Check if there are metals in the system but outside the cluster
    all_metals = []
    cmd.iterate('system',
                'all_metals.append(elem) if elem.strip().upper() in METAL_ELEMS else None',
                space={'all_metals': all_metals, 'METAL_ELEMS': METAL_ELEMS})
    if all_metals:
        print(f"Note: {len(all_metals)} metal(s) in system but none within {CUTOFF} Å of ligand — excluded from cluster")
```

### Block: Collect residue info

```python
cluster_res = set()
cmd.iterate(f'inner and not resn {LIGAND_RESN} and name CA',
            'cluster_res.add((chain, int(resi)))',
            space={'cluster_res': cluster_res})
all_res = set()
cmd.iterate(f'system and not resn {LIGAND_RESN} and name CA',
            'all_res.add((chain, int(resi)))',
            space={'all_res': all_res})
print(f"Cluster residues: {len(cluster_res)}")
```

### Pre-flight: Protonation diagnostic

**Run this standalone script before writing the extraction script.** It collects every
hydrogen atom on every titratable residue across the full system and prints them as a
table. Read the output, apply your knowledge of each amino acid's standard hydrogen
complement, and determine the actual protonation state of each residue. Record your
conclusions as `actual_charges` in the extraction script — do not let the script infer
charges automatically, as hardcoded rules cannot cover every naming convention.

#### PyMOL coding notes for this script

**Collecting atom names** — use `cmd.iterate` with a list captured via the `space` dict.
Any variable referenced inside the expression string must appear in `space`, otherwise
PyMOL uses an empty namespace and the result is silently empty:

```python
names = []
cmd.iterate('system and chain A and resi 25',
            'names.append(name.strip())',
            space={'names': names})
```

**Always call `.strip()`** — CHARMM PDB atom names are stored with leading/trailing
spaces (e.g., `' HD2'`). Without `.strip()`, the string `'HD2'` will not match.

**Selection scoping** — always prefix selections with the object name (`system`, `inner`)
to avoid cross-contamination when multiple objects are loaded simultaneously:

```python
cmd.iterate('system and chain A and resi 25', ...)  # safe
cmd.iterate('chain A and resi 25', ...)             # risky
```

#### Diagnostic script

```python
# Run with: pymol -cq -r hcount.py
from pymol import cmd

TITRATABLE = {'ASP', 'GLU', 'LYS', 'ARG',
              'HIS', 'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP',
              'CYS', 'CYM', 'ASH', 'ASPP', 'GLH', 'GLUP', 'LYN'}

cmd.load('FILL_PDB_PATH', 'system')

tit_res = []
cmd.iterate('system and name CA',
            'tit_res.append((chain, resi, resn)) if resn in TITRATABLE else None',
            space={'tit_res': tit_res, 'TITRATABLE': TITRATABLE})

print(f"{'Resn':<6} {'Ch':>2} {'Resi':>4}  {'nH':>3}  H atoms")
print("-" * 70)
for chain, resi, resn in sorted(tit_res, key=lambda x: (x[0], int(x[1]))):
    names = []
    cmd.iterate(f'system and chain {chain} and resi {resi}',
                'names.append(name.strip())', space={'names': names})
    h_names = sorted(n for n in names if n.startswith('H'))
    print(f"{resn:<6} {chain:>2} {resi:>4}  {len(h_names):>3}  {h_names}")
```

#### Interpreting the output

For each row, compare the H atom list against the expected hydrogen complement for that
residue type. Any deviation from the standard count indicates a non-standard protonation
state — regardless of what the residue name says.

Use your built-in knowledge of amino acid chemistry to make this determination. As a
reference, the sidechain H atoms that change with protonation state are:

| Residue | Extra H when protonated | Missing H when deprotonated |
|---------|------------------------|------------------------------|
| ASP (std −1) | HD1 or HD2 on carboxylate → charge 0 | — |
| GLU (std −1) | HE1 or HE2 on carboxylate → charge 0 | — |
| LYS (std +1) | — | fewer than 3 HZ on NZ → charge 0 |
| ARG (std +1) | — | rarely deprotonated; flag if HH/HE missing |
| HIS (std 0) | HD1 = Nδ protonated (HSD); HE2 = Nε protonated (HSE); both = HIP (+1) | |
| CYS (std 0) | — | HG absent on SG → −1 or disulfide (check SG–SG < 2.5 Å) |

#### After running the diagnostic

Read the printed table row by row. For each titratable residue, compare its H atom list
against the expected complement for that residue type (reference table above) and assign
a charge. Then write your conclusions to `protonation.json` using the Write tool — do
not keep them only in conversation context. The extraction script will load this file.

`protonation.json` schema:

```json
{
  "ligand": {"resn": "LIG", "charge": 0},
  "residues": [
    {"chain": "B", "resi": 25, "resn": "ASP", "charge":  0, "h_atoms": ["HA","HB1","HB2","HD2","HN"], "note": "protonated — has HD2"},
    {"chain": "B", "resi": 29, "resn": "ASP", "charge": -1, "h_atoms": ["HA","HB1","HB2","HN"],       "note": "deprotonated"},
    {"chain": "B", "resi": 45, "resn": "LYS", "charge": +1, "h_atoms": ["HA","HB1","HB2","HG","HZ1","HZ2","HZ3","HN"], "note": "3 HZ — charged"}
  ],
  "inner_charge": -2,
  "system_charge": 7
}
```

Include every titratable residue in the system (not only the inner region) so total
system charge can also be verified. Include `h_atoms` verbatim from the diagnostic
output so the file is self-documenting and auditable.

### Block: Find backbone cut bonds

```python
cuts = []
for (chain, resi) in sorted(cluster_res):
    if (chain, resi + 1) in all_res and (chain, resi + 1) not in cluster_res:
        cuts.append(('C_terminal', chain, resi, resi + 1))
    if (chain, resi - 1) in all_res and (chain, resi - 1) not in cluster_res:
        cuts.append(('N_terminal', chain, resi, resi - 1))
print(f"Backbone cuts: {len(cuts)}")
```

### Block: Place capping hydrogens

One H per cut bond, along the original bond vector.

```python
def _get_xyz(sel, atom_name):
    coords = []
    cmd.iterate_state(1, f'({sel}) and name {atom_name}',
                      'coords.append((x, y, z))', space={'coords': coords})
    if len(coords) != 1:
        raise ValueError(f"Expected 1 atom for '{sel} name {atom_name}', got {len(coords)}")
    return np.array(coords[0])

CH_BOND = 1.09  # C–H
NH_BOND = 1.01  # N–H

cap_atoms = []
for cut_type, chain, inside_resi, outside_resi in cuts:
    base = f'chain {chain} and resi'
    if cut_type == 'C_terminal':
        inside_xyz  = _get_xyz(f'{base} {inside_resi}',  'C')
        outside_xyz = _get_xyz(f'{base} {outside_resi}', 'N')
        bond_len = CH_BOND
    else:
        inside_xyz  = _get_xyz(f'{base} {inside_resi}',  'N')
        outside_xyz = _get_xyz(f'{base} {outside_resi}', 'C')
        bond_len = NH_BOND

    vec = outside_xyz - inside_xyz
    h_xyz = inside_xyz + bond_len * vec / np.linalg.norm(vec)

    resn_buf = []
    cmd.iterate(f'{base} {inside_resi} and name CA', 'resn_buf.append(resn)',
                space={'resn_buf': resn_buf})
    cap_atoms.append(('HC', resn_buf[0], chain, inside_resi, h_xyz))
print(f"Capping H atoms: {len(cap_atoms)}")
```

### Block: Save cluster PDB with caps

```python
OUTPUT_PDB = 'cluster.pdb'  # FILL

cmd.save('_cluster_raw.pdb', 'inner')
with open('_cluster_raw.pdb') as f:
    lines = [l for l in f if not l.startswith('END')]

serials = [int(l[6:11]) for l in lines if l.startswith(('ATOM', 'HETATM'))]
next_serial = (max(serials) + 1) if serials else 1

with open(OUTPUT_PDB, 'w') as f:
    f.writelines(lines)
    for atom_name, resn, chain, resi, (x, y, z) in cap_atoms:
        f.write(
            f"ATOM  {next_serial:5d}  {atom_name:<3s} {resn:3s} "
            f"{chain:1s}{int(resi):4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00           H\n"
        )
        next_serial += 1
    f.write('END\n')

import os; os.remove('_cluster_raw.pdb')
print(f"Wrote {OUTPUT_PDB}")
```

### Block: Calculate and report cluster charge

Loads protonation states from `protonation.json` written during the pre-flight
diagnostic. Only residues that appear in the inner region contribute to the cluster
charge; the JSON may contain the full system so this filters by `(chain, resi)`.

```python
import json

with open('protonation.json') as f:
    prot = json.load(f)

LIGAND_CHARGE = prot['ligand']['charge']

# Build lookup from JSON; only count residues present in the inner region
all_prot = {(r['chain'], str(r['resi'])): r for r in prot['residues']}

inner_res_data = []
cmd.iterate(f'inner and not resn {LIGAND_RESN} and name CA',
            'inner_res_data.append((chain, resi, resn))',
            space={'inner_res_data': inner_res_data})

actual_charges = {}
missing = []
for chain, resi, resn in inner_res_data:
    key = (chain, str(resi))
    if key in all_prot:
        actual_charges[key] = all_prot[key]['charge']
    # Non-titratable residues (ALA, GLY, etc.) carry no charge — skip silently

protein_charge = sum(actual_charges.values())
total_charge = protein_charge + LIGAND_CHARGE

print(f"\nCluster charge breakdown:")
for (ch, ri), q in sorted(actual_charges.items(), key=lambda x: (x[0][0], int(x[0][1]))):
    resn = all_prot[(ch, ri)]['resn']
    note = all_prot[(ch, ri)].get('note', '')
    print(f"  {resn} chain {ch} resi {ri}: {q:+d}  ({note})")
print(f"  Ligand ({prot['ligand']['resn']}): {LIGAND_CHARGE:+d}")
print(f"  ─────────────────")
print(f"  Inner total: {total_charge:+d}")
```

---

## Notes
- `byres` expands distance selections to complete residues.
- Atom names are case-sensitive: `name CA` not `name ca`.
- If the PDB already has hydrogens, skip capping and save directly.
- Always present the per-residue protonation audit table to the user before reporting the final charge.
