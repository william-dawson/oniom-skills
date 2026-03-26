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

**Cysteines** — Measure SG–SG distances with `cmd.get_distance`. Also check SG–metal distances for metal coordination.

**Metals** — Search for Zn, Fe, Cu, Mg, Mn, Co, Ni. Force-include coordinating residues regardless of cutoff.

---

## Runtime

**NEVER use `python3` or `python` to run these scripts.** The `pymol` module is only available inside PyMOL's bundled Python. Always use:

```bash
PYMOL_PYTHON=$(head -1 "$(which pymol)" | sed 's/#!//')
"$PYMOL_PYTHON" your_script.py
```

This is the only Python interpreter that will work.

---

## Code Building Blocks

### Block: Boilerplate

```python
import pymol
from pymol import cmd
import numpy as np

pymol.finish_launching(['pymol', '-c'])
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

# Force-include metal-coordinating residues if needed:
# cmd.select('inner', 'inner or (system and resn ZN)')
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

Two methods as a sanity check: (1) residue-name table, (2) counting +/− atom name annotations (CHARMM PDBs encode charges as `N1+`, `O1-` in atom names).

```python
LIGAND_CHARGE = 0  # FILL

RESIDUE_CHARGES = {
    'ASP': -1, 'GLU': -1,
    'ARG': +1, 'LYS': +1,
    'HIS':  0, 'HID':  0, 'HIE':  0, 'HIP': +1,
    'HSD':  0, 'HSE':  0, 'HSP': +1,
    'ASH':  0, 'ASPP': 0, 'GLH':  0, 'GLUP': 0, 'LYN': 0,
    'CYM': -1,
}

# Method 1: residue-name table
resnames = []
cmd.iterate(f'inner and not resn {LIGAND_RESN} and name CA',
            'resnames.append(resn)', space={'resnames': resnames})
protein_charge = sum(RESIDUE_CHARGES.get(r, 0) for r in resnames)
total_m1 = protein_charge + LIGAND_CHARGE

# Method 2: count +/- in atom names (CHARMM annotation)
charge_data = {'pos': 0, 'neg': 0}
cmd.iterate('inner',
            'charge_data["pos"] += (1 if "+" in name.strip() else 0); '
            'charge_data["neg"] += (1 if "-" in name.strip() else 0)',
            space={'charge_data': charge_data})
total_m2 = charge_data['pos'] - charge_data['neg']

print(f"\nCharge (residue names): {total_m1:+d}")
print(f"Charge (atom annotations, excl. ligand): {total_m2:+d}")
if total_m1 != total_m2:
    print(f"*** WARNING: methods disagree — check protonation states ***")
```

### Block: Teardown

```python
cmd.quit()
```

---

## Notes
- Always call `cmd.quit()` or the process hangs.
- `byres` expands distance selections to complete residues.
- Atom names are case-sensitive: `name CA` not `name ca`.
- If the PDB already has hydrogens, skip capping and save directly.
