---
name: pymol
description: Extract a QM cluster model from a protein-ligand PDB file. Guides ligand identification, residue selection, amide capping, and charge determination. Use when preparing a cluster for QM or QM/MM calculations.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# Cluster Model Extraction

Your job is to write a Python script tailored to the user's specific PDB, assembled from the building blocks below. Do not use a one-size-fits-all script — inspect the structure first, ask questions, then compose only what is needed.

## Before Writing Any Code

Read the PDB file and reason through the following. Engage the user on anything ambiguous.

**Detect CHARMM PDB format**
CHARMM-GUI PDBs use segment IDs (columns 73–76, e.g. `PROA`, `PROB`, `HETA`) to distinguish chains, but assign them the same chain letter (e.g. both protein chains are `P`). This causes residue number collisions. Check for this and, if detected, use the "Normalize CHARMM PDB" block immediately after loading.

**Identify the ligand**
Look at HETATM records (excluding HOH/WAT). If more than one non-water HETATM group is present, ask the user which is the ligand of interest.

**Assess the ligand charge**
Look at the ligand's element composition and any CONECT records. Reason about functional groups:
- Carboxylate → likely −1
- Ammonium/guanidinium → likely +1
- Phosphate → likely −2
- If ambiguous, ask: *"The ligand [RESNAME] looks like [description]. What is its formal charge?"*

**Check histidines**
Note any HIS residues near the ligand. If the PDB uses generic `HIS` rather than `HID`/`HIE`/`HIP`, the protonation state is unknown. Ask the user or note it as an assumption.

**Check for metals**
If a metal ion is present (Zn, Fe, Cu, Mg, Mn, etc.), identify its coordinating residues. These must be in the cluster regardless of distance cutoff. Flag the metal's likely oxidation state and ask the user to confirm — this affects both charge and spin multiplicity.

**Check for open-shell character**
If a metal with unpaired d-electrons or a radical intermediate is present, flag it explicitly and ask for the spin multiplicity before proceeding. Do not assume singlet.

---

## Runtime

PyMOL ships its own Python. Always discover it at runtime:

```bash
PYMOL_PYTHON=$(head -1 "$(which pymol)" | sed 's/#!//')
"$PYMOL_PYTHON" your_script.py
```

---

## Code Building Blocks

Assemble these into a single script. Each block has clearly marked inputs you must fill in.

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

### Block: Normalize CHARMM PDB (use when CHARMM format detected)

CHARMM PDBs put multiple segments on the same chain letter, so `chain` alone
cannot distinguish them. This block detects the collision and remaps segment IDs
to unique chain letters, so all downstream code using `chain` works correctly.

```python
# Check if any chain letter has multiple segment IDs
chain_segs = {}
cmd.iterate('system and name CA',
            'chain_segs.setdefault(chain, set()).add(segi)',
            space={'chain_segs': chain_segs})

if any(len(segs) > 1 for segs in chain_segs.values()):
    # Collect unique segment IDs in order
    segi_order = []
    _seen = set()
    cmd.iterate('system', 'lst.append(segi) if segi not in _seen else None; _seen.add(segi)',
                space={'lst': segi_order, '_seen': _seen})

    # Map each segment ID to a unique chain letter
    chain_pool = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
    segi_to_chain = {s: chain_pool[i] for i, s in enumerate(segi_order)}
    cmd.alter('system', 'chain = segi_to_chain[segi]', space={'segi_to_chain': segi_to_chain})
    cmd.sort()  # re-sort after chain reassignment

    print("CHARMM PDB detected — remapped segment IDs to chain letters:")
    for segi, ch in segi_to_chain.items():
        print(f"  {segi} → {ch}")
```

### Block: Select inner region

```python
# FILL: ligand residue name, cutoff distance
LIGAND_RESN = 'LIG'
CUTOFF = 5.0

cmd.select('ligand', f'system and resn {LIGAND_RESN}')
cmd.select('inner', f'byres (system within {CUTOFF} of ligand) or ligand')

# If a metal must be included regardless of distance, add:
# cmd.select('inner', 'inner or (system and resn ZN)')
```

### Block: Collect residue info

```python
cluster_res = set()
cmd.iterate(
    f'inner and not resn {LIGAND_RESN} and name CA',
    'cluster_res.add((chain, int(resi)))',
    space={'cluster_res': cluster_res}
)
all_res = set()
cmd.iterate(
    f'system and not resn {LIGAND_RESN} and name CA',
    'all_res.add((chain, int(resi)))',
    space={'all_res': all_res}
)
print(f"Cluster residues: {sorted(cluster_res)}")
```

### Block: Find backbone cut bonds

```python
# Each entry: (cut_type, chain, inside_resi, outside_resi)
cuts = []
for (chain, resi) in sorted(cluster_res):
    # C-terminal boundary: C(resi)–N(resi+1), outside is resi+1
    if (chain, resi + 1) in all_res and (chain, resi + 1) not in cluster_res:
        cuts.append(('C_terminal', chain, resi, resi + 1))
    # N-terminal boundary: C(resi-1)–N(resi), outside is resi-1
    if (chain, resi - 1) in all_res and (chain, resi - 1) not in cluster_res:
        cuts.append(('N_terminal', chain, resi, resi - 1))

print(f"Backbone cuts: {len(cuts)}")
for c in cuts:
    print(f"  {c}")
```

### Block: Place capping hydrogens

One H per cut bond, placed along the original bond vector. No stub atoms.

```python
def _get_xyz(sel, atom_name):
    coords = []
    cmd.iterate_state(1, f'({sel}) and name {atom_name}',
                      'coords.append((x, y, z))', space={'coords': coords})
    if len(coords) != 1:
        raise ValueError(f"Expected 1 atom for '{sel} name {atom_name}', got {len(coords)}")
    return np.array(coords[0])

CH_BOND = 1.09  # Å, C–H
NH_BOND = 1.01  # Å, N–H

cap_atoms = []  # list of (atom_name, resn, chain, resi, xyz)
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

print(f"Capping H atoms placed: {len(cap_atoms)}")
```

### Block: Save cluster PDB with caps

```python
# FILL: output path
OUTPUT_PDB = 'cluster.pdb'

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

```python
# FILL: LIGAND_CHARGE (from user), and override any residue charges
# that differ from the defaults below (e.g. a protonated ASP)
LIGAND_CHARGE = 0  # FILL

RESIDUE_CHARGES = {
    'ASP': -1, 'GLU': -1,
    'ARG': +1, 'LYS': +1,
    'HIS':  0, 'HID':  0, 'HIE':  0, 'HIP': +1,
    'CYM': -1,
}

resnames = []
cmd.iterate(
    f'inner and not resn {LIGAND_RESN} and name CA',
    'resnames.append(resn)', space={'resnames': resnames}
)
protein_charge = sum(RESIDUE_CHARGES.get(r, 0) for r in resnames)
total_charge = protein_charge + LIGAND_CHARGE

print(f"\nCharge breakdown:")
for r in resnames:
    q = RESIDUE_CHARGES.get(r, 0)
    if q != 0:
        print(f"  {r}: {q:+d}")
print(f"  Ligand ({LIGAND_RESN}): {LIGAND_CHARGE:+d}")
print(f"  Total: {total_charge:+d}")
```

### Block: Teardown

```python
cmd.quit()
```

---

## Notes
- Always call `cmd.quit()` at the end or the process hangs.
- `byres` expands distance selections to complete residues — never use a raw distance selection as a cluster boundary.
- Selection atom names are case-sensitive: `name CA` not `name ca`.
- If the PDB already has hydrogens, skip capping hydrogen placement and just save the cluster directly.
