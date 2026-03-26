---
name: oniom
description: Set up a two-layer ONIOM calculation using XTB or ORCA. Guides method, charge, multiplicity, and inner region decisions. Use after running the pymol skill.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# ONIOM Setup

Write input files tailored to the user's system from the building blocks below. **Ask each question one at a time**, waiting for the answer before moving on.

## Prerequisites

If the pymol skill (`/pymol:pymol`) has NOT been run in this conversation, stop and tell the user:

> Run `/pymol:pymol` first to normalize the PDB, identify the ligand, and determine charges. Then come back here.

Do not attempt PDB cleanup in this skill.

---

## Questions (ask one at a time)

### 1 — Charges

```
From the cluster extraction:
  Inner region charge:  ___
  Total system charge:  ___

Does this look right?
```

### 2 — Multiplicity

State your assessment and confirm. Only ask if ambiguous (metals, radicals).

```
This system appears closed-shell → multiplicity = 1 (singlet).
Is that correct?
```

### 3 — Method

```
Which approach?

  ① XTB built-in ONIOM (fast, single command)
     GFN2-xTB (inner) / GFN-FF (outer)

  ② 3-point manual (three separate calculations)
     You pick the high and low methods independently.
     Generates three run scripts — can use different programs.

  ③ ORCA native QM/XTB (electrostatic embedding)
     You specify the QM method (e.g. wB97X-D3 def2-TZVP)

I'd recommend ① for a first look. Which do you prefer?
```

If ② or ③, follow up asking which high-level method (GFN2, r2SCAN-3c, wB97X-D3, etc.).

### 4 — Inner region cutoff

```
How large should the inner (high-level) region be?

  3.5 Å — ligand + immediate contacts
  4.5 Å — compact active site (recommended)
  5.0 Å — standard active site

What cutoff would you like?
```

### 5 — Execution

```
How do you want to run this?

  a) Locally (generates run.sh)
  b) On a SLURM cluster (I'll help set up a batch script)
```

If (b), use `/pymol:slurm`.

### Then generate

Only after all questions are answered, assemble the script.

---

## ONIOM Energy

```
E_oniom = E(whole, low) − E(inner, low) + E(inner, high)
```

XTB places link atoms at cut bonds automatically. **Only cut single bonds.** Use `xtb ... --cut` to verify before running.

---

## Code Building Blocks

**NEVER use `python3` or `python`.** Always run via PyMOL:

```bash
pymol -cq -r your_script.py
```

### Block: Boilerplate

```python
from pymol import cmd
import os
```

### Block: Load and select inner region

```python
# FILL
PDB_PATH    = 'protein_ligand.pdb'
LIGAND_RESN = 'LIG'
CUTOFF      = 4.5

cmd.load(PDB_PATH, 'system')
cmd.select('ligand', f'system and resn {LIGAND_RESN}')
cmd.select('inner',  f'byres (system within {CUTOFF} of ligand) or ligand')
```

### Block: Collect atom order and inner indices

```python
all_atoms = []
cmd.iterate('system', 'all_atoms.append(index)', space={'all_atoms': all_atoms})

inner_pymol = set()
cmd.iterate('inner', 'inner_pymol.add(index)', space={'inner_pymol': inner_pymol})

print(f"Total atoms:  {len(all_atoms)}")
print(f"Inner region: {len(inner_pymol)} atoms")
```

### Block: Compute charges

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

inner_resnames, all_resnames = [], []
cmd.iterate(f'inner and not resn {LIGAND_RESN} and name CA',
            'inner_resnames.append(resn)', space={'inner_resnames': inner_resnames})
cmd.iterate(f'system and not resn {LIGAND_RESN} and name CA',
            'all_resnames.append(resn)', space={'all_resnames': all_resnames})

inner_charge = sum(RESIDUE_CHARGES.get(r, 0) for r in inner_resnames) + LIGAND_CHARGE
total_charge  = sum(RESIDUE_CHARGES.get(r, 0) for r in all_resnames)  + LIGAND_CHARGE
print(f"Inner charge: {inner_charge:+d}   Total charge: {total_charge:+d}")
```

---

### Option ① — XTB Built-in ONIOM Blocks

#### Format XTB index ranges (1-based)

```python
def fmt_xtb(indices):
    indices = sorted(indices)
    ranges, start, prev = [], indices[0], indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
        else:
            ranges.append(f'{start}-{prev}' if start != prev else str(start))
            start = prev = i
    ranges.append(f'{start}-{prev}' if start != prev else str(start))
    return ','.join(ranges)

inner_xtb = [i + 1 for i, idx in enumerate(all_atoms) if idx in inner_pymol]
inner_range = fmt_xtb(inner_xtb)
```

#### Save XYZ and write XTB run script

```python
OUTDIR = 'oniom_xtb'  # FILL
HIGH   = 'gfn2'
LOW    = 'gfnff'

os.makedirs(OUTDIR, exist_ok=True)
cmd.save(f'{OUTDIR}/system.xyz', 'system')

with open(f'{OUTDIR}/inner.txt', 'w') as f:
    f.write(inner_range + '\n')

xtb_cmd = (f'xtb system.xyz --oniom {HIGH}:{LOW} {inner_range} '
           f'--chrg {inner_charge}:{total_charge}')

with open(f'{OUTDIR}/run.sh', 'w') as f:
    f.write('#!/bin/sh\n')
    f.write(f'# {HIGH.upper()} (inner) / {LOW.upper()} (whole)\n')
    f.write(f'# Inner: {len(inner_xtb)} atoms, charge {inner_charge:+d}\n')
    f.write(f'# Total: {len(all_atoms)} atoms, charge {total_charge:+d}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write(xtb_cmd + '\n')
os.chmod(f'{OUTDIR}/run.sh', 0o755)
print(f"Run: cd {OUTDIR} && ./run.sh")
```

#### xcontrol for ORCA high level (optional within option ①)

Only if using ORCA as the high level inside XTB's ONIOM. The ORCA input must include `! engrad`. Add `--input xcontrol` to `xtb_cmd`.

```python
ORCA_INP_PATH = 'orca.inp'  # FILL

with open(f'{OUTDIR}/xcontrol', 'w') as f:
    f.write('$external\n')
    f.write(f'   orca input file={os.path.abspath(ORCA_INP_PATH)}\n')
    f.write('$end\n')
```

---

### Option ② — 3-Point Manual ONIOM Blocks

Three separate calculations, combined as: `E_oniom = E(whole,low) − E(inner,low) + E(inner,high)`.

This uses the **capped cluster PDB** from the pymol skill as the inner region. The cluster already has link atoms (capping H) placed at backbone cuts.

#### Save geometries

```python
OUTDIR      = 'oniom_3pt'    # FILL
CLUSTER_PDB = 'cluster.pdb'  # FILL: capped cluster from pymol skill

os.makedirs(OUTDIR, exist_ok=True)

# Full system XYZ
cmd.save(f'{OUTDIR}/whole.xyz', 'system')

# Inner region XYZ (load the capped cluster, save as XYZ)
cmd.load(CLUSTER_PDB, 'cluster')
cmd.save(f'{OUTDIR}/inner.xyz', 'cluster')
cmd.delete('cluster')
```

#### Generate three run scripts

```python
# FILL
HIGH_METHOD = 'gfn2'    # e.g. 'gfn2', or 'orca' for DFT
LOW_METHOD  = 'gfnff'
MULTIPLICITY = 1

# --- Run 1: E(whole, low) ---
with open(f'{OUTDIR}/run_whole_low.sh', 'w') as f:
    f.write('#!/bin/sh\n')
    f.write(f'# Point 1: whole system at {LOW_METHOD.upper()}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write(f'xtb whole.xyz --{LOW_METHOD} --chrg {total_charge} > whole_low.out 2>&1\n')
    f.write('grep "TOTAL ENERGY" whole_low.out\n')
os.chmod(f'{OUTDIR}/run_whole_low.sh', 0o755)

# --- Run 2: E(inner, low) ---
with open(f'{OUTDIR}/run_inner_low.sh', 'w') as f:
    f.write('#!/bin/sh\n')
    f.write(f'# Point 2: inner region at {LOW_METHOD.upper()}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write(f'xtb inner.xyz --{LOW_METHOD} --chrg {inner_charge} > inner_low.out 2>&1\n')
    f.write('grep "TOTAL ENERGY" inner_low.out\n')
os.chmod(f'{OUTDIR}/run_inner_low.sh', 0o755)

# --- Run 3: E(inner, high) ---
if HIGH_METHOD.startswith('gfn') or HIGH_METHOD == 'gfnff':
    # XTB high level
    with open(f'{OUTDIR}/run_inner_high.sh', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write(f'# Point 3: inner region at {HIGH_METHOD.upper()}\n')
        f.write('cd "$(dirname "$0")"\n')
        f.write(f'xtb inner.xyz --{HIGH_METHOD} --chrg {inner_charge} > inner_high.out 2>&1\n')
        f.write('grep "TOTAL ENERGY" inner_high.out\n')
    os.chmod(f'{OUTDIR}/run_inner_high.sh', 0o755)
else:
    # ORCA high level — generate .inp file
    with open(f'{OUTDIR}/inner.xyz') as xf:
        coord_lines = xf.readlines()[2:]
    with open(f'{OUTDIR}/inner_high.inp', 'w') as f:
        f.write(f'! {HIGH_METHOD}\n\n')
        f.write(f'* XYZ {inner_charge} {MULTIPLICITY}\n')
        f.writelines(coord_lines)
        f.write('*\n')
    with open(f'{OUTDIR}/run_inner_high.sh', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write(f'# Point 3: inner region at {HIGH_METHOD} (ORCA)\n')
        f.write('cd "$(dirname "$0")"\n')
        f.write('orca inner_high.inp > inner_high.out 2>&1\n')
    os.chmod(f'{OUTDIR}/run_inner_high.sh', 0o755)

print(f"Wrote 3 run scripts in {OUTDIR}/")
```

#### Run all three and combine

```bash
cd oniom_3pt
./run_whole_low.sh
./run_inner_low.sh
./run_inner_high.sh
```

After all three finish, extract the `TOTAL ENERGY` from each `.out` file:

```
E_oniom = E(whole,low) − E(inner,low) + E(inner,high)
```

The LLM should parse the output files, extract the energies, and compute `E_oniom` for the user.

---

### Option ③ — ORCA Native QM/XTB Blocks

Indices are **0-based**. Coordinates embedded in the input file. Low level is always XTB2.

#### Format ORCA QMATOMS ranges (0-based)

```python
def fmt_orca(indices):
    indices = sorted(indices)
    ranges, start, prev = [], indices[0], indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
        else:
            ranges.append(f'{start}:{prev}' if start != prev else str(start))
            start = prev = i
    ranges.append(f'{start}:{prev}' if start != prev else str(start))
    return '{' + ' '.join(ranges) + '}'

inner_orca = [i for i, idx in enumerate(all_atoms) if idx in inner_pymol]
qmatoms_str = fmt_orca(inner_orca)
```

#### Write ORCA QM/XTB input

```python
OUTDIR       = 'oniom_orca'  # FILL
QM_METHOD    = 'r2SCAN-3c'   # FILL
MULTIPLICITY = 1              # FILL

os.makedirs(OUTDIR, exist_ok=True)

cmd.save(f'{OUTDIR}/_tmp.xyz', 'system')
with open(f'{OUTDIR}/_tmp.xyz') as f:
    coord_lines = f.readlines()[2:]
os.remove(f'{OUTDIR}/_tmp.xyz')

with open(f'{OUTDIR}/orca_oniom.inp', 'w') as f:
    f.write(f'! QM/XTB {QM_METHOD}\n\n')
    f.write(f'%QMMM\n   QMATOMS {qmatoms_str}\nEND\n\n')
    f.write(f'* XYZ {total_charge} {MULTIPLICITY}\n')
    f.writelines(coord_lines)
    f.write('*\n')

with open(f'{OUTDIR}/run.sh', 'w') as f:
    f.write('#!/bin/sh\n')
    f.write(f'# ORCA QM/XTB: {QM_METHOD} (inner) / XTB2 (outer)\n')
    f.write(f'# QM region: {len(inner_orca)} atoms, total: {len(all_atoms)}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write('orca orca_oniom.inp\n')
os.chmod(f'{OUTDIR}/run.sh', 0o755)
print(f"Run: cd {OUTDIR} && ./run.sh")
```

