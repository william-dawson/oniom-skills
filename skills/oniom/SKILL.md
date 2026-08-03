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
Which ONIOM method?

  ① XTB built-in ONIOM (fast, single command)
     GFN2-xTB (inner) / GFN-FF (outer)

  ② ORCA native QM/XTB (electrostatic embedding)
     You specify the QM method (e.g. r2SCAN-3c, wB97X-D3 def2-TZVP)

I'd recommend ① for a first look. Which do you prefer?
```

If ②, follow up asking which DFT method.

### 4 — Inner region cutoff

```
How large should the inner (high-level) region be?

  3.5 Å — ligand + immediate contacts
  4.5 Å — compact active site (recommended)
  5.0 Å — standard active site

What cutoff would you like?
```

### 5 — Binding energy (3-point calculation)

```
Do you want a binding energy calculation?

  a) No — just a single-point energy of the complex
  b) Yes — compute ΔE_bind = E(complex) − E(protein) − E(ligand)
     (generates three calculations using the same ONIOM setup)
```

If (b), three inputs are generated using the same method and inner region definition:
- **Complex**: protein + ligand (the full system)
- **Protein**: ligand removed
- **Ligand**: ligand alone (no ONIOM needed — just a single-point)

### 6 — Execution

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

## XTB Binary: Fork & Build

Upstream xTB main does **not** yet apply `$fix`/`$constrain` at the ONIOM wrapper level. Use this fork:

- **Fork**: `https://github.com/william-dawson/xtb`  (branch `main`)
- **Commit**: `a3d1be7` — "Apply constraints and fixed atoms at the ONIOM wrapper level"
- **Build** (CMake):
  ```bash
  mkdir build && cd build
  cmake .. -DWITH_TBLITE=true
  make -j$(nproc)
  ```

This fix zeroes gradients on fixed atoms once at the ONIOM level instead of inside each sub-calculator (where atom indices are wrong).

---

## Fixed Atoms & Constraints with ONIOM

To freeze part of the protein during ONIOM (e.g. keep the outer region rigid while optimizing the active site):

### 1. Write `xcontrol` (or `xtb.inp`)

```
$fix
atoms: 30-289
$end
```

**CRITICAL**: **NO leading whitespace/indentation** before keywords. The parser does `trim(line(:ie-1))` which only strips trailing spaces. Lines like `  atoms: 30-289` silently fail.

### 2. Pass `--input` explicitly

```bash
xtb --input xcontrol --oniom gfn2:gfnff "1-8" --grad cluster.pdb
```

xTB does **not** autoload `xcontrol` from the working directory; `--input` is mandatory.

### 3. Verify in the `gradient` file

After the run, fixed atoms must show exactly `0.00000000` for all x, y, z components:

```bash
awk 'NR>=293+29 && NR<=293+288{print}' gradient | head
```

### 4. Environment for large systems

If you encounter stack-overflow crashes (common with 5000+ atoms and OpenMP), set:

```bash
export OMP_STACKSIZE=4G     # or 2G, 8G as needed
export OMP_NUM_THREADS=2
```

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

### 3-Point Binding Energy Blocks

Use when the user selects binding energy in Question 5. Generates three ONIOM calculations using the **same method and inner region definition**:

- **Complex**: protein + ligand (the full system as-is)
- **Protein**: ligand removed from the system
- **Ligand**: ligand alone (no ONIOM — just a single-point)

`ΔE_bind = E(complex) − E(protein) − E(ligand)`

#### Prepare three geometries

```python
# FILL
OUTDIR      = 'binding_energy'
LIGAND_RESN = 'LIG'
LIGAND_CHARGE = 0  # FILL

os.makedirs(OUTDIR, exist_ok=True)

# 1. Complex (full system)
cmd.save(f'{OUTDIR}/complex.xyz', 'system')

# 2. Protein (ligand removed)
cmd.save(f'{OUTDIR}/protein.xyz', f'system and not resn {LIGAND_RESN}')

# 3. Ligand alone
cmd.save(f'{OUTDIR}/ligand.xyz', f'system and resn {LIGAND_RESN}')
```

#### Recompute inner indices for complex and protein

The complex uses the same inner region as before. For the protein-only system, the inner region indices shift because the ligand atoms are gone — recompute them.

```python
# Complex inner indices (already computed above)
# inner_xtb / inner_orca are still valid for the complex

# Protein inner indices — ligand atoms removed, so indices shift
prot_atoms = []
cmd.iterate(f'system and not resn {LIGAND_RESN}',
            'prot_atoms.append(index)', space={'prot_atoms': prot_atoms})
prot_inner_pymol = inner_pymol - set()  # same inner set
# But remove ligand atoms from inner
lig_indices = set()
cmd.iterate(f'system and resn {LIGAND_RESN}',
            'lig_indices.add(index)', space={'lig_indices': lig_indices})
prot_inner_pymol = inner_pymol - lig_indices

# Map to 1-based XYZ positions in the protein-only file
prot_inner_xtb = [i + 1 for i, idx in enumerate(prot_atoms)
                  if idx in prot_inner_pymol]
prot_inner_range = fmt_xtb(prot_inner_xtb) if prot_inner_xtb else ''

# Protein charge (no ligand)
protein_charge = total_charge - LIGAND_CHARGE
prot_inner_charge = inner_charge - LIGAND_CHARGE
```

#### Write run scripts (XTB example — adapt for ORCA)

```python
# Uses the same HIGH, LOW from the method selection
for label, xyz, chrg, inner_chrg, ir in [
    ('complex', 'complex.xyz', total_charge, inner_charge, inner_range),
    ('protein', 'protein.xyz', protein_charge, prot_inner_charge, prot_inner_range),
]:
    xtb_cmd = (f'xtb {xyz} --oniom {HIGH}:{LOW} {ir} '
               f'--chrg {inner_chrg}:{chrg}')
    with open(f'{OUTDIR}/run_{label}.sh', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write(f'# {label}: {HIGH.upper()} / {LOW.upper()}\n')
        f.write('cd "$(dirname "$0")"\n')
        f.write(f'{xtb_cmd} > {label}.out 2>&1\n')
        f.write(f'grep "TOTAL ENERGY" {label}.out\n')
    os.chmod(f'{OUTDIR}/run_{label}.sh', 0o755)

# Ligand — no ONIOM, just a single-point
with open(f'{OUTDIR}/run_ligand.sh', 'w') as f:
    f.write('#!/bin/sh\n')
    f.write(f'# ligand alone at {HIGH.upper()}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write(f'xtb ligand.xyz --{HIGH} --chrg {LIGAND_CHARGE} > ligand.out 2>&1\n')
    f.write('grep "TOTAL ENERGY" ligand.out\n')
os.chmod(f'{OUTDIR}/run_ligand.sh', 0o755)

print(f"Wrote 3 run scripts in {OUTDIR}/")
print("After running, compute: ΔE_bind = E(complex) − E(protein) − E(ligand)")
```

The LLM should parse the output files after the runs complete and report `ΔE_bind`.

---

### Option ② — ORCA Native QM/XTB Blocks

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

