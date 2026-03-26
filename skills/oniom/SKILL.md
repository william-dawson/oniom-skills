---
name: oniom
description: Set up and run a two-layer ONIOM calculation on a protein-ligand system, using XTB or ORCA as backends. Guides method selection, charge/multiplicity decisions, and inner region definition. Use after the cluster has been prepared (see the pymol skill).
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# ONIOM Setup (XTB or ORCA backends)

Your job is to write input files tailored to the user's system, assembled from the building blocks below. Reason through the chemistry first — do not silently fill in defaults.

## Before Writing Any Input

Present the user with this setup form, pre-filled with what you know from the cluster extraction step. Let them review and edit before generating input files.

```
═══════════════════════════════════════════════════════════
  ONIOM CALCULATION — SETUP
═══════════════════════════════════════════════════════════

─── Charges (from cluster extraction) ─────────────────────

  Inner region charge:  ___
  Total system charge:  ___

─── Multiplicity ──────────────────────────────────────────

  [ ] Singlet (1) — closed-shell, no unpaired electrons
  [ ] Doublet (2)
  [ ] Triplet (3)
  [ ] Other: ___
  (reasoning: _________________________)

─── Method ────────────────────────────────────────────────

  Driver:  ① XTB   ② XTB+ORCA   ③ ORCA native   →  ___

  ┌─────────────────────────────────────────────────────┐
  │  ① XTB (fast semiempirical)                        │
  │     High: GFN2-xTB    Low: GFN-FF                  │
  │                                                     │
  │  ② XTB + ORCA (DFT inner, XTB drives)              │
  │     High: _______________ (e.g. r2SCAN-3c)         │
  │     Low:  GFN-FF                                    │
  │                                                     │
  │  ③ ORCA native QM/XTB (electrostatic embedding)    │
  │     QM:   _______________ (e.g. wB97X-D3 def2-TZVP)│
  │     Low:  XTB2 (automatic)                          │
  └─────────────────────────────────────────────────────┘

─── Inner region cutoff ───────────────────────────────────

  Cutoff: _______ Å

  (guide: 3.5 Å = ligand + immediate contacts
          4.5 Å = compact active site
          5.0 Å = standard active site)

─── Execution ─────────────────────────────────────────────

  [ ] Run locally (generates run.sh)
  [ ] Submit via SLURM (use /pymol:slurm to configure)

═══════════════════════════════════════════════════════════
```

### How to fill it in

**Charges** — Should be known from the cluster extraction step. Confirm with the user.

**Multiplicity** — Default to singlet only if confident the system is closed-shell. If a transition metal is present, reason about its oxidation state and d-electron count. If uncertain, ask.

**Driver** — Recommend based on the situation:
- Fast scan or large system → XTB driver, GFN2:GFN-FF
- High accuracy needed, ORCA available → ORCA native QM/XTB
- DFT inner but want XTB to handle the ONIOM mechanics → XTB + ORCA

**Inner region cutoff** — Can be the same as the cluster cutoff, or tighter. A tighter cutoff (e.g. 3.5 Å) for the high-level region is common when the cluster cutoff was generous (e.g. 8 Å).

---

## Runtime

The setup script below uses PyMOL's bundled Python:

```bash
PYMOL_PYTHON=$(head -1 "$(which pymol)" | sed 's/#!//')
"$PYMOL_PYTHON" setup_oniom.py
```

---

## ONIOM Energy Scheme

```
E_oniom = E(whole, low) − E(inner, low) + E(inner, high)
```

XTB places link atoms at cut bonds automatically. **Only cut single bonds.**

---

## Code Building Blocks

Assemble these into a script. Fill in the marked values before running.

### Block: Boilerplate

```python
import pymol, os
from pymol import cmd
pymol.finish_launching(['pymol', '-c'])
```

### Block: Load and select inner region

If the PDB is CHARMM format, insert the "Normalize CHARMM PDB" block from the pymol skill immediately after `cmd.load` and before any selections.

```python
# FILL
PDB_PATH    = 'protein_ligand.pdb'
LIGAND_RESN = 'LIG'
CUTOFF      = 5.0  # inner region cutoff in Å

cmd.load(PDB_PATH, 'system')
# >>> INSERT CHARMM NORMALIZATION HERE IF NEEDED (see pymol skill) <<<
cmd.select('ligand', f'system and resn {LIGAND_RESN}')
cmd.select('inner',  f'byres (system within {CUTOFF} of ligand) or ligand')
```

### Block: Collect atom order and inner indices

PyMOL iterates atoms in index order, which matches XYZ file line order.

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
# FILL
LIGAND_CHARGE = 0

RESIDUE_CHARGES = {
    'ASP': -1, 'GLU': -1,
    'ARG': +1, 'LYS': +1,
    'HIS':  0, 'HID':  0, 'HIE':  0, 'HIP': +1,  # AMBER naming
    'HSD':  0, 'HSE':  0, 'HSP': +1,              # CHARMM naming
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

### XTB Driver Blocks

Use these when `driver = 'xtb'`.

#### Block: Format XTB index ranges (1-based)

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

#### Block: Save XYZ and write XTB run script

```python
# FILL
OUTDIR   = 'oniom_xtb'
HIGH     = 'gfn2'   # or 'orca' if using ORCA high level
LOW      = 'gfnff'

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
print(f"Wrote {OUTDIR}/run.sh\nRun: cd {OUTDIR} && ./run.sh")
```

#### Block: xcontrol for ORCA as XTB high level (optional)

Only needed when using ORCA inside the XTB driver. The ORCA input must include `! engrad`.

```python
ORCA_INP_PATH = 'orca.inp'  # FILL: path to your ORCA input file

with open(f'{OUTDIR}/xcontrol', 'w') as f:
    f.write('$external\n')
    f.write(f'   orca input file={os.path.abspath(ORCA_INP_PATH)}\n')
    f.write('$end\n')

# Append --input xcontrol to xtb_cmd before writing run.sh
```

---

### ORCA Driver Blocks

Use these when `driver = 'orca'`. Indices are **0-based**. The full system XYZ is embedded in the input file. Low level is always XTB2.

#### Block: Format ORCA QMATOMS ranges (0-based)

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

#### Block: Write ORCA QM/XTB input with embedded coordinates

```python
# FILL
OUTDIR     = 'oniom_orca'
QM_METHOD  = 'r2SCAN-3c'  # e.g. 'wB97X-D3 def2-TZVP'
MULTIPLICITY = 1  # FILL: 1 = singlet, 2 = doublet, 3 = triplet, etc.

os.makedirs(OUTDIR, exist_ok=True)

# Save XYZ temporarily to get coordinate lines
cmd.save(f'{OUTDIR}/_tmp.xyz', 'system')
with open(f'{OUTDIR}/_tmp.xyz') as f:
    coord_lines = f.readlines()[2:]  # skip atom count and comment
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
    f.write(f'# QM region: {len(inner_orca)} atoms\n')
    f.write(f'# Total: {len(all_atoms)} atoms, charge {total_charge:+d}, mult {MULTIPLICITY}\n')
    f.write('cd "$(dirname "$0")"\n')
    f.write('orca orca_oniom.inp\n')
os.chmod(f'{OUTDIR}/run.sh', 0o755)
print(f"Wrote {OUTDIR}/orca_oniom.inp\nRun: cd {OUTDIR} && ./run.sh")
```

### Block: Teardown

```python
cmd.quit()
```

---

## Verify Before Running

Remind the user to:
1. Check the inner/outer boundary does not cut any double or aromatic bond — use `xtb ... --cut` to inspect.
2. Confirm multiplicity if any metal or radical is present.
3. For ORCA driver: confirm `! engrad` is present if using ORCA inside XTB driver mode.

## Running the Calculation

After generating the input files, ask the user how they want to run:

```
─── Execution ─────────────────────────────────────────────

  [ ] Run locally (just use the generated run.sh)
  [ ] Submit to a SLURM cluster
```

If they choose SLURM, use the **slurm** skill (`/pymol:slurm`) to walk them through partition, resources, and modules, then generate a batch script that calls `./run.sh` (or the appropriate run command) from the output directory.
