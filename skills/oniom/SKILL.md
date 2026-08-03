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
- **Commit**: `d17b2ea` — "Apply constraints and fixed atoms at the ONIOM wrapper level"
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

---

## Metadynamics

Metadynamics is **iterative** — each geometry step evaluates ONIOM energy + gradient once. The GFN-FF outer region runs on the **entire system** every step regardless of inner size. System truncation is therefore mandatory.

**Upstream xtb main does not apply metadynamics at the ONIOM wrapper level.** Use the special fork (see Technical Reference below).

For plain energy calculations (single-point or one-shot geometry optimization without bias), use the sections above.

---

### Physics of `static=true` metadynamics

Reading `src/metadynamic.f90` and `src/dynamic.f90` reveals how the bias actually works:

#### `static=true` with `--opt` (repulsive wall)

All `save=N` "structures" are copies of the same reference geometry. The bias is a **single Gaussian** with amplitude `N × factor`:

```
bias(rmsd) = (N × factor) × exp(-width × rmsd²)
```

This is a **massive repulsive wall** at the reference geometry. The optimizer pushes the system away until physical forces balance the bias. `save` and `factor` are **not independent controls** — they just multiply.

| Parameter | Role | To push farther |
|-----------|------|-----------------|
| `save × factor` | **Amplitude** at rmsd = 0 (Eh) | Increase |
| `width` | **Range** — inverse Gaussian variance (Å⁻²) | Decrease |

| Variant | `save` | `factor` | `width` | Amplitude (Eh) | Push to ~1 Eh |
|---------|--------|----------|---------|---------------|---------------|
| `weak` | 50 | 0.5 | 1.0 | 25 | ~1.7 Å |
| `med` | 200 | 0.5 | 0.5 | 100 | ~2.2 Å |
| `strong` | 500 | 1.0 | 0.5 | 500 | ~2.8 Å |
| `long` | 1000 | 1.0 | 0.2 | 1000 | ~4.3 Å |

*Push distance ≈ RMSD where bias ≈ 1×10⁻²⁰ Eh. Formula: `rmsd = sqrt(ln(amplitude) / width)`.*

#### `static=false` with `--md` (dynamic metadynamics)

Hills are deposited during MD in a FIFO buffer. Each new hill is ramped gradually. This is true history-dependent metadynamics that can fill free-energy basins.

**ONIOM + dynamic MD metadynamics is compatible** — confirmed on trypsin-BEN (472 atoms, GFN2/GFN-FF) with commit `d17b2ea`. The ONIOM wrapper calls `metadynamic()` on the full-system coordinates.

**Parameter mapping:** In dynamic mode, `factor` maps to `kpush` and `width` maps to `alpha`. The log printout shows:
```
--- metadynamics parameter ---
 kpush  :    0.500   (from factor)
 alpha  :    1.000   (from width)
 update :     50     (from save)
```

**No `$fix` during MD:** xTB explicitly disables `$fix` for MD runs (`src/prog/main.F90:1142` resets `fixset%n = 0` before calling the MD driver). Do not use `$fix` with `--md` — the shell atoms will move and the trajectory is invalid.

**Use `$wall` confinement instead:** To keep the truncated system from drifting, add a logfermi sphere centered on the ligand (or binding site):
```
$wall
potential=logfermi
sphere: 15.0, all   # adjust radius to cover inner region + buffer
temp=1000
$end
```
The `$wall` interacts with ALL atoms and applies a counter-force at the boundary, which is stable during dynamics. If shell atoms still drift too much, increase the sphere radius or use the full system without truncation.

| Strategy | Mode | Controls | Use case | Trajectory length |
|----------|------|----------|----------|-------------------|
| **Repulsive wall** | `--opt` + `static=true` | `save×factor` = amplitude, `width` = range | Rapidly find new minima | ~50–500 frames |
| **Well-tempered MTD** | `--md` (omit `static`) | `factor`= hill height, `save` = deposit interval, `width` = sigma | True free-energy mapping | 1,000–10,000+ frames |

---

### Phase 0 — Analyze the site before hypothesizing

Do not guess which residues matter. Look at the structure first.

#### Step 0.1: Identify the ligand

| Property | How to determine |
|----------|-----------------|
| Residue name | First non-standard resn in the PDB (e.g. BEN, NAG, ATP) |
| Formal charge | Sum protonation states of ionizable groups |
| Key functional groups | `cmd.iterate` atom names; look for carboxylate, ammonium, guanidinium, phosphate, amide |

#### Step 0.2: Map the pocket with heavy-atom distances

For every protein residue, compute the **minimum distance from any ligand heavy atom to any residue heavy atom** (exclude hydrogen).

> **Why heavy atoms only?** `cmd.h_add()` places hydrogens algorithmically; they can spuriously pull residues closer. On trypsin-BEN, His57 is 3.98 Å with H but 5.53 Å with heavy atoms only. Asp102 is 7.39 Å with H, 7.96 Å with heavy atoms.

```python
from pymol import cmd
import math

# ... after loading, removing solvent, and h_add ...
lig_model = cmd.get_model('system and resn BEN')
lig_heavy = [(a.coord[0], a.coord[1], a.coord[2]) for a in lig_model.atom
             if a.name[0] not in ('H', 'D')]

atom_elems = {}
cmd.iterate('system', 'atom_elems[index] = elem', space={'atom_elems': atom_elems})

def residue_min_heavy_dist(resi):
    model = cmd.get_model(f'system and resi {resi}')
    heavy = [a for a in model.atom if atom_elems.get(a.index, 'H') != 'H']
    if not heavy:
        return float('inf')
    min_d = float('inf')
    for a in heavy:
        for lx, ly, lz in lig_heavy:
            d = math.sqrt((a.coord[0]-lx)**2 + (a.coord[1]-ly)**2 + (a.coord[2]-lz)**2)
            if d < min_d:
                min_d = d
    return min_d
```

Classify residues into zones:

| Zone | Distance criterion | Purpose in system |
|------|--------------------|-------------------|
| **Inner (QM)** | ≤ 5 Å from ligand OR selected for a hypothesis | High-level GFN2 |
| **Active (MM, moving)** | ≤ 6 Å from ligand, not inner | GFN-FF, unrestrained |
| **Frozen (MM, fixed)** | > 6 Å and ≤ truncation cutoff | GFN-FF, `$fix` to prevent collapse |
| **Excluded** | > truncation cutoff | Removed from system entirely |

**Use atom-level distances, NOT `byres` distances.** `cmd.select('byres (within 12 of ...)')` keeps whole residues when a single sidechain atom is within 12 Å. This overcounts by ~2.5×. On trypsin-BEN, `byres` 12 Å keeps 1311 atoms; atom-distance-based truncation at 8 Å keeps **472 atoms**.

> **Example — trypsin with benzamidine:**
> - ≤5 Å (inner candidates): 17 residues, 215 atoms
> - 5–6 Å (active): 5 residues, 79 atoms (e.g. Lys224 5.14 Å, His57 5.53 Å, Tyr172 5.69 Å)
> - 6–8 Å (frozen): 10 residues, 160 atoms (e.g. Asp102 7.96 Å, Leu158 7.98 Å)
> - >8 Å (exclude): everything else = 3106 atoms
> - Truncated total at 8 Å: **472 atoms** (vs. 3238 full, vs. 1311 `byres` 12 Å)

#### Step 0.3: Build the pocket map

```python
# After loading with h_add, compute for every residue:
pocket = []  # list of (resi, resn, min_dist, atom_count)
# ... fill from atom-level heavy-atom distances ...
pocket.sort(key=lambda x: x[2])

for resi, resn, dist, n_atoms in pocket:
    if dist <= 10.0:
        zone = "inner" if dist <= 5.0 else ("active" if dist <= 6.0 else "buffer")
        print(f"  {resn:>5}{resi:>5}  {dist:>6.2f} Å  ({n_atoms} atoms)  [{zone}]")
```

---

### Phase 1 — Generate hypotheses from the pocket map

Hypotheses must be grounded in **observed proximity** and **known chemistry**. Do NOT invoke mechanistic labels unless the involved residues are actually within interaction distance.

#### Required hypotheses

| Hypothesis | What to look for in pocket map | Inner residues (example: trypsin-BEN) |
|-----------|--------------------------------------|---------------------------------------|
| **H1 Salt bridge** | Charged ligand group + Asp/Glu within 5–8 Å | BEN + Asp189 (2.87 Å) |
| **H2 Hydrophobic pocket** | Aromatic / aliphatic wall around ligand | BEN + Asp189 + Trp215 (3.84 Å) |
| **H3 Nearby catalytic** | Ser/Thr/His/backbone within 5–8 Å that could H-bond | BEN + Asp189 + Ser195 + Gly193 (≤5.5 Å) |
| **H4 Sidechain rotation** | Sidechain reaching toward ligand, might rotate | Optional — pick from pocket map |

#### Validation rule

After generating hypotheses, **verify**: every residue in the inner list must be ≤ 10 Å from the ligand. If a residue is farther:
- Either it cannot physically interact — **drop it**
- Or expand truncation radius and accept the cost

**Never silently include a residue just because it's in a textbook mechanism.**

> **Bad example:** "BEN + His57 + Asp102 + Asp189 + Ser195" as "catalytic triad." Asp102 is 7.96 Å from BEN, His57 is 5.53 Å — both are outside a compact inner-QM region.
>
> **Good example:** "BEN + Asp189 + Ser195 + Gly193" — all within 2.9–5.5 Å.

---

### Phase 2 — Build the truncated system per hypothesis

#### Truncation radius

Set by the **farthest hypothesis residue** plus a **2 Å buffer** (minimum 8 Å):

```python
max_hypothesis_dist = max(dist[resi] for resi in hypothesis_residues)
TRUNCATE_CUTOFF = max(max_hypothesis_dist + 2.0, 8.0)
```

In practice this is almost always **8–10 Å**.

#### Atom-level selection (not byres)

Select atoms individually by distance — residue boundaries do not matter for GFN-FF.

```python
keep_resis = {resi for resi, resn in all_residues
              if residue_distances[resi] <= TRUNCATE_CUTOFF}
```

#### Freeze shell

Freeze atoms in the buffer zone to prevent pocket collapse during metadynamics. Residues whose closest heavy atom is **> 6.0 Å** from the ligand are frozen.

```python
freeze_resis = {r for r in keep_resis
                if residue_distances[r] > 6.0 and r not in hypothesis_residues}
```

---

### Phase 3 — Write per-hypothesis inputs

#### inner.txt

Same format as the ONIOM energy blocks above: `fmt_xtb(sorted(inner_atoms))`

#### xcontrol

**For `--opt` (biased geometry optimization):**

```
$fix
atoms: {freeze_atoms}
$end

$metadyn
save={save}
factor={factor}
width={width}
coord=reference.xyz
atoms: {bias_atoms}
static=true
$end
```

**For `--md` (dynamic metadynamics — preferred for trajectory length):**

```
$wall
potential=logfermi
sphere: 15.0, all   # adjust radius to cover inner region + buffer
temp=1000
$end

$metadyn
save={save}
factor={factor}
width={width}
$end

$md
temp=300
time=100.0
step=1.0
dump=100.0
$end
```

> **No `$fix` for `--md`:** xTB disables `$fix` during MD (`main.F90:1142` resets `fixset%n = 0`). Use `$wall` instead to prevent the truncated shell from drifting.
>
> **Time and dump:** `time=100.0` = 100 ps, `dump=100.0` = 100 fs snapshot interval → ~1,000 frames. Scale `time` up to 1,000 ps (1 ns) for thorough exploration. `dump` should be 50–100 fs to keep file sizes reasonable.
>
> **No `coord` or `atoms`:** Dynamic metadynamics uses the full system geometry (or a default RMSD collective variable). The `coord=` and `atoms:` keywords are only meaningful with `static=true`.

#### Bias-atom selection

The `$metadyn atoms:` defines the RMSD Gaussian reference. **Always include the full ligand** unless studying sidechain-only motion.

| Hypothesis type | Extra bias atoms | Rationale |
|-----------------|-----------------|-----------|
| Salt bridge | Asp/Glu carboxylate O atoms | Track Asp rotation relative to ligand |
| Hydrophobic pocket | Aromatic ring heavy atoms (~6-membered) | Track π-stacking displacement |
| Nearby H-bond | H-bond donor/acceptor sidechain atoms | Track H-bond distance/orientation |
| Sidechain rotation | Full sidechain of the rotating residue | Track torsional reorganization |

**Bias fewer than ~25 atoms total.** More atoms = broader Gaussian = weaker push.

---

### Charges

- Pass `--chrg N` on the CLI (e.g. `--chrg 0` for the truncated trypsin-BEN system).
- **Recompute total charge after truncation.** Removal of charged surface residues changes the total. On full trypsin the charge is +7; after 8 Å truncation it drops to **0**.
- Inner charge: auto-computed by xtb ONIOM via GFN-FF EEQ partial charges.
- **Never use `--chrg inner:total`** — colon syntax triggers argument parser error.

---

### Threading rules

xtb links **OpenBLAS** (not MKL). OpenBLAS spawns its own pthreads independently of OpenMP, causing oversubscription when OpenMP threads > 1.

```bash
# For geometry optimization / metadynamics (topology cached after step 1)
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
# In Slurm: -c 4 --mem=16G

# For single-point energy (no geopt)
export OMP_NUM_THREADS=1
# Leave OPENBLAS_NUM_THREADS unset (default)
# In Slurm: -c 1 --mem=8G
```

#### Why 4 cores?

After topology caching, each step evaluates GFN-FF (threaded gradients) + GFN2 SCC (threaded integral/H1 build). On systems < 1500 atoms, >4 OpenMP threads causes oversubscription with OpenBLAS and slows the run. Measured on trypsin (1311 atoms, cached topology): n=4 = 13.4 s, n=8 = 24.0 s.

---

### Core allocation (measured, HBW2 mpc nodes)

#### Per-evaluation cost (cached topology)

One ONIOM energy + gradient evaluation after GFN-FF topology is cached on disk:

| Total atoms | Inner atoms | Truncation | GFN-FF outer | GFN-FF inner | GFN2 SCC | Total / eval | Threads |
|-------------|-------------|------------|--------------|--------------|----------|--------------|---------|
| 3238 (full) | 30 | none | ~0.5 s* | ~0.0 s* | ~0.1 s* | ~0.6 s | 1 |
| 1311 (`byres` 12Å) | 30 | `byres` 12Å | **0.312 s** | **0.000 s** | **0.132 s** | **0.444 s** | 8 (default OpenBLAS) |
| 472 (atom 8Å) | 54 | atom 8Å | **0.035 s** | **0.001 s** | **0.06–0.13 s** | **0.10–0.17 s** | 4 |
| 472 (atom 8Å) | 30 | atom 8Å | **0.032–0.070 s** | **0.000 s** | **0.015–2.27 s** | **0.05–2.34 s** | 4 |

\* Estimated from scaling; not directly measured.

> **Note on h1_saltbridge (30 inner) anomaly:** The first GFN2 SCC for a very small inner region can take **2.27 s** if the initial density guess is poor. Subsequent cycles reuse cached charges and drop to **~0.015 s**. In production metadynamics this only affects step 1; the remaining steps are fast.

#### Single-point benchmark (includes topology build)

| Total atoms | Inner atoms | Truncation | Total wall | SCF wall | Threads |
|-------------|-------------|------------|------------|----------|---------|
| 3238 (full) | 30 | none | 4 m 41 s | 4.9 s | 1 |
| 1311 (`byres` 12Å) | 30 | `byres` 12Å | **20.3 s** | 0.7 s | 8 |
| 472 (atom 8Å) | 54 | atom 8Å | **~3.6 s** | 0.19 s | 4 |

#### Metadynamics total-time estimate

A 50-step static metadynamics evaluates energy+gradient once per step. First step builds and caches GFN-FF topology. ANC optimizer overhead is included.

| Truncation | Total atoms | Est. total time |
|------------|-------------|-----------------|
| None (full) | ~3238 | ~5–7 min |
| `byres` 12Å | ~1311 | ~1.5 min |
| Atom 8Å | ~472 | **~1 min** |

**Key insight:** cost scales with **total atoms**, not inner atoms. The atom-level truncation is **~2.8× smaller** than `byres` 12Å and runs **~4–8× faster per step**.

---

### Phase 4 — Production run: one directory per hypothesis

#### Directory convention

```
my_project/
├── h1_saltbridge/
│   ├── system.xyz
│   ├── inner.txt
│   ├── xcontrol
│   ├── reference.xyz
│   ├── charge.txt
│   └── info.txt
├── h2_pocket/
│   └── ...
└── submit_all.sh
```

#### Submit script

**Recommended: `--md` for trajectory-based exploration (longer wall time)**

```bash
#!/bin/bash -l
#SBATCH -J my_project-mtd
#SBATCH -p mpc
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH --account FILL

module load gcc
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
XTB=/path/to/xtb
BASE=/path/to/my_project

cd "$BASE"

for cfg in h1_saltbridge h2_pocket h3_local_func; do
    CHRG=$(cat ${cfg}/charge.txt)
    INNER=$(cat ${cfg}/inner.txt)
    mkdir -p "run_${cfg}"
    cp ${cfg}/{system.xyz,inner.txt,xcontrol,reference.xyz} "run_${cfg}/"
    cd "run_${cfg}"
    $XTB system.xyz --oniom gfn2:gfnff "$INNER" --input xcontrol --chrg "$CHRG" --md --gfnff 2>&1 | tee "${cfg}.log"
    cd ..
done
```

**Quick test: `--opt` for repulsive-wall exploration (short wall time)**

```bash
#SBATCH -t 00:20:00
# ... same setup ...
$XTB system.xyz --oniom gfn2:gfnff "$INNER" --input xcontrol --chrg "$CHRG" --opt 2>&1 | tee "${cfg}.log"
```

> **CRITICAL:** Each hypothesis gets its own `run_*` directory. This avoids the GFN-FF topology cache bug that hangs sequential runs in the same directory.
>
> **Wall-time scaling ( `--md` ):**
> | `time` in xcontrol | Steps | Est. wall time (472 atoms, 4 cores) |
> |--------------------|-------|---------------------------------------|
> | 10 ps | 10,000 | ~10 min |
> | 100 ps | 100,000 | ~1.5 h |
> | 1,000 ps (1 ns) | 1,000,000 | ~15 h |
>
> Use `mpc` for quick tests (1 h max), `mpc_l` for production (up to 72 h).

---

### Phase 5 — Fetch results

#### `--opt` outputs

| File | Size | Purpose |
|------|------|---------|
| `run_*/xtbopt.log` | ~100 KB | **Trajectory** (multi-frame XYZ) — load in PyMOL |
| `run_*/xtbopt.xyz` | ~35 KB | Final optimized geometry |
| `run_*/*.log` | ~180 KB | Full xtb stdout (energies, timing) |
| `run_*/.xtboptok` | 0 B | Completion marker |

#### `--md` outputs

| File | Size | Purpose |
|------|------|---------|
| `run_*/xtb.trj` | ~1–5 MB | **MD trajectory** (multi-frame XYZ, 1,000–10,000 frames) |
| `run_*/mdt.log` | ~200 KB–2 MB | Full stdout (energies, metadynamics parameters, timing) |
| `run_*/gfnff_charges` | ~1 KB | Cached EEQ charges |
| `run_*/gfnff_topo` | ~1–2 MB | Cached topology (optional) |

> **No `.xtboptok` for `--md`**: Since MD runs for a fixed simulation time, completion is inferred from the job state, not a disk marker.

```bash
mkdir -p outputs/my_project/{h1_saltbridge,h2_pocket,h3_local_func}
for h in h1_saltbridge h2_pocket h3_local_func; do
    # For --opt:
    # rsync -avz user@cluster:/path/to/my_project/run_${h}/{xtbopt.log,xtbopt.xyz,*.log,.xtboptok} \
    #     outputs/my_project/${h}/
    # For --md:
    rsync -avz user@cluster:/path/to/my_project/run_${h}/{xtb.trj,mdt.log,*.log,gfnff_charges,gfnff_topo} \
        outputs/my_project/${h}/
done
```

---

### Phase 6 — Visualize

#### Load trajectory and compare to reference

**For `--opt` (optimization trajectory):**

```python
from pymol import cmd

PROJECT = '/path/to/outputs/my_project'
H = 'h1_saltbridge'

cmd.load(f'{PROJECT}/{H}/reference.xyz', 'ref')
cmd.load(f'{PROJECT}/{H}/xtbopt.log', 'traj')   # multi-frame XYZ
cmd.load(f'{PROJECT}/{H}/xtbopt.xyz', 'final')
cmd.align('traj', 'ref')
cmd.align('final', 'ref')

cmd.show('sticks', 'ref and not elem H')
cmd.color('gray', 'ref')
cmd.show('sticks', 'traj and not elem H')
cmd.color('cyan', 'traj')
cmd.show('spheres', 'traj and resn BEN')

# Check if ligand moved
print("=== Did shit happen? ===")
cmd.rms_cur('traj', 'ref', mobile=f'traj and resn BEN', target=f'ref and resn BEN')
```

**For `--md` (MD trajectory, much longer):**

```python
from pymol import cmd

PROJECT = '/path/to/outputs/my_project'
H = 'h1_saltbridge'

cmd.load(f'{PROJECT}/{H}/reference.xyz', 'ref')
cmd.load(f'{PROJECT}/{H}/xtb.trj', 'md_traj')   # 1,000–10,000 frames
cmd.align('md_traj', 'ref')

cmd.show('sticks', 'ref and not elem H')
cmd.color('gray', 'ref')
cmd.show('sticks', 'md_traj and not elem H')
cmd.color('cyan', 'md_traj')
cmd.show('spheres', 'md_traj and resn BEN')

# Frame 1 = start of MD; last frame = end of simulation
print(f"Total frames: {cmd.count_states('md_traj')}")
cmd.rms_cur('md_traj', 'ref', mobile=f'md_traj and resn BEN', target=f'ref and resn BEN')
```

#### What to look for

| Question | How to check |
|----------|--------------|
| Did the ligand move? | `rms_cur` between `traj` and `ref` on ligand atoms |
| Did the salt bridge break? | `distance` between amidinium N and Asp189 OD across frames |
| Did the pocket wall shift? | `rms_cur` on Trp215 (or other aromatic) heavy atoms |
| Did the bias push too hard? | Energy trend in `*.log` — should oscillate, not diverge |
| Is the pocket collapsing? | `rms_cur` on frozen buffer atoms should be ~0; if >0.3 Å, truncation too aggressive |

---

### Phase 7 — Parameter sweep (aggressive exploration)

To discover how far a given hypothesis can push the system, run a **parameter sweep** across amplitudes.

#### Static bias (`--opt`)

Use `coord=` + `atoms:` to bias toward a specific geometry:

```python
variants = [
    ('ext_weak',  {'save': 50,  'factor': 0.5, 'width': 1.0}),
    ('ext_med',   {'save': 200, 'factor': 0.5, 'width': 0.5}),
    ('ext_strong',{'save': 500, 'factor': 1.0, 'width': 0.5}),
    ('ext_long',  {'save': 1000,'factor': 1.0, 'width': 0.2}),
]
```

Each variant generates an `xcontrol` with the same `atoms:` bias but different `save`/`factor`/`width`.

#### Dynamic metadynamics (`--md`)

Omit `coord` and `atoms:` so hills deposit on the fly during MD:

```python
variants = [
    # Gentle exploration (mcp-reactor default)
    ('md_gentle',   {'save': 50,  'factor': 0.01, 'width': 1.0}),
    # Moderate pushing
    ('md_mod',      {'save': 50,  'factor': 0.1,  'width': 1.0}),
    # Strong exploration
    ('md_strong',   {'save': 50,  'factor': 0.5,  'width': 1.0}),
    # Very aggressive
    ('md_aggro',    {'save': 50,  'factor': 1.0,  'width': 0.5}),
]
```

**Mapping**: `factor` → `kpush`, `width` → `alpha`, `save` → `update` (hill deposit interval in steps).

Group into batches based on expected wall time (see Phase 4 wall-time table).

---

### Phase 8 — Path chaining and substrate swap

#### Path chaining

After Phase 1, use endpoint geometries as **new references** and bias toward the NEXT expected intermediate:

1. Save `xtbopt.xyz` as `reference_step2.xyz`
2. Design new bias atoms that push toward the next expected state
3. Run metadynamics with the new reference

This creates a **chain of metadynamics runs**: each link pushes from one intermediate to the next.

#### Substrate swap (for full catalytic cycles)

If the ligand is an inhibitor (e.g. BEN in trypsin), the full catalytic cycle requires a substrate analog. Options:
- Modify the inhibitor in PyMOL to add a scissile peptide bond
- Search for PDBs with tetrahedral intermediate analogs (e.g. trypsin-DIP, trypsin-TAME)
- Use existing covalent inhibitor structures as TS-like references

After substrate preparation, regenerate `system.xyz`, recompute charges, and re-run the full hypothesis pipeline.

---

## Technical Reference

### Special xtb Fork

Upstream xtb main does **not** correctly apply `$fix`/`$constrain` and metadynamics at the ONIOM wrapper level with correct atom indices. This causes silent failures or index-out-of-bounds errors.

Use this fork:

- **Fork**: `https://github.com/william-dawson/xtb` (branch `main`)
- **Commit**: `d17b2ea` — "Apply constraints and fixed atoms at the ONIOM wrapper level"
- **Build** (CMake):
  ```bash
  git clone https://github.com/william-dawson/xtb.git
  cd xtb
  git checkout d17b2ea
  mkdir build && cd build
  cmake .. -DWITH_TBLITE=true
  cmake --build . -j$(nproc)
  ```

This fix moves `constrain_pot`, `constrpot`, `cavity_egrad`, and `metadynamic` calls to the ONIOM wrapper (`src/oniom.f90:~497`) where they operate on the full-system coordinate frame with correct atom indices. Fixed-atom zeroing also happens once at the ONIOM level after all sub-calculations.

---

### Fixed Atoms & Constraints

#### 1. Write xcontrol (or xtb.inp)

```
$fix
atoms: 30-289
$end
```

**CRITICAL**: **NO leading whitespace** before keywords. The parser does `trim(line(:ie-1))` which only strips trailing spaces. Lines like `  atoms: 30-289` silently fail.

#### 2. Pass --input explicitly

```bash
xtb --input xcontrol --oniom gfn2:gfnff "1-8" --grad cluster.pdb
```

xTB does **not** autoload `xcontrol` from the working directory; `--input` is mandatory.

#### 3. Verify in the gradient file

After the run, fixed atoms must show exactly `0.00000000` for all x, y, z components:

```bash
awk 'NR>=293+29 && NR<=293+288{print}' gradient | head
```

---

### Internal metadynamics mechanism

`metadynamic(metaset, nat, at, xyz, ebias, gradient)` (`src/metadynamic.f90:21-89`) computes:

```
ebias += Σ_i factor(i) × exp(-width(i) × rmsd_i²)
```

The gradient contribution is:

```
g += -2 × width × factor × exp(-width × rmsd²) × rmsd × (∂rmsd/∂xyz)
```

### Validation test pattern

The unit test `test_oniom_metadynamics` in `test/unit/test_oniom.f90` validates the ONIOM wrapper path: with reference = current geometry and `factor=0.5`, the energy increment must equal `0.5` Eh within `1e-6`.

### Programmatic cleanup

After using metadynamics programmatically, clean up with:

```fortran
call clear_metadyn
metaset%nstruc = 0
metaset%maxsave = 0
```

`clear_metadyn` deallocates arrays. Zeroing `nstruc` and `maxsave` prevents residual state from being picked up by subsequent calculations, because `metadynamic()` returns early when `nstruc < 1`.

---

### Environment for Large Systems

If you encounter stack-overflow crashes (common with 5000+ atoms and OpenMP), set:

```bash
export OMP_STACKSIZE=4G     # or 2G, 8G as needed
export OMP_NUM_THREADS=2
```

---

### Common Gotchas

- **Index conventions**: XTB uses **1-based** atom indices in CLI and `inner.txt`. ORCA `QMATOMS` uses **0-based** indices.
- **Link atoms**: XTB places link atoms at cut bonds automatically. **Only cut single bonds.** Use `xtb ... --cut` to verify before running.
- **xcontrol autoload**: xTB never autoloads `xcontrol`. `--input` is mandatory.
- **Whitespace sensitivity**: Keywords in `$fix`, `$constrain`, and `$metadyn` must start at column 1. Leading spaces silently fail.
- **`$fix` is ignored during `--md`:** xTB resets `fixset%n = 0` before calling the MD driver (`src/prog/main.F90:1142`). `$fix` atoms will move freely. Use `$wall` confinement instead.
- **Colon syntax**: Never use `--chrg inner:total` — it may trigger the argument parser's help screen.


