---
name: slurm
description: Generate SLURM batch scripts for HPC job submission. Ask the user about partition, resources, and software environment one question at a time.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# SLURM Script Setup

**Ask each question one at a time.** Suggest defaults based on the calculation type, but let the user confirm.

## Questions

### 1 — Partition and account

Both are cluster-specific — you cannot guess them.

```
What partition and account should the job use?
(If unsure, run `sinfo -s` for partitions and `sacctmgr show assoc user=$USER` for accounts.)
```

### 2 — CPUs and memory

```
How many CPUs and how much memory?
I'd suggest: CPUs: ___    Memory: ___ GB
```

Defaults: XTB → 4 CPUs, 4 GB. ORCA → 8 CPUs, 16 GB.

### 3 — Wall time

```
How much wall time? (HH:MM:SS)
I'd suggest: ___
```

Defaults: XTB → 02:00:00. ORCA → 24:00:00.

### 4 — Software environment

Do not guess module names or paths.

```
How is [XTB / ORCA] set up on your cluster?

  a) Module — what's the module name?
  b) Conda — what's the env name?
  c) Custom path — what's the install directory?
  d) Already in $PATH
```

For (c), add to the script:
```bash
export PATH="/path/to/bin:$PATH"
# For ORCA also:
export LD_LIBRARY_PATH="/path/to/orca:$LD_LIBRARY_PATH"
```

### 5 — GPU (skip unless relevant)

Only ask if the user mentioned a GPU partition or the calculation benefits from it.

### Then generate

Use the template below. Default to 1 node, 1 task, output `%j.out`.

## Template

```bash
#!/bin/bash
#SBATCH --job-name=FILL
#SBATCH --output=%j.out
#SBATCH --partition=FILL
#SBATCH --account=FILL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=FILL
#SBATCH --mem=FILL
#SBATCH --time=FILL

# --- Environment ---
FILL_MODULE_OR_PATH_SETUP

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- Run ---
cd $SLURM_SUBMIT_DIR
FILL_RUN_COMMAND
```

## Notes
- Use `$SLURM_CPUS_PER_TASK` for thread counts, never hardcode.
- ORCA: keep `--ntasks-per-node` in sync with `%pal nprocs` in the input.
- XTB: reads `OMP_NUM_THREADS`, no MPI needed.
