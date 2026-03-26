---
name: slurm
description: Generate SLURM batch scripts for submitting calculations to HPC clusters. Guides the user through partition, resource, and module selection. Use when the user wants to run a calculation on a cluster via SLURM.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# SLURM Batch Script Setup

Your job is to write a SLURM submission script tailored to the user's cluster and calculation. Present the setup form below, fill in what you can, and ask the user for the rest.

## Setup Form

```
═══════════════════════════════════════════════════════════
  SLURM SUBMISSION — SETUP
═══════════════════════════════════════════════════════════

─── Job identity ──────────────────────────────────────────

  Job name:       _______________
  Output file:    [ ] %j.out (default)    [ ] _______________
  Error file:     [ ] same as output      [ ] _______________

─── Resources ─────────────────────────────────────────────

  Partition:      _______________  (ask user — cluster-specific)
  Nodes:          [ ] 1           [ ] ___
  Tasks per node: [ ] 1           [ ] ___
  CPUs per task:  [ ] 1           [ ] ___
  Memory:         [ ] default     [ ] ___ GB
  GPU:            [ ] none        [ ] ___ (e.g. gpu:1, gpu:a100:2)
  Wall time:      ___:___:___     (HH:MM:SS)

─── Environment ───────────────────────────────────────────

  Modules:        _______________  (e.g. xtb/6.6, orca/5.0)
  Conda env:      [ ] none        [ ] _______________
  Extra env vars: _______________

─── Execution ─────────────────────────────────────────────

  Working dir:    [ ] submission directory  [ ] _______________
  Run command:    _______________

═══════════════════════════════════════════════════════════
```

### How to fill it in

**Partition** — This is always cluster-specific. Ask the user. Common names: `batch`, `normal`, `compute`, `gpu`, `short`, `long`.

**Resources** — Sensible defaults depend on the calculation:
- XTB ONIOM: 1 node, 1 task, 4–8 CPUs, ~4 GB, no GPU, 1–4 hours
- ORCA ONIOM: 1 node, 1 task, 8–16 CPUs, ~16 GB, no GPU, 12–48 hours
- Pure XTB single-point: 1 node, 1 task, 4 CPUs, ~2 GB, no GPU, < 1 hour

Suggest defaults based on the calculation but always let the user override.

**Modules** — Ask the user what module system their cluster uses and what the module names are. Do not guess module names — they are cluster-specific.

**Conda env** — If the user installed XTB or other tools via conda, ask for the environment name.

**Run command** — This comes from the calculation setup (e.g. `./run.sh` from the oniom skill, or a direct `xtb` / `orca` command).

## Script Template

```bash
#!/bin/bash
#SBATCH --job-name=FILL_JOB_NAME
#SBATCH --output=FILL_OUTPUT
#SBATCH --error=FILL_ERROR
#SBATCH --partition=FILL_PARTITION
#SBATCH --nodes=FILL_NODES
#SBATCH --ntasks-per-node=FILL_NTASKS
#SBATCH --cpus-per-task=FILL_CPUS
#SBATCH --mem=FILL_MEM
#SBATCH --time=FILL_WALLTIME
# SBATCH --gres=FILL_GPU            # uncomment if GPU needed

# --- Environment ---
module purge
module load FILL_MODULES

# If using conda:
# source activate FILL_CONDA_ENV

# Set parallelism
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- Run ---
cd $SLURM_SUBMIT_DIR
FILL_RUN_COMMAND
```

## Notes

- Always use `$SLURM_CPUS_PER_TASK` for `OMP_NUM_THREADS` rather than hardcoding — this ensures the thread count matches what SLURM allocated.
- ORCA requires `--ntasks-per-node` to match the `%pal nprocs` setting in the ORCA input. If generating both, keep them in sync.
- XTB reads `OMP_NUM_THREADS` for parallelism. No MPI needed.
- For ORCA parallel runs, the run command is `orca input.inp` (ORCA handles MPI internally via its own `orca` wrapper).
- If the user doesn't know their partition or modules, suggest they run `sinfo -s` (partitions) and `module avail xtb` or `module avail orca` on the cluster.
