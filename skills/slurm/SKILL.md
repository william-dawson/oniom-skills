---
name: slurm
description: Generate SLURM batch scripts for submitting calculations to HPC clusters. Guides the user through partition, resource, and module selection one question at a time. Use when the user wants to run a calculation on a cluster via SLURM.
user-invocable: true
allowed-tools: Read, Write, Bash, Glob
---

# SLURM Batch Script Setup

Your job is to write a SLURM submission script tailored to the user's cluster and calculation. **Ask each question one at a time**, waiting for the user's answer before moving to the next. Do not present all options at once.

## Questions

Ask these in order. Suggest a default where you can, but let the user confirm or change it.

### Question 1 — Partition

This is always cluster-specific. You cannot guess it.

```
What partition should the job run on?

(Common names: batch, normal, compute, gpu, short, long.
 If unsure, you can run `sinfo -s` on the cluster to see available partitions.)
```

### Question 2 — CPUs and memory

Suggest a default based on the calculation type, then ask:

```
How many CPUs and how much memory?

For this [XTB ONIOM / ORCA / etc.] calculation I'd suggest:
  CPUs: ___    Memory: ___ GB

Does that work, or would you like different values?
```

Guidelines for your suggestion:
- XTB ONIOM: 4–8 CPUs, ~4 GB
- ORCA ONIOM: 8–16 CPUs, ~16 GB
- Pure XTB single-point: 4 CPUs, ~2 GB

### Question 3 — Wall time

```
How much wall time should I request? (HH:MM:SS)

For this calculation I'd suggest: ___
```

Guidelines:
- XTB ONIOM: 1–4 hours
- ORCA ONIOM: 12–48 hours
- Pure XTB: < 1 hour

### Question 4 — Software environment

Do not guess module names or install paths — they are cluster-specific. XTB and ORCA are often installed manually rather than as modules.

```
How is [XTB / ORCA / etc.] set up on your cluster?

  a) Module system — what's the module name? (e.g. xtb/6.6)
  b) Conda environment — what's the env name?
  c) Installed at a specific path — what's the path?
     (I'll add it to $PATH in the script)
  d) Already in $PATH — no setup needed
```

If the user gives a path (option c), add it to `$PATH` in the script:
```bash
export PATH="/path/to/xtb/bin:$PATH"
```

For ORCA, also set `$LD_LIBRARY_PATH` if the user provides a path:
```bash
export PATH="/path/to/orca:$PATH"
export LD_LIBRARY_PATH="/path/to/orca:$LD_LIBRARY_PATH"
```

### Question 5 — GPU (only if relevant)

Skip this question unless the calculation could benefit from GPU acceleration or the user mentioned a GPU partition.

```
Do you need a GPU? If so, what type and how many?
(e.g. gpu:1, gpu:a100:2)
```

### Then generate

After collecting the answers, write the SLURM script using the template below. Use sensible defaults for anything not explicitly asked (1 node, 1 task per node, output to `%j.out`, error to same file, working directory = submission directory).

## Script Template

```bash
#!/bin/bash
#SBATCH --job-name=FILL_JOB_NAME
#SBATCH --output=%j.out
#SBATCH --partition=FILL_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
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

- Always use `$SLURM_CPUS_PER_TASK` for `OMP_NUM_THREADS` rather than hardcoding.
- ORCA requires `--ntasks-per-node` to match the `%pal nprocs` setting in the ORCA input. If generating both, keep them in sync.
- XTB reads `OMP_NUM_THREADS` for parallelism. No MPI needed.
- For ORCA parallel runs, the run command is `orca input.inp` (ORCA handles MPI internally).
- If the user doesn't know their partition or modules, suggest they run `sinfo -s` (partitions) and `module avail xtb` or `module avail orca` on the cluster.
