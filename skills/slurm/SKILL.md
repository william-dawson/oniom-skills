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

## Stacksize

A lot of codes require you to increase the openmp stacksize:
```
export OMP_STACKSIZE=512M
```
XTB is a good example of this. Confirm with the user if you are not sure. 

## Notes
- Use `$SLURM_CPUS_PER_TASK` for thread counts, never hardcode.
- ORCA: keep `--ntasks-per-node` in sync with `%pal nprocs` in the input.
- XTB: reads `OMP_NUM_THREADS`, no MPI needed.

---

## Periodic Fetch and Visualization Workflow

For iterative metadynamics runs, the real workflow is **rsync + local analysis**, not just checking `get_job_status`.

### One-shot sync

```bash
mkdir -p outputs/my_project
cd outputs/my_project

# Pull only output files (not inputs) from all run_* directories
rsync -avz --prune-empty-dirs \
  --include='run_*/' \
  --include='run_*/xtbopt.log' \
  --include='run_*/xtbopt.xyz' \
  --include='run_*/*.log' \
  --include='run_*/.xtboptok' \
  --include='run_*/reference.xyz' \
  --exclude='*' \
  user@cluster:/path/to/my_project/ .
```

### Continuous poller

Save as `poll_and_fetch.py` and run locally:

```python
#!/usr/bin/env python3
"""Poll cluster for completed metadynamics runs and rsync results."""
import argparse, subprocess, sys, time

REMOTE_USER = "user"
REMOTE_HOST = "cluster"
REMOTE_BASE = "/path/to/my_project"
LOCAL_BASE = "./outputs/my_project"

def rsync_outputs():
    cmd = [
        "rsync", "-avz", "--prune-empty-dirs",
        "--include=run_*/",
        "--include=run_*/xtbopt.log",
        "--include=run_*/xtbopt.xyz",
        "--include=run_*/*.log",
        "--include=run_*/.xtboptok",
        "--include=run_*/reference.xyz",
        "--exclude=*",
        f"{REMOTE_USER}@{REMOTE_HOST}:{REMOTE_BASE}/",
        LOCAL_BASE,
    ]
    subprocess.run(cmd)

def count_completed():
    import glob, os
    runs = glob.glob(f"{LOCAL_BASE}/run_*")
    done = sum(1 for d in runs if os.path.isfile(f"{d}/.xtboptok"))
    return len(runs), done

def count_frames(run_dir):
    log = f"{run_dir}/xtbopt.log"
    if not os.path.isfile(log):
        return 0
    result = subprocess.run(["grep", "-c", "^\\s*\\d\\+\\s*$", log],
                            capture_output=True, text=True)
    return int(result.stdout.strip()) if result.returncode == 0 else 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--watch", action="store_true", help="loop forever")
    parser.add_argument("--interval", type=int, default=60)
    args = parser.parse_args()

    while True:
        rsync_outputs()
        total, done = count_completed()
        print(f"Runs: {done}/{total} completed")
        for d in sorted(glob.glob(f"{LOCAL_BASE}/run_*")):
            n = count_frames(d)
            status = "DONE" if os.path.isfile(f"{d}/.xtboptok") else "running"
            print(f"  {os.path.basename(d)}: {n} frames [{status}]")
        if not args.watch:
            break
        time.sleep(args.interval)
```

Run:
```bash
python3 poll_and_fetch.py          # one-shot
python3 poll_and_fetch.py --watch  # every 60s
```

### Load trajectory in PyMOL

```bash
pymol outputs/my_project/run_h1_saltbridge/xtbopt.log
```

Then in PyMOL:
```
load outputs/my_project/run_h1_saltbridge/reference.xyz, ref
align traj, ref
show sticks, ref and not elem H
color gray, ref
show sticks, traj and not elem H
color cyan, traj
```

### Key insight

`get_job_status` and `read_job_output` are for diagnostics, but the real feedback loop is:
1. `rsync` pulls `xtbopt.log` trajectories
2. PyMOL visualizes them locally
3. You decide which hypotheses worked and design the next batch

