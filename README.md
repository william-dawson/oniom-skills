# oniom-skills

A Claude Code plugin for setting up QM and QM/MM calculations from protein-ligand PDB files.
The protein file you are working on should be in CHARMM format (generate with CHARMM-GUI)

## Installation

### Claude Code (CLI)

```
/plugin marketplace add william-dawson/oniom-skills
/plugin marketplace update oniom-skills
/plugin install pymol@oniom-skills
/reload-plugins
```


## What it does

This plugin provides three skills that guide you through the full workflow from raw PDB to a running calculation on your cluster:

1. **`/pymol:pymol`** — Extract a QM cluster model. Normalizes CHARMM PDBs, identifies the ligand, selects residues within a cutoff, applies amide capping, handles metal coordination shells, and determines the cluster charge.

2. **`/pymol:oniom`** — Set up a two-layer ONIOM calculation. Three approaches: XTB built-in ONIOM, 3-point manual (three separate calculations you control), or ORCA native QM/XTB. Walks you through method, cutoff, charge, and multiplicity one question at a time.

3. **`/pymol:slurm`** — Generate a SLURM batch script. Asks about partition, account, CPUs, memory, wall time, and how your software is installed (module, conda, or custom path). Produces a ready-to-submit script.

The skills are conversational — they inspect your structure, ask about ambiguous chemistry (ligand charge, protonation states, metal oxidation states, open-shell character), and write scripts tailored to your specific system.

## Requirements

- [PyMOL](https://pymol.org) (open-source or commercial) — must be on `$PATH`
- [XTB](https://github.com/grimme-lab/xtb) ≥ 6.6 — for ONIOM calculations with the XTB driver
- [ORCA](https://orcasoftware.de) — optional, for DFT-level inner region calculations
- [CHARMM-GUI](https://www.charmm-gui.org) - recommended way to generate the PDB files.

## Usage

Say something like:
```
Using your pymol skills, prepare an oniom calculation with XTB on the pdb file in this directory
```

## Codex

You can install these skills directly into Codex from GitHub.

For a project-local install:

```bash
mkdir -p .codex/skills
python3 ~/.codex/skills/.system/skill-installer/scripts/install-skill-from-github.py \
    --repo william-dawson/oniom-skills \
    --path skills/pymol skills/oniom skills/slurm \
    --dest "$PWD/.codex/skills"
```

## Credits

PyMOL is developed by [Schrödinger, LLC](https://pymol.org). XTB is developed by the [Grimme group](https://github.com/grimme-lab/xtb). ORCA is developed by the [Neese group](https://orcasoftware.de). This plugin provides only the Claude Code skill layer on top of these tools.
