# oniom-skills

A Claude Code plugin for setting up QM and QM/MM calculations from protein-ligand PDB files.

## Installation

### Claude Code (CLI)

```
/plugin marketplace add william-dawson/oniom-skills
/plugin install pymol@oniom-skills
/reload-plugins
```

### Claude Desktop

1. Download `oniom-skills.zip` from the [latest release](https://github.com/william-dawson/oniom-skills/releases/latest)
2. Open the **Cowork** tab at the top of Claude Desktop
3. Select **Customize** on the left sidebar
4. Select **Skills**
5. Click the **+** button and upload both `.md` files from the zip

## What it does

This plugin provides two skills that guide you through the full workflow from raw PDB to a ready-to-run quantum chemistry calculation:

1. **`/pymol:pymol`** — Extract a QM cluster model. Selects residues around a ligand, applies amide capping (link-atom approach), and determines the cluster charge.

2. **`/pymol:oniom`** — Set up a two-layer ONIOM calculation. Supports XTB as the driver (GFN2/GFN-FF) and ORCA as the driver (QM/XTB embedding). Handles inner region selection, charge/multiplicity, and generates ready-to-run input files.

The skills are conversational — they inspect your structure, ask about ambiguous chemistry (ligand charge, histidine protonation states, metal oxidation states, open-shell character), and write a script tailored to your specific system.

## Requirements

- [PyMOL](https://pymol.org) (open-source or commercial) — must be on `$PATH`
- [XTB](https://github.com/grimme-lab/xtb) ≥ 6.6 — for ONIOM calculations with the XTB driver
- [ORCA](https://orcasoftware.de) — optional, for DFT-level inner region calculations

## Usage

```
/pymol:pymol    # start cluster extraction
/pymol:oniom    # set up ONIOM calculation
```

Both skills write Python scripts using PyMOL's bundled interpreter, discovered automatically at runtime — no hardcoded paths.

## Credits

PyMOL is developed by [Schrödinger, LLC](https://pymol.org). XTB is developed by the [Grimme group](https://github.com/grimme-lab/xtb). ORCA is developed by the [Neese group](https://orcasoftware.de). This plugin provides only the Claude Code skill layer on top of these tools.
