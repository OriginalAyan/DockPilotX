# Third-party notices

## Search-box editor and mouse interaction

The search-box editor (box model, undo/redo, panel) is derived from ChimeraX Vina Grid Box
(https://github.com/OriginalAyan/VinaGridBox, same author). Its mouse move/resize code in turn adapts
`mousemode.py` from ChimeraX-Vinardock (https://github.com/rquiroga7/ChimeraX-Vinardock) by Rodrigo
Quiroga, licensed LGPL-2.1. Keep this notice and `LICENSE.txt` with all redistributions.

## AutoDock Vina (not bundled)

AutoDock Vina (https://github.com/ccsb-scripps/AutoDock-Vina) is a separate program under the Apache-2.0
licence. Lattice VinaPilotX does not include it; it runs a copy the user supplies or downloads from the official
release page. If you publish results, cite Trott & Olson, J. Comput. Chem. 2010, 31, 455 and
Eberhardt et al., J. Chem. Inf. Model. 2021, 61, 3891.

## Meeko, RDKit, micromamba, ProDy (not bundled)

The Prepare tab installs, on request and into a private environment, Meeko (https://github.com/forlilab/Meeko,
LGPL-2.1), RDKit (BSD-3-Clause), ProDy (MIT), gemmi (MPL-2.0), SciPy/NumPy (BSD-3-Clause) and, for the
conda route, micromamba (BSD-3-Clause, https://github.com/mamba-org/mamba) with packages from conda-forge. None of
them is included in this bundle. Cite Meeko (Forli lab) and RDKit if you publish results prepared with them.
