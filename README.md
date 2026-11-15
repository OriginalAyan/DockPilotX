# VinaPilotX (working name) — AutoDock Vina in UCSF ChimeraX

Runs **AutoDock Vina** from ChimeraX: search box, single and batch docking, results table. Inputs must already be prepared `.pdbqt` files (Meeko-based preparation is a later step).

**Documentation: [DETAILS.md](DETAILS.md)** — highlights, how every feature works, and the command list.

## Install / build (inside ChimeraX)
```
devel clean   <path-to-this-folder>
devel install <path-to-this-folder>     # then Tools > Binding Analysis > Lattice VinaPilotX, or:  vinax gui
devel build   <path-to-this-folder>     # wheel for the Toolshed
```

## Tests (outside ChimeraX)
```
python -m pytest -q tests               # core logic, uses a fake Vina
python tests/gui_smoke.py               # offscreen GUI run with stub ChimeraX (needs PySide6-Essentials)
```

## Licence
LGPL-2.1-only (derived from VinaGridBox v2.0.0, which adapts LGPL code from ChimeraX-Vinardock). See `LICENSE.txt` and `THIRD_PARTY_NOTICES.md`.
