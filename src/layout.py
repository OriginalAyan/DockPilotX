"""Names of the folders and files that Lattice VinaPilotX manages inside a working directory.

Every folder the tool owns starts with ``vinax_`` and holds a short note, so you can tell them from your own
files. Folders made by older versions (``receptor``, ``ligands`` ...) are still read.
"""

import os

PROJECT_FILE = "project.vinapx"
NOTE_NAME = "_MANAGED_BY_VinaPilotX.txt"
NOTE_TEXT = ("This folder is managed by Lattice VinaPilotX.\n"
             "Please do not rename, move or edit the files in it: the project file (project.vinapx) and the\n"
             "Results tab rely on these names. Copy files out if you want to change them.\n")

NAMES = {"receptor": "vinax_receptor", "ligands": "vinax_ligands", "results": "vinax_results",
         "logs": "vinax_logs", "prepared": "vinax_prepared", "cache": "vinax_pose_cache",
         "selected": "vinax_selected_poses"}
LEGACY = {"receptor": "receptor", "ligands": "ligands", "results": "results", "logs": "logs",
          "prepared": "prepared", "cache": "pose_cache", "selected": "selected_poses"}


def path(workdir, key):
    """Where the folder *key* is written."""
    return os.path.join(workdir, NAMES[key])


def existing(workdir, key):
    """Where to read folder *key* from: the new name, else the name used by older versions."""
    new = path(workdir, key)
    old = os.path.join(workdir, LEGACY[key])
    if os.path.isdir(new) or not os.path.isdir(old):
        return new
    return old


def make(workdir, key):
    """Create the folder (with its note file) and return its path."""
    p = path(workdir, key)
    os.makedirs(p, exist_ok=True)
    note = os.path.join(p, NOTE_NAME)
    if not os.path.exists(note):
        try:
            with open(note, "w", encoding="utf-8") as f:
                f.write(NOTE_TEXT)
        except OSError:
            pass
    return p
