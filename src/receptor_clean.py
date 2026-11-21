"""Clean a PDB/mmCIF receptor in ChimeraX (Dock Prep) and save it as a PDB file ready for Meeko.

Runs on ChimeraX's main thread. The user's own models are never modified: the source is opened as a
new model ("<name>_prep"), cleaned with ``dockprep`` (plus ``delete`` for chains and ligands) and kept in
the scene so it can be minimized and inspected before it is saved and passed to Meeko.
"""

import os
import re

from chimerax.core.commands import run
from chimerax.core.errors import UserError

# Dock Prep's own wording for the incomplete side chain methods
SIDE_CHAIN_CHOICES = [
    ("Replace using Dunbrack rotamer library", "Dunbrack"),
    ("Replace using Dynameomics rotamer library", "Dynameomics"),
    ("Replace using Richardson common rotamer library", "Richardson.common"),
    ("Replace using Richardson most probable rotamer library", "Richardson.mode"),
    ("Change to alanine", "ala"),
    ("Change to glycine", "gly"),
]
STANDARDIZE_NOTE = [
    "selenomethionine (MSE) \u2192 methionine (MET)",
    "bromo-UMP (5BU) \u2192 UMP (U)",
    "methylselenyl-dUMP (UMS) \u2192 UMP (U)",
    "methylselenyl-dCMP (CSL) \u2192 CMP (C)",
]


def _q(path):
    return '"%s"' % os.path.abspath(path).replace("\\", "/")


def _spec(model):
    return "#" + ".".join(str(i) for i in model.id)


def _b(v):
    return "true" if v else "false"


def count_hydrogens(model):
    """Number of hydrogen atoms in an atomic model (0 if it cannot be told)."""
    atoms = getattr(model, "atoms", None)
    if atoms is None:
        return 0
    try:
        return int(sum(1 for e in atoms.element_names if e == "H"))
    except Exception:
        pass
    try:
        return sum(1 for a in atoms if getattr(getattr(a, "element", None), "name", "") == "H")
    except Exception:
        return 0


def _dockprep(session, spec, base, optional, say):
    """Run dockprep. *optional* keywords are dropped (with a note) if this ChimeraX rejects them."""
    last = None
    for drop in (False, True):
        opts = list(base) + ([] if drop else list(optional))
        try:
            run(session, "dockprep %s %s" % (spec, " ".join(opts)))
            if drop and optional:
                say("This ChimeraX did not accept %s; Dock Prep's defaults were used for them." % " / ".join(optional))
            return
        except UserError as err:
            last = err
            if not optional:
                break
    raise last


def clean_receptor(session, source, name, *, chains="", del_water=True, del_ions=True, del_ligands=True,
                   del_altlocs=True, del_incomplete=False, standardize=True,
                   side_chains="Dunbrack", log=None):
    """Open *source* as a new model named ``<name>_prep`` and clean it. Returns (model, summary).

    side_chains: a rotamer method / 'ala' / 'gly', or '' to leave incomplete side chains alone.
    Hydrogens are not touched here (see hydrogens.py and the Hydrogens step).
    """
    say = log or (lambda m: None)
    top = lambda: [m for m in session.models.list() if len(m.id) == 1]
    before = set(top())
    run(session, "open %s" % _q(source))
    new = [m for m in top() if m not in before]
    if not new:
        raise UserError("ChimeraX could not open %s" % source)
    model = new[0]
    ok = False
    try:
        kids = model.child_models()
        target = kids[0] if kids and not hasattr(model, "atoms") else model
        for extra in new[1:]:
            extra.delete()
        spec = _spec(target)
        n0 = len(target.atoms) if hasattr(target, "atoms") else 0
        if chains.strip():
            letters = ",".join(c.strip() for c in chains.replace(";", ",").split(",") if c.strip())
            run(session, "delete %s & ~/%s" % (spec, letters))
            say("Kept chain(s) %s." % letters)
        if del_ligands:
            run(session, "delete %s & ~(protein|nucleic|solvent|ions)" % spec)
        base = ["delSolvent %s" % _b(del_water), "delIons %s" % _b(del_ions),
                "delAltLocs %s" % _b(del_altlocs), "delMissingBackbone %s" % _b(del_incomplete),
                "ah false", "ac false"]
        optional = ["completeSideChains %s" % (side_chains if side_chains else "false")]
        if not standardize:
            optional.append("standardizeResidues none")
        _dockprep(session, spec, base, optional, say)
        try:
            target.name = "%s_prep" % name
        except Exception:
            pass
        n1 = len(target.atoms) if hasattr(target, "atoms") else 0
        ok = True
        return target, "Cleaned %s: %d -> %d atoms." % (os.path.basename(source), n0, n1)
    finally:
        if not ok:
            try:
                session.models.close([m for m in new if not getattr(m, "deleted", False)])
            except Exception:
                pass


def open_prep_model(session, path, name):
    """Open *path* as a new model named ``<name>_prep`` (used for the PDB2PQR result)."""
    top = lambda: [m for m in session.models.list() if len(m.id) == 1]
    before = set(top())
    run(session, "open %s" % _q(path))
    new = [m for m in top() if m not in before]
    if not new:
        raise UserError("ChimeraX could not open %s" % path)
    model = new[0]
    kids = model.child_models()
    target = kids[0] if kids and not hasattr(model, "atoms") else model
    for extra in new[1:]:
        extra.delete()
    try:
        target.name = "%s_prep" % name
    except Exception:
        pass
    return target


def save_model_pdb(session, model, out_pdb):
    spec = _spec(model)
    try:
        run(session, "save %s %s relModel %s" % (_q(out_pdb), spec, spec))
    except UserError:
        run(session, "save %s %s" % (_q(out_pdb), spec))


class MinimizeWatcher:
    """Listens to ChimeraX's status line to tell when a minimization has finished.

    ``minimize`` runs in the background and only reports on the status line (transient). This hook
    notices "Minimization ... completed", writes a permanent line (with the last energy) to the Log and
    calls *on_done*. It removes itself afterwards, or when stopped.
    """

    _ENERGY = re.compile(r"energy[^\d+-]{0,12}(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", re.I)

    def __init__(self, session, on_done, steps_note=""):
        self.session, self.on_done, self.steps_note = session, on_done, steps_note
        self.energies = []
        self.finished = False
        self._added = False
        self.excludes_log = False

    # -- logger interface ------------------------------------------------
    def log(self, level, msg, *args, **kw):
        return False

    def status(self, msg, color="black", secondary=False):
        try:
            text = str(msg)
            low = text.lower()
            if "minimi" not in low and "energy" not in low:
                return
            m = self._ENERGY.search(text)
            if m:
                self.energies.append(float(m.group(1)))
            if "minimi" in low and ("complete" in low or "converged" in low or "finished" in low):
                self._defer_finish(text)
        except Exception:
            pass

    def _defer_finish(self, text):
        """Finish outside the logger's own loop (removing a log while it iterates would skip entries)."""
        if self.finished:
            return
        try:
            from Qt.QtCore import QTimer
            QTimer.singleShot(0, lambda: self._finish(text))
        except Exception:
            self._finish(text)

    def start(self):
        try:
            self.session.logger.add_log(self)
            self._added = True
        except Exception:
            self._added = False
        return self._added

    def stop(self):
        if self._added:
            try:
                self.session.logger.remove_log(self)
            except Exception:
                pass
            self._added = False

    def _finish(self, text):
        if self.finished:
            return
        self.finished = True
        self.stop()
        energy = self.energies[-1] if self.energies else None
        msg = "Minimization completed"
        if energy is not None:
            msg += " (final energy %.2f)" % energy
        if self.steps_note:
            msg += ". " + self.steps_note
        self.session.logger.info(msg)
        try:
            self.on_done(msg, energy)
        except Exception:
            pass


def minimize_model(session, model, steps=100, until_converged=False, on_done=None):
    """Start ChimeraX energy minimization (runs in the background).

    Energies are written to the Log (``logEnergy true``) and, when it finishes, "Minimization completed" is
    written there too and *on_done(message, last_energy)* is called. Dock Prep is not repeated, so
    hydrogens, protonation and charges stay as they are. Returns the watcher (call ``.stop()`` to cancel).
    """
    max_steps = 0 if until_converged else int(steps)
    watcher = None
    if on_done is not None:
        watcher = MinimizeWatcher(session, on_done)
        watcher.start()
    base = "minimize %s dockPrep false maxSteps %d" % (_spec(model), max_steps)
    try:
        try:
            run(session, base + " logEnergy true")
        except UserError:           # older ChimeraX without logEnergy
            run(session, base)
    except UserError:
        if watcher:
            watcher.stop()
        raise
    return watcher
