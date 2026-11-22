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
    """Follows a ChimeraX minimization through its status-line messages.

    ``minimize`` runs in the background and reports on the status line ("Starting minimization",
    "step N: energy X", "Minimization complete"). Those messages are intercepted by wrapping the logger's
    ``status`` method for the duration of the run (the original is called as before and put back at the end).
    *on_energy(step, energy)* is called for every energy report and *on_done(message, last_energy)* once when
    the minimization has completed.
    """

    DONE = re.compile(r"minimi[sz]ation\s+(complete|completed|converged|finished)", re.I)
    ENERGY = re.compile(r"step\s+(\d+)[^\d+-]{0,6}energy[^\d+-]{0,6}(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", re.I)

    def __init__(self, session, on_done, on_energy=None):
        self.session, self.on_done, self.on_energy = session, on_done, on_energy
        self.energies = []
        self.started = False
        self.finished = False
        self._orig = None
        self._wrapper = None

    def start(self):
        logger = self.session.logger
        self._orig = logger.status
        orig, watcher = self._orig, self

        def status(msg, *args, **kw):
            try:
                watcher._seen(str(msg))
            except Exception:
                pass
            return orig(msg, *args, **kw)
        self._wrapper = status
        logger.status = status
        return True

    def stop(self):
        logger = self.session.logger
        if self._wrapper is not None:
            try:
                if logger.__dict__.get("status") is self._wrapper:
                    del logger.status          # back to the class's own method
                elif getattr(logger, "status", None) is self._wrapper:
                    logger.status = self._orig
            except Exception:
                pass
            self._wrapper = None

    def _seen(self, text):
        low = text.lower()
        if "starting minimization" in low:
            self.started = True
        m = self.ENERGY.search(text)
        if m:
            step, e = int(m.group(1)), float(m.group(2))
            self.energies.append((step, e))
            if self.on_energy:
                self.on_energy(step, e)
        if self.DONE.search(text) and not self.finished:
            self.finished = True
            self._defer(text)

    def _defer(self, text):
        """Finish outside the status call itself (the logger is still busy with it)."""
        try:
            from Qt.QtCore import QTimer
            QTimer.singleShot(0, self._finish)
        except Exception:
            self._finish()

    def _finish(self):
        self.stop()
        energy = self.energies[-1][1] if self.energies else None
        msg = "Minimization completed"
        if energy is not None:
            msg += " (last reported energy %.2f at step %d)" % (energy, self.energies[-1][0])
        self.session.logger.info(msg)
        try:
            self.on_done(msg, energy)
        except Exception:
            pass


def assign_charges(session, model, say=None):
    """Assign force-field types and charges before minimizing (standard residues get Amber ff14SB values;
    anything else gets Gasteiger charges, which is fast). Meeko computes its own charges later, so these do
    not reach the docking inputs."""
    run(session, "addcharge %s method gasteiger" % _spec(model))
    if say:
        say("Assigned force-field types and charges (addcharge).")


TEMPLATE_HINT = ("A residue could not be matched to the force field. This usually means non-polar hydrogens are "
                 "missing (add them in step 2) or the model still has a ligand, ion or non-standard residue.")


def minimize_model(session, model, steps=100, until_converged=False, on_done=None, on_energy=None, log=None):
    """Assign charges, then start ChimeraX energy minimization (it runs in the background).

    ``logEnergy true`` makes ChimeraX report energies in its Log (every 100 steps; shorter runs report none).
    Hydrogens, protonation and the structure are not changed by the dock-prep step (``dockPrep false``).
    Returns the watcher (call ``.stop()`` to stop following).
    """
    say = log or (lambda m: None)
    max_steps = 0 if until_converged else int(steps)
    try:
        assign_charges(session, model, say)
    except UserError as err:
        raise UserError("Charges could not be assigned: %s\n%s" % (err, TEMPLATE_HINT))
    watcher = None
    if on_done is not None:
        watcher = MinimizeWatcher(session, on_done, on_energy)
        watcher.start()
    base = "minimize %s dockPrep false maxSteps %d" % (_spec(model), max_steps)
    try:
        try:
            run(session, base + " logEnergy true")
        except UserError as err:
            if "AMBER/GAFF" in str(err) or "partial charges" in str(err):
                raise UserError("%s\n%s" % (err, TEMPLATE_HINT))
            run(session, base)           # older ChimeraX without logEnergy
    except UserError:
        if watcher:
            watcher.stop()
        raise
    return watcher
