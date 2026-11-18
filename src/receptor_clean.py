"""Clean a PDB/mmCIF receptor in ChimeraX (Dock Prep) and save it as a PDB file ready for Meeko.

Runs on ChimeraX's main thread. The user's own models are never modified: the source is opened as a
new model ("<name>_prep"), cleaned with ``dockprep`` (plus ``delete`` for chains and ligands) and kept in
the scene so it can be minimized and inspected before it is saved and passed to Meeko.
"""

import os

from chimerax.core.commands import run
from chimerax.core.errors import UserError

HYDROGEN_CHOICES = [
    ("Add if missing (recommended)", "auto"),
    ("Always add missing hydrogens", "always"),
    ("Don't add", "never"),
]
SIDE_CHAIN_CHOICES = [
    ("Complete with Dunbrack rotamers", "Dunbrack"),
    ("Complete with Dynameomics rotamers", "Dynameomics"),
    ("Complete with Richardson (common)", "Richardson.common"),
    ("Complete with Richardson (most probable)", "Richardson.mode"),
    ("Truncate to alanine", "ala"),
    ("Truncate to glycine", "gly"),
    ("Leave as they are", "none"),
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


def _dockprep(session, spec, options, side_chains, say):
    """Run dockprep. The side-chain keyword is the least certain one, so it is retried on failure."""
    sc_options = []
    if side_chains in ("none", "false", ""):
        sc_options = ["completeSideChains false"]
    elif side_chains:
        sc_options = ["completeSideChains %s" % side_chains]
    candidates = [sc_options] if sc_options else [[]]
    if side_chains in ("none", "false", ""):
        candidates.append(["completeSideChains none"])
    candidates.append([])
    last = None
    for sc in candidates:
        try:
            run(session, "dockprep %s %s" % (spec, " ".join(options + sc)))
            if sc != candidates[0] and sc == []:
                say("This ChimeraX did not accept the side-chain option; Dock Prep's default was used.")
            return
        except UserError as err:
            last = err
            if not sc:
                break
    raise last


def clean_receptor(session, source, name, *, chains="", del_water=True, del_ions=True, del_ligands=True,
                   del_altlocs=True, side_chains="Dunbrack", hydrogens="auto", log=None):
    """Open *source* as a new model named ``<name>_prep`` and clean it. Returns (model, summary).

    hydrogens: 'auto' = add only if the structure has none (keeps protonation you already set),
    'always' = add missing hydrogens, 'never' = leave as is.
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
        n_h = count_hydrogens(target)
        add_h = hydrogens == "always" or (hydrogens == "auto" and n_h == 0)
        options = ["delSolvent %s" % _b(del_water), "delIons %s" % _b(del_ions),
                   "delAltLocs %s" % _b(del_altlocs), "ah %s" % _b(add_h), "ac false"]
        _dockprep(session, spec, options, side_chains, say)
        try:
            target.name = "%s_prep" % name
        except Exception:
            pass
        n1 = len(target.atoms) if hasattr(target, "atoms") else 0
        if hydrogens == "never":
            h_note = "hydrogens left as they are (%d present)" % n_h
        elif add_h:
            h_note = "hydrogens added"
        else:
            h_note = "existing hydrogens kept (%d), none added" % n_h
        ok = True
        return target, "Cleaned %s: %d -> %d atoms; %s." % (os.path.basename(source), n0, n1, h_note)
    finally:
        if not ok:
            try:
                session.models.close([m for m in new if not getattr(m, "deleted", False)])
            except Exception:
                pass


def save_model_pdb(session, model, out_pdb):
    spec = _spec(model)
    try:
        run(session, "save %s %s relModel %s" % (_q(out_pdb), spec, spec))
    except UserError:
        run(session, "save %s %s" % (_q(out_pdb), spec))


def minimize_model(session, model, steps=100, until_converged=False):
    """Start ChimeraX energy minimization (runs in the background; ChimeraX shows live updates).

    Dock Prep is not repeated, so hydrogens, protonation and charges stay as they are. Standard residues
    use the force field's own charges.
    """
    max_steps = 0 if until_converged else int(steps)
    run(session, "minimize %s dockPrep false maxSteps %d" % (_spec(model), max_steps))
