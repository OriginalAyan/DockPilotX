"""Clean a PDB/mmCIF receptor in ChimeraX and save it as a PDB file ready for Meeko.

Runs on ChimeraX's main thread. The user's own models are never modified: the source is opened as a
new temporary model, cleaned, saved and closed again.
"""

import os

from chimerax.core.commands import run
from chimerax.core.errors import UserError


def _q(path):
    return '"%s"' % os.path.abspath(path).replace("\\", "/")


def clean_receptor(session, source, out_pdb, *, chains="", remove_water=True, remove_ions=True,
                   remove_ligands=True, add_h=True, log=None):
    """Open *source*, clean it and write *out_pdb*. Returns a one-line summary string."""
    say = log or (lambda m: None)
    top = lambda: [m for m in session.models.list() if len(m.id) == 1]
    before = set(top())
    run(session, "open %s" % _q(source))
    new = [m for m in top() if m not in before]
    if not new:
        raise UserError("ChimeraX could not open %s" % source)
    model = new[0]
    try:
        kids = model.child_models()
        target = kids[0] if kids and not hasattr(model, "atoms") else model
        spec = "#" + ".".join(str(i) for i in target.id)
        for extra in new[1:]:
            extra.delete()
        n0 = len(target.atoms) if hasattr(target, "atoms") else 0
        if chains.strip():
            letters = ",".join(c.strip() for c in chains.replace(";", ",").split(",") if c.strip())
            run(session, "delete %s & ~/%s" % (spec, letters))
            say("Kept chain(s) %s." % letters)
        if remove_water:
            run(session, "delete %s & solvent" % spec)
        if remove_ions:
            run(session, "delete %s & ions" % spec)
        if remove_ligands:
            keep = ["protein", "nucleic"] + ([] if remove_water else ["solvent"]) + ([] if remove_ions else ["ions"])
            run(session, "delete %s & ~(%s)" % (spec, "|".join(keep)))
        if add_h:
            try:
                run(session, "addh %s hbond true useHisName true" % spec)
            except UserError:
                run(session, "addh %s hbond true" % spec)
        try:
            run(session, "save %s %s relModel %s" % (_q(out_pdb), spec, spec))
        except UserError:
            run(session, "save %s %s" % (_q(out_pdb), spec))
        n1 = len(target.atoms) if hasattr(target, "atoms") else 0
    finally:
        try:
            session.models.close([m for m in new if not getattr(m, "deleted", False)])
        except Exception:
            pass
    return "Cleaned %s: %d -> %d atoms%s." % (os.path.basename(source), n0, n1, " (hydrogens added)" if add_h else "")
