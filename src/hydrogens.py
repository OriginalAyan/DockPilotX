"""Add polar / non-polar hydrogens to a receptor model in ChimeraX.

Polar hydrogens sit on N, O, S (and other non-carbon atoms); non-polar ones sit on carbon. ``addh`` is
restricted to one kind of heavy atom ("only hydrogens bonded to the specified part will be added"), and the
result is checked afterwards: a hydrogen of the wrong kind that still got added is deleted again, so a polar
hydrogen you removed by hand is never brought back by the non-polar step.
"""

from chimerax.core.commands import run
from chimerax.core.errors import UserError


def _spec(model):
    return "#" + ".".join(str(i) for i in model.id)


def _is_h(atom):
    try:
        return atom.element.name == "H"
    except Exception:
        return False


def _on_carbon(h):
    try:
        return any(n.element.name == "C" for n in h.neighbors)
    except Exception:
        return False


def count_hydrogens(model):
    """(polar, non-polar) hydrogen counts of an atomic model."""
    polar = nonpolar = 0
    try:
        for a in model.atoms:
            if _is_h(a):
                if _on_carbon(a):
                    nonpolar += 1
                else:
                    polar += 1
    except Exception:
        pass
    return polar, nonpolar


def _delete_atoms(model, atoms):
    if not atoms:
        return 0
    try:
        from chimerax.atomic import Atoms
        Atoms(atoms).delete()
    except Exception:
        for a in atoms:
            try:
                a.delete()
            except Exception:
                pass
    return len(atoms)


def _added_wrong_kind(model, before, want_nonpolar):
    """Hydrogens added since *before* that are not of the wanted kind."""
    bad = []
    for a in model.atoms:
        if a in before or not _is_h(a):
            continue
        if _on_carbon(a) != want_nonpolar:
            bad.append(a)
    return bad


def add_hydrogens(session, model, kind, log=None):
    """kind: 'polar' (hydrogens on non-carbon atoms) or 'nonpolar' (on carbon). Returns (added, removed_wrong)."""
    say = log or (lambda m: None)
    nonpolar = kind == "nonpolar"
    spec = _spec(model)
    before = set(model.atoms)
    n0 = len(before)
    part = "@@element_name=C" if nonpolar else "~@@element_name=C"
    try:
        run(session, "addh %s & %s" % (spec, part))
    except UserError as err:
        say("addh with an atom restriction was not accepted (%s); adding all hydrogens and removing the other kind." % err)
        run(session, "addh %s" % spec)
    wrong = _delete_atoms(model, _added_wrong_kind(model, before, nonpolar))
    added = len(model.atoms) - n0
    return added, wrong


def remove_nonpolar(model):
    """Delete the hydrogens on carbon atoms. Returns how many were removed."""
    return _delete_atoms(model, [a for a in model.atoms if _is_h(a) and _on_carbon(a)])


def lacks_nonpolar(model):
    """True if the model has carbon atoms but no non-polar hydrogens (minimization needs them)."""
    polar, nonpolar = count_hydrogens(model)
    return nonpolar == 0
