"""Add polar / non-polar hydrogens to a receptor model in ChimeraX.

Polar hydrogens sit on N, O, S (and other non-carbon atoms); non-polar ones sit on carbon. ``addh`` cannot be
restricted to one kind of atom, so the new hydrogens are sorted afterwards: those of the wrong kind are deleted
again, so a polar hydrogen you removed by hand is never brought back by the non-polar step.
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
    """kind: 'polar' (hydrogens on non-carbon atoms) or 'nonpolar' (on carbon).

    ChimeraX ``addh`` takes whole structures only (an atom restriction is refused), so it always places every
    missing hydrogen and the result is sorted afterwards: hydrogens of the other kind are deleted again.
    Returns {'placed': all hydrogens ChimeraX added, 'kept': of the wanted kind, 'removed': other kind}.
    """
    say = log or (lambda m: None)
    nonpolar = kind == "nonpolar"
    before = set(model.atoms)
    run(session, "addh %s" % _spec(model))
    placed = sum(1 for a in model.atoms if a not in before and _is_h(a))
    removed = _delete_atoms(model, _added_wrong_kind(model, before, nonpolar))
    kept = placed - removed
    say("ChimeraX placed %d hydrogens (addh cannot be limited to one kind); kept %d %s; removed %d %s (not asked for here)." % (
        placed, kept, "on carbon" if nonpolar else "on N/O/S", removed, "on N/O/S" if nonpolar else "on carbon"))
    return {"placed": placed, "kept": kept, "removed": removed}


def add_hydrogens_to_atoms(session, model, targets):
    """Let ChimeraX protonate *targets* only: run addh, then delete every new hydrogen that is not bonded to
    one of them. Returns the number of hydrogens kept."""
    before = set(model.atoms)
    tset = set(targets)
    run(session, "addh %s" % _spec(model))
    bad = []
    for a in model.atoms:
        if a in before or not _is_h(a):
            continue
        try:
            ok = any(n in tset for n in a.neighbors)
        except Exception:
            ok = False
        if not ok:
            bad.append(a)
    _delete_atoms(model, bad)
    return sum(1 for t in targets for n in t.neighbors if _is_h(n))


def remove_nonpolar(model):
    """Delete the hydrogens on carbon atoms. Returns how many were removed."""
    return _delete_atoms(model, [a for a in model.atoms if _is_h(a) and _on_carbon(a)])


def lacks_nonpolar(model):
    """True if the model has carbon atoms but no non-polar hydrogens (minimization needs them)."""
    polar, nonpolar = count_hydrogens(model)
    return nonpolar == 0


# ---------------------------------------------------------------------------------------- chain ends
def _bonded(res_a, res_b):
    try:
        c, n = res_a.find_atom("C"), res_b.find_atom("N")
        return c is not None and n is not None and n in c.neighbors
    except Exception:
        return False


def chain_ends(model):
    """Chain ends of the protein chains of *model* as dicts {chain, kind 'N'|'C', number, real}.

    An end is *real* only at the true start/end of the full sequence. Ends next to missing residues (known
    from SEQRES, which ChimeraX reads) or next to a break in the chain are not real: they are artificial
    cuts, like the ends of a chain taken from a larger protein. Without SEQRES nothing is known about missing
    ends, so the ends of the resolved chain count as real (ChimeraX does the same).
    """
    ends = []
    try:
        chains = list(model.chains)
    except Exception:
        return ends
    for ch in chains:
        try:
            rs = list(ch.residues)
            cid = ch.chain_id
        except Exception:
            continue
        n, i = len(rs), 0
        while i < n:
            if rs[i] is None:
                i += 1
                continue
            j = i
            while j + 1 < n and rs[j + 1] is not None and _bonded(rs[j], rs[j + 1]):
                j += 1
            ends.append({"chain": cid, "kind": "N", "number": rs[i].number, "real": i == 0})
            ends.append({"chain": cid, "kind": "C", "number": rs[j].number, "real": j == n - 1})
            i = j + 1
    return ends


def cooh_hydrogen_position(o, c, oxt, length=0.97, angle=107.0):
    """Position for the acid hydrogen on OXT of a protonated C-terminus: cis to the carbonyl oxygen
    (dihedral O=C-OXT-H = 0). *o*, *c*, *oxt* are xyz triples."""
    import numpy as np
    a, b, cc = (np.asarray(x, dtype=float) for x in (o, c, oxt))
    bc = cc - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    th = np.radians(angle)
    d = cc + length * (-np.cos(th) * bc - np.sin(th) * m)      # cis to the carbonyl oxygen
    return d


def _dihedral(p0, p1, p2, p3):
    import numpy as np
    b0, b1, b2 = (np.asarray(p1) - np.asarray(p0)), (np.asarray(p2) - np.asarray(p1)), (np.asarray(p3) - np.asarray(p2))
    b1n = b1 / np.linalg.norm(b1)
    v, w = b0 - np.dot(b0, b1n) * b1n, b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


def _coord(atom):
    return [float(x) for x in atom.coord]


def _keep_one_amide_h(model, nat, res):
    """Keep one hydrogen on a cut N-terminal nitrogen: the one whose H-N-CA angle is closest to 120 degrees (planar
    amide), named H. The others are deleted."""
    import numpy as np
    hs = sorted([a for a in nat.neighbors if _is_h(a)], key=lambda a: a.name)
    best = hs[0]
    try:
        ca = res.find_atom("CA")
        n = np.array(_coord(nat))
        u = np.array(_coord(ca)) - n

        def err(h):
            v = np.array(_coord(h)) - n
            c = float(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))
            return abs(np.degrees(np.arccos(max(-1.0, min(1.0, c)))) - 120.0)
        scored = [(err(h), h.name, h) for h in hs]
        best = min(scored, key=lambda t: (t[0] if t[0] == t[0] else 1e9, t[1]))[2]
    except Exception:
        pass
    _delete_atoms(model, [h for h in hs if h is not best])
    try:
        best.name = "H"
    except Exception:
        pass


def fix_termini(session, model, ends, rows, ph, log=None):
    """Set the charge state of the chain ends after PDB2PQR (which keeps every end charged).

    Real ends follow PROPKA: N-terminus NH3+ stays unless its pKa is below the pH (then NH2); C-terminus
    COO- stays unless its pKa is above the pH (then COOH). Artificial ends (chain cut from a larger protein,
    or next to a gap) become neutral like a chain-internal residue: the NH3+ hydrogens are replaced by
    ChimeraX's amide hydrogen, and OXT is removed. *rows* are the PROPKA table rows. Returns
    {(chain, number, 'N'|'C'): note}.
    """
    say = log or (lambda m: None)
    notes = {}
    by_key = {(r["chain"], r["number"], r["group"]): r for r in rows}
    resmap = {}
    try:
        for r in model.residues:
            resmap[(r.chain_id, r.number)] = r
    except Exception:
        return notes
    for e in ends:
        res = resmap.get((e["chain"], e["number"]))
        if res is None:
            continue
        key = (e["chain"], e["number"], e["kind"])
        label = "%s-terminus %s%d" % ("N" if e["kind"] == "N" else "C", e["chain"], e["number"])
        try:
            if e["kind"] == "N":
                nat = res.find_atom("N")
                if nat is None:
                    continue
                hs = sorted([a for a in nat.neighbors if _is_h(a)], key=lambda a: a.name)
                if not e["real"]:
                    _delete_atoms(model, hs)
                    kept = add_hydrogens_to_atoms(session, model, [nat])
                    if kept > 1:
                        _keep_one_amide_h(model, nat, res)
                        say("%s: ChimeraX placed %d hydrogens on this cut end; kept the one closest to a planar "
                            "amide position and named it H." % (label, kept))
                    notes[key] = "not a real end: neutral amide (ChimeraX)"
                else:
                    row = by_key.get((e["chain"], e["number"], "N+"))
                    if row and row["pka"] < ph and len(hs) >= 3:
                        _delete_atoms(model, hs[-1:])
                        notes[key] = "real end: neutral NH2 (pKa below pH)"
                    else:
                        notes[key] = "real end: charged NH3+ (PROPKA)"
            else:
                oxt = res.find_atom("OXT")
                if not e["real"]:
                    if oxt is not None:
                        _delete_atoms(model, [a for a in oxt.neighbors if _is_h(a)] + [oxt])
                    notes[key] = "not a real end: neutral (OXT removed)"
                else:
                    row = by_key.get((e["chain"], e["number"], "C-"))
                    if row and row["pka"] > ph and oxt is not None and not any(_is_h(a) for a in oxt.neighbors):
                        o, c = res.find_atom("O"), res.find_atom("C")
                        pos = cooh_hydrogen_position(_coord(o), _coord(c), _coord(oxt))
                        from chimerax.atomic.struct_edit import add_atom
                        add_atom("HO", "H", res, pos, bonded_to=oxt)
                        notes[key] = "real end: COOH (pKa above pH)"
                    else:
                        notes[key] = "real end: charged COO- (PROPKA)"
        except Exception as err:
            say("%s: could not adjust (%s)" % (label, err))
            continue
        say("%s: %s" % (label, notes[key]))
    return notes
