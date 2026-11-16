"""Runs inside the preparation environment (needs RDKit): read ligands, add hydrogens / 3D coordinates,
and write one clean SDF per input file for Meeko.

Usage: python stage_ligands.py --out STAGE_DIR [--no-addh] [--seed N] FILE [FILE ...]
Prints one line per molecule starting with VPX followed by JSON.
"""

import argparse
import json
import os
import re
import sys


def emit(**kw):
    print("VPX" + json.dumps(kw), flush=True)


def safe(name):
    s = re.sub(r"[^A-Za-z0-9._+-]+", "_", str(name).strip()).strip("._")
    return s[:80] or "ligand"


def read_mols(path):
    """Yield (name_or_None, mol_or_None, message)."""
    from rdkit import Chem
    ext = os.path.splitext(path)[1].lower()
    if ext == ".sdf":
        for m in Chem.SDMolSupplier(path, removeHs=False, sanitize=True):
            yield (m.GetProp("_Name") if m is not None and m.HasProp("_Name") else None), m, "could not be read"
    elif ext == ".mol":
        m = Chem.MolFromMolFile(path, removeHs=False)
        yield None, m, "could not be read"
    elif ext == ".mol2":
        m = Chem.MolFromMol2File(path, removeHs=False)
        yield None, m, "could not be read"
    elif ext in (".smi", ".smiles", ".txt"):
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split()
                if not parts or parts[0].startswith("#"):
                    continue
                yield (parts[1] if len(parts) > 1 else None), Chem.MolFromSmiles(parts[0]), "invalid SMILES " + parts[0]
    else:
        yield None, None, "unsupported file type %s (use .sdf, .mol, .mol2 or .smi)" % ext


def prepare(mol, add_h, seed):
    """Return (mol, note). Adds hydrogens when missing and 3D coordinates when absent or flat."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    note = []
    has_h = any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())
    conf = mol.GetConformer() if mol.GetNumConformers() else None
    need3d = conf is None or not conf.Is3D()
    if add_h and (not has_h or need3d):
        mol = Chem.AddHs(mol, addCoords=not need3d)
        note.append("added H")
    if need3d:
        mol.RemoveAllConformers()
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                raise RuntimeError("3D embedding failed")
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass
        note.append("3D embedded")
    return mol, ", ".join(note)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-addh", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    from rdkit import Chem
    os.makedirs(a.out, exist_ok=True)
    used = set()
    for n_file, path in enumerate(a.files, start=1):
        stem = safe(os.path.splitext(os.path.basename(path))[0])
        out = os.path.join(a.out, "%03d_%s.sdf" % (n_file, stem))
        writer, n_ok, idx = None, 0, 0
        try:
            for raw_name, mol, msg in read_mols(path):
                idx += 1
                label = raw_name or ("%s_%d" % (stem, idx) if idx > 1 else stem)
                name = safe(label)
                base, k = name, 2
                while name in used:
                    name = "%s_%d" % (base, k)
                    k += 1
                if mol is None:
                    emit(file=path, name=name, ok=False, msg=msg)
                    continue
                try:
                    mol, note = prepare(mol, not a.no_addh, a.seed)
                except Exception as err:
                    emit(file=path, name=name, ok=False, msg=str(err))
                    continue
                used.add(name)
                mol.SetProp("_Name", name)
                if writer is None:
                    writer = Chem.SDWriter(out)
                writer.write(mol)
                n_ok += 1
                emit(file=path, name=name, ok=True, msg=note)
        except Exception as err:
            emit(file=path, name=stem, ok=False, msg=str(err))
        if writer is not None:
            writer.close()
        emit(file=path, staged=out if n_ok else None, count=n_ok)
    print("VPXDONE", flush=True)


if __name__ == "__main__":
    sys.exit(main())
