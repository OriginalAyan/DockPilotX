"""Receptor and ligand preparation with Meeko, run in the preparation environment.

Pure Python. Receptor *cleaning* (waters, hydrogens, ...) is done in ChimeraX by the tool before
``prepare_receptor`` is called; this module only drives the external programs.
"""

import glob
import json
import os
import re

from . import prep_env
from .prep_env import PrepError, run_stream, tool_argv, tool_help

LIGAND_SUFFIXES = (".sdf", ".mol", ".mol2", ".smi", ".smiles")
_STAGE_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stage_ligands.py")


def list_ligand_files(paths, recursive=True):
    """Expand files and folders into ligand input files (sorted, unique)."""
    out, seen = [], set()
    for p in paths:
        p = os.path.abspath(os.path.expanduser(str(p).strip().strip('"')))
        found = []
        if os.path.isdir(p):
            for root, dirs, files in os.walk(p):
                found += [os.path.join(root, f) for f in sorted(files)]
                if not recursive:
                    break
        elif os.path.isfile(p):
            found = [p]
        for f in found:
            if f.lower().endswith(LIGAND_SUFFIXES) and f not in seen:
                seen.add(f)
                out.append(f)
    return out


# ---------------------------------------------------------------------------------------------
def receptor_args(help_text, info, clean_pdb, out_base, delete_bad=True):
    """Command-line options for mk_prepare_receptor, adapted to the installed Meeko version."""
    h = help_text or ""
    args = []
    if "--read_with_prody" in h:
        args += ["--read_with_prody" if info.get("prody") else "--read_pdb", clean_pdb]
        if not info.get("prody") and "--read_pdb" not in h:
            raise PrepError("This Meeko needs ProDy to read PDB files, but ProDy is not installed in the environment.")
    elif "--read_pdb" in h:
        args += ["--read_pdb", clean_pdb]
    else:
        args += ["-i", clean_pdb]
    args += ["-o", out_base, "-p"]
    if delete_bad:
        if "--delete_bad_res" in h:
            args.append("--delete_bad_res")
        elif "--allow_bad_res" in h:
            args.append("--allow_bad_res")
    return args


def prepare_receptor(env, info, clean_pdb, out_base, delete_bad=True, log=None, cancelled=None):
    """Run mk_prepare_receptor on a cleaned PDB. Returns (path of the written .pdbqt, [keys of residues Meeko skipped, e.g. 'A:7'])."""
    help_text = tool_help(env, info, "mk_prepare_receptor")
    cmd = tool_argv(env, info, "mk_prepare_receptor") + receptor_args(help_text, info, clean_pdb, out_base, delete_bad)
    if log:
        log("mk_prepare_receptor: " + " ".join(_quote(c) for c in cmd[-8:]))
    os.makedirs(os.path.dirname(out_base) or ".", exist_ok=True)
    for old in glob.glob(out_base + "*.pdbqt"):
        os.remove(old)
    warnings = []

    def on_line(line):
        if "Template matching failed for" in line:
            warnings.extend(re.findall(r"'([^']+)'", line.split("for:", 1)[-1]))
        if log:
            log(line)

    code = run_stream(cmd, on_line, cancelled, cwd=os.path.dirname(out_base) or None)
    if code == -1:
        raise PrepError("Cancelled.")
    produced = sorted(glob.glob(out_base + "*.pdbqt"))
    exact = out_base + ".pdbqt"
    if code != 0 or not produced:
        raise PrepError("mk_prepare_receptor failed (exit code %d). See the log above." % code)
    return (exact if exact in produced else produced[0]), warnings


def _quote(s):
    return '"%s"' % s if " " in s else s


# ---------------------------------------------------------------------------------------------
def prepare_ligands(env, info, inputs, out_dir, stage_dir, add_h=True, log=None, progress=None, cancelled=None):
    """Prepare ligands: RDKit staging (H / 3D), then mk_prepare_ligand.

    Returns {'files': [pdbqt paths], 'failed': [(name, message)], 'n_inputs': int}.
    """
    files = list_ligand_files(inputs)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(stage_dir, exist_ok=True)
    result = {"files": [], "failed": [], "n_inputs": len(files)}
    if not files:
        raise PrepError("No ligand files found (supported: %s)." % ", ".join(LIGAND_SUFFIXES))
    stage_lines = []

    def on_line(line):
        if line.startswith("VPX"):
            try:
                stage_lines.append(json.loads(line[3:]))
            except ValueError:
                pass
        elif log and not line.startswith("VPXDONE"):
            log(line)

    cmd = env.python_command([_STAGE_SCRIPT, "--out", stage_dir] + (["--no-addh"] if not add_h else []) + files)
    if log:
        log("Reading %d ligand file(s) with RDKit ..." % len(files))
    code = run_stream(cmd, on_line, cancelled)
    if code == -1:
        raise PrepError("Cancelled.")
    staged = [d for d in stage_lines if "staged" in d and d["staged"]]
    for d in stage_lines:
        if "ok" in d and not d["ok"]:
            result["failed"].append((d.get("name", "?"), d.get("msg", "failed")))
            if log:
                log("  skipped %s: %s" % (d.get("name"), d.get("msg")))
    if code != 0 and not staged:
        raise PrepError("RDKit could not read the ligands (exit code %d)." % code)
    total = sum(d.get("count", 0) for d in staged)
    done = 0
    for d in staged:
        if cancelled is not None and cancelled():
            raise PrepError("Cancelled.")
        sdf, n = d["staged"], d.get("count", 0)
        before = set(glob.glob(os.path.join(out_dir, "*.pdbqt")))
        base = tool_argv(env, info, "mk_prepare_ligand")
        if n == 1:
            name = [x["name"] for x in stage_lines if x.get("ok") and x.get("file") == d["file"]][0]
            cmd = base + ["-i", sdf, "-o", os.path.join(out_dir, name + ".pdbqt")]
        else:
            cmd = base + ["-i", sdf, "--multimol_outdir", out_dir]
        if log:
            log("mk_prepare_ligand: %s (%d molecule%s)" % (os.path.basename(sdf), n, "" if n == 1 else "s"))
        code = run_stream(cmd, lambda l: log(l) if log else None, cancelled)
        if code == -1:
            raise PrepError("Cancelled.")
        new = sorted(set(glob.glob(os.path.join(out_dir, "*.pdbqt"))) - before)
        result["files"] += new
        if code != 0:
            result["failed"].append((os.path.basename(sdf), "mk_prepare_ligand exit code %d" % code))
        done += n
        if progress:
            progress(done, total)
    return result
