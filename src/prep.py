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


def _res_key(line):
    """'A:1170' style key of a PDB ATOM/HETATM line (insertion code appended when present)."""
    return "%s:%s%s" % (line[21].strip(), line[22:26].strip(), line[26].strip())


def residue_atoms(pdb_path, keys):
    """{key: [atom names]} for the residues *keys* of a PDB file (for diagnostics)."""
    out = {k: [] for k in keys}
    try:
        with open(pdb_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 27:
                    k = _res_key(line)
                    if k in out:
                        out[k].append(line[12:16].strip())
    except OSError:
        pass
    return out


def strip_hydrogens(pdb_in, pdb_out, keys):
    """Copy a PDB file, dropping the hydrogens of the residues *keys*. Returns the number removed."""
    keys, removed = set(keys), 0
    with open(pdb_in, encoding="utf-8", errors="replace") as f, open(pdb_out, "w", encoding="utf-8") as g:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 27 and _res_key(line) in keys:
                elem = line[76:78].strip() if len(line) >= 78 else ""
                name = line[12:16].strip()
                if elem == "H" or (not elem and name.startswith("H")):
                    removed += 1
                    continue
            g.write(line)
    return removed


def _run_receptor(env, info, help_text, pdb, base, delete_bad, log, cancelled):
    cmd = tool_argv(env, info, "mk_prepare_receptor") + receptor_args(help_text, info, pdb, base, delete_bad)
    for old in glob.glob(base + "*.pdbqt"):
        os.remove(old)
    skipped = []

    def on_line(line):
        if "Template matching failed for" in line:
            skipped.extend(re.findall(r"'([^']+)'", line.split("for:", 1)[-1]))
        if log:
            log(line)

    code = run_stream(cmd, on_line, cancelled, cwd=os.path.dirname(base) or None, env=env.environ())
    if code == -1:
        raise PrepError("Cancelled.")
    produced = sorted(glob.glob(base + "*.pdbqt"))
    if code != 0 or not produced:
        raise PrepError("mk_prepare_receptor failed (exit code %d). See the log above." % code)
    exact = base + ".pdbqt"
    return (exact if exact in produced else produced[0]), skipped, cmd


def prepare_receptor(env, info, clean_pdb, out_base, delete_bad=True, log=None, cancelled=None, rescue=True):
    """Run mk_prepare_receptor on a cleaned PDB. Returns (path of the .pdbqt, [skipped residue keys], [rescued keys]).

    Residues Meeko cannot match are dropped (*delete_bad*). With *rescue*, such residues get a second try
    without their hydrogens (Meeko then adds standard ones from its templates), which fixes residues whose
    added hydrogens do not fit any template.
    """
    if info.get("prody") or info.get("prody_installed") is None:
        # a ProDy that cannot be imported (e.g. no pkg_resources) is repaired, else the plain PDB reader is used
        info = dict(prep_env.probe(env, repair=True, log=log, cancelled=cancelled).get("info") or info)
    help_text = tool_help(env, info, "mk_prepare_receptor")
    os.makedirs(os.path.dirname(out_base) or ".", exist_ok=True)
    path, skipped, cmd = _run_receptor(env, info, help_text, clean_pdb, out_base, delete_bad, log, cancelled)
    if log:
        log("mk_prepare_receptor: " + " ".join(_quote(c) for c in cmd[-8:]))
    rescued = []
    if skipped and rescue:
        if log:
            for k, names in residue_atoms(clean_pdb, skipped).items():
                log("Residue %s did not match a Meeko template; its atoms: %s" % (k, " ".join(names)))
            log("Retrying %d residue(s) without their hydrogens ..." % len(skipped))
        retry_pdb = out_base + "_retry.pdb"
        strip_hydrogens(clean_pdb, retry_pdb, skipped)
        retry_base = out_base + "_retry"
        try:
            path2, skipped2, _ = _run_receptor(env, info, help_text, retry_pdb, retry_base, delete_bad, log, cancelled)
        except PrepError as err:
            if str(err) == "Cancelled.":
                raise
            path2, skipped2 = None, skipped
        if path2 is not None and len(skipped2) < len(skipped):
            os.replace(path2, out_base + ".pdbqt")
            path = out_base + ".pdbqt"
            rescued = [k for k in skipped if k not in skipped2]
            skipped = skipped2
        for f in glob.glob(retry_base + "*") + [retry_pdb]:
            try:
                os.remove(f)
            except OSError:
                pass
    return path, skipped, rescued


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
    code = run_stream(cmd, on_line, cancelled, env=env.environ())
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
        code = run_stream(cmd, lambda l: log(l) if log else None, cancelled, env=env.environ())
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
