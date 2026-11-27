import os, sys, threading, functools, http.server, pytest
import conftest
from ldsrc import prep_env as pe, prep

NEW_HELP = "usage: mk_prepare_receptor.py [--read_pdb X] [-i M] [--read_with_prody] [--delete_bad_res] [-d D]"
OLD_HELP = "usage: mk_prepare_receptor.py [-i RECEPTOR] [-a] [--allow_bad_res]"


def test_receptor_args_versions():
    a = prep.receptor_args(NEW_HELP, {"prody": True}, "c.pdb", "out/rec", True)
    assert a[:2] == ["--read_with_prody", "c.pdb"] and "--delete_bad_res" in a and a[2:5] == ["-o", "out/rec", "-p"]
    assert prep.receptor_args(NEW_HELP, {"prody": False}, "c.pdb", "o", False)[:1] == ["--read_pdb"]
    assert "--delete_bad_res" not in prep.receptor_args(NEW_HELP, {"prody": True}, "c.pdb", "o", False)
    old = prep.receptor_args(OLD_HELP, {}, "c.pdb", "o", True)
    assert old[:2] == ["-i", "c.pdb"] and "--allow_bad_res" in old


def test_probe_flags_unimportable_prody(tmp_path):
    """A ProDy that is installed but fails to import must count as unavailable (real import, not find_spec)."""
    site = tmp_path / "site"; (site / "prody").mkdir(parents=True)
    (site / "prody" / "__init__.py").write_text("import pkg_resources_not_there\n")
    py = pe.external_env(sys.executable)
    py.python_command = lambda args: [sys.executable] + list(args)
    env = pe.PrepEnv("external", python=sys.executable)
    os.environ["PYTHONPATH"] = str(site)
    try:
        info = pe._probe_once(env)["info"]
    finally:
        os.environ.pop("PYTHONPATH", None)
    assert info.get("prody") is False and info.get("prody_installed") is True and "pkg_resources_not_there" in info["prody_error"]


def test_setuptools_pinned():
    assert pe.SETUPTOOLS_PIN in pe.CONDA_PACKAGES and pe.SETUPTOOLS_PIN in pe.PIP_OPTIONAL


def test_list_ligand_files(tmp_path):
    (tmp_path / "sub").mkdir()
    for n in ("a.sdf", "b.MOL2", "c.txt", "sub/d.smi", "e.pdb"):
        (tmp_path / n).write_text("x")
    got = [os.path.basename(p) for p in prep.list_ligand_files([str(tmp_path), str(tmp_path / "a.sdf")])]
    assert sorted(got) == ["a.sdf", "b.MOL2", "d.smi"]
    assert [os.path.basename(p) for p in prep.list_ligand_files([str(tmp_path)], recursive=False)] == ["a.sdf", "b.MOL2"]


def test_env_commands_and_roundtrip(tmp_path):
    mm = pe.PrepEnv("micromamba", prefix="/e", mamba="/mm", root="/r")
    assert mm.command(["python", "-m", "x"])[1:] == ["-m", "x"] and "python" in mm.command(["python"])[0]
    assert mm.command(["tool.py"]) == ["/mm", "run", "-r", "/r", "-p", "/e", "tool.py"]
    assert mm.environ()["CONDA_PREFIX"] == "/e" and mm.environ()["PATH"].split(os.pathsep)[0] in ("/e", "/e/bin")
    assert pe.PrepEnv("external", python="/p").environ() is None
    v = pe.PrepEnv("venv", prefix="/v", python="/v/bin/python")
    assert v.python_command(["-c", "1"]) == ["/v/bin/python", "-c", "1"] and v.command(["tool.py"]) == ["tool.py"]
    assert pe.PrepEnv.from_dict(mm.to_dict()).to_dict() == mm.to_dict() and pe.PrepEnv.from_dict(None) is None
    ext = pe.external_env(sys.executable)
    assert ext.kind == "external" and ext.exists()
    assert pe.external_env(str(tmp_path / "nope")) is None
    assert pe.probe(None)["ok"] is False


def test_probe_reports_missing_meeko():
    r = pe.probe(pe.external_env(sys.executable))
    assert not r["ok"] and ("Meeko" in r["message"] or "RDKit" in r["message"]) or r["ok"]


def test_download_micromamba(tmp_path, monkeypatch):
    srv = tmp_path / "srv"; srv.mkdir()
    plat, ext = pe.mamba_platform()
    (srv / ("micromamba-%s%s" % (plat, ext))).write_bytes(b"\0" * 200_000)
    h = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(srv))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), h)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    monkeypatch.setattr(pe, "MICROMAMBA_URLS", ["http://127.0.0.1:%d/micromamba-{plat}{ext}" % httpd.server_port])
    try:
        path = pe.download_micromamba(str(tmp_path / "dl"))
    finally:
        httpd.shutdown()
    assert os.path.getsize(path) == 200_000 and os.access(path, os.X_OK) or ext


def test_find_python_for_venv():
    assert pe.find_python_for_venv()


REAL = os.environ.get("VPX_TEST_PYTHON")


@pytest.mark.skipif(not REAL, reason="set VPX_TEST_PYTHON to a Python that has Meeko + RDKit")
def test_real_meeko_pipeline(tmp_path):
    env = pe.external_env(REAL)
    res = pe.probe(env); assert res["ok"], res
    info = res["info"]
    gen = """
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
d = sys.argv[1]
m = Chem.AddHs(Chem.MolFromSequence("AGSDKHF"), addResidueInfo=True)
ps = AllChem.ETKDGv3(); ps.useRandomCoords = True; ps.randomSeed = 3; AllChem.EmbedMolecule(m, ps)
Chem.MolToPDBFile(m, d + "/rec.pdb")
w = Chem.SDWriter(d + "/lib.sdf")
for n, s in (("asp", "CC(=O)Oc1ccccc1C(=O)O"), ("caf", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")):
    mm = Chem.MolFromSmiles(s); mm.SetProp("_Name", n); w.write(mm)
w.close()
"""
    assert pe.run_stream(env.python_command(["-c", gen, str(tmp_path)])) == 0
    pdb = tmp_path / "rec.pdb"
    path, warns, rescued = prep.prepare_receptor(env, info, str(pdb), str(tmp_path / "out" / "rec"), True)
    assert path.endswith("rec.pdbqt") and os.path.getsize(path) > 100 and warns == [] and rescued == ["A:7"]
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "lib.sdf").write_text((tmp_path / "lib.sdf").read_text())
    (tmp_path / "in" / "x.smi").write_text("CCO ethanol\nXXX broken\n")
    (tmp_path / "in" / "same").mkdir(); (tmp_path / "in" / "same" / "lib.sdf").write_text((tmp_path / "in" / "lib.sdf").read_text())
    r = prep.prepare_ligands(env, info, [str(tmp_path / "in")], str(tmp_path / "L"), str(tmp_path / "S"), True)
    names = sorted(os.path.basename(f) for f in r["files"])
    assert "ethanol.pdbqt" in names and "asp.pdbqt" in names and "asp_2.pdbqt" in names, names   # same stem in two folders: no overwrite
    assert [b[0] for b in r["failed"]] == ["broken"]
    # the outputs are read by the docking input checks
    from ldsrc import pdbqt_tools as pt
    assert pt.check_receptor_pdbqt(path) is not None
    recs = pt.scan_ligands([str(tmp_path / "L")], False, set())
    assert len(recs) == len(r["files"]) and all(not x.problems for x in recs)


def test_strip_hydrogens_and_residue_atoms(tmp_path):
    def atom(n, name, res, resno, elem):
        return "ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s" % (n, name, res, resno, 1.0, 2.0, 3.0, elem)
    lines = [atom(1, "N", "HIS", 7, "N"), atom(2, "H", "HIS", 7, "H"), atom(3, "HE2", "HIS", 7, "H"), atom(4, "CA", "HIS", 7, "C"),
             atom(5, "N", "ALA", 8, "N"), atom(6, "H", "ALA", 8, "H")]
    src = tmp_path / "a.pdb"; src.write_text("\n".join(lines) + "\n")
    assert prep.residue_atoms(str(src), ["A:7"]) == {"A:7": ["N", "H", "HE2", "CA"]}
    out = tmp_path / "b.pdb"
    assert prep.strip_hydrogens(str(src), str(out), ["A:7"]) == 2
    kept = [l[12:16].strip() for l in out.read_text().splitlines()]
    assert kept == ["N", "CA", "N", "H"]      # residue 8 keeps its hydrogen


# ------------------------------------------------------------------------------ PDB2PQR
from ldsrc import pdb2pqr_tools as pq

SUMMARY = """junk
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP   4 A     3.94       3.80                      
   HIS   6 A     6.31       6.50                      
   C-   12 A     3.30       3.20                      
   N+    1 B     7.68       8.00                      

INFO:PROPKA information:
   LYS   9 A    10.00      10.50
"""


def test_pdb2pqr_summary_and_flags():
    rows = pq.parse_pka_summary(SUMMARY)
    assert [(r["group"], r["number"], r["chain"]) for r in rows] == [("ASP", 4, "A"), ("HIS", 6, "A"), ("C-", 12, "A"), ("N+", 1, "B")]
    assert pq.classify(6.4, 7.4) == "close" and pq.classify(8.4, 7.4) == "close" and pq.classify(6.39, 7.4) == "near"
    assert pq.classify(5.4, 7.4) == "near" and pq.classify(9.4, 7.4) == "near" and pq.classify(5.39, 7.4) == "" and pq.classify(9.5, 7.4) == ""
    t = pq.build_table(rows, {("A", 6): "HIS"}, {("A", 6): "HID"}, 7.4)
    assert t[1]["renamed"] == "HIS → HID" and t[3]["residue"].startswith("N-terminus") and t[2]["predicted"] == "deprotonated"



def test_pdb2pqr_args_and_hydrogen_strip(tmp_path):
    a = pq.pdb2pqr_args(7.4, "o.pdb", "o.pqr", "i.pdb")
    assert a[:4] == ["--ff", "AMBER", "--ffout", "AMBER"] and "--with-ph" in a and a[a.index("--with-ph") + 1] == "7.4" and a[-2:] == ["i.pdb", "o.pqr"]
    src = tmp_path / "a.pdb"
    src.write_text("ATOM      1  N   ALA A   1       0.0   0.0   0.0  1.00  0.00           N  \n"
                   "ATOM      2  H   ALA A   1       0.0   0.0   1.0  1.00  0.00           H  \n"
                   "ATOM      3 1HB  ALA A   1       0.0   0.0   2.0  1.00  0.00           H  \n"
                   "HETATM    4 HG    HG A   2       0.0   0.0   3.0  1.00  0.00          HG  \n"
                   "ATOM      5  CA  ALA A   1       0.0   0.0   4.0  1.00  0.00           C  \n")
    n = pq.strip_all_hydrogens(str(src), str(tmp_path / "b.pdb"))
    kept = (tmp_path / "b.pdb").read_text().splitlines()
    assert n == 2 and len(kept) == 3 and any(" HG " in k for k in kept)      # mercury is not a hydrogen
    assert pq.residue_names(str(tmp_path / "b.pdb")) == {("A", 1): "ALA", ("A", 2): "HG"}


def test_setup_dir_layout(tmp_path, monkeypatch):
    import ldsrc.controller as ctl
    monkeypatch.setattr(ctl, "data_dir", lambda: str(tmp_path / "data"))
    os.makedirs(tmp_path / "data")
    c = ctl.DockingController.__new__(ctl.DockingController)
    c.settings = ctl.Settings(str(tmp_path / "s.json"))
    assert c.setup_dir == ctl.default_setup_dir()
    ok, path = c.set_setup_dir(str(tmp_path / "vpx"))
    assert ok and c.setup_dir == str(tmp_path / "vpx") and c.vina_root == str(tmp_path / "vpx" / "vina") and c.prep_dir == str(tmp_path / "vpx" / "meeko")
    assert not c.set_setup_dir(str(tmp_path / "s.json" / "x"))[0]


@pytest.mark.skipif(not REAL, reason="set VPX_TEST_PYTHON to a Python that has Meeko + RDKit")
def test_real_pdb2pqr_pipeline(tmp_path):
    env = pe.external_env(REAL)
    info = pe.probe(env)["info"]
    if not info.get("pdb2pqr"):
        pytest.skip("pdb2pqr is not installed in VPX_TEST_PYTHON")
    gen = """
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
m = Chem.AddHs(Chem.MolFromSequence("AGSDKHFEYCNQ"), addResidueInfo=True)
ps = AllChem.ETKDGv3(); ps.useRandomCoords = True; ps.randomSeed = 3; AllChem.EmbedMolecule(m, ps)
Chem.MolToPDBFile(m, sys.argv[1])
"""
    g = tmp_path / "gen.py"; g.write_text(gen)
    pdb = tmp_path / "with_h.pdb"
    assert pe.run_stream(env.python_command([str(g), str(pdb)]), env=env.environ()) == 0
    logs = []
    r = pq.run_pdb2pqr(env, str(pdb), str(tmp_path / "out"), "pep", 7.4, log=logs.append)       # input has hydrogens: stripped first
    assert os.path.isfile(r["pdb"]) and os.path.isfile(r["log_path"]) and os.path.isfile(r["csv_path"])
    groups = {t["residue"] for t in r["table"]}
    assert any(g.startswith("HIS") for g in groups) and any(g.startswith("ASP") for g in groups)
    his = [t for t in r["table"] if t["residue"].startswith("HIS")][0]
    assert his["renamed"].startswith("HIS →") or pq.residue_names(r["pdb"])[("A", 6)] in ("HID", "HIE", "HIP")
    low = pq.run_pdb2pqr(env, str(pdb), str(tmp_path / "low"), "pep", 3.0)                     # pH changes the states
    assert pq.residue_names(low["pdb"])[("A", 6)] == "HIP" and any(n.startswith("ASP") for n in [x.split(":")[0] for x in low["renamed"]]) is not None
    assert "ASH" in set(pq.residue_names(low["pdb"]).values()) or "GLH" in set(pq.residue_names(low["pdb"]).values())


@pytest.mark.skipif(not REAL, reason="set VPX_TEST_PYTHON to a Python that has Meeko + RDKit")
def test_renames_include_residues_without_pka(tmp_path):
    """The first His of a chain gets no PROPKA pKa, but its new name (HIP/HID) must still be reported."""
    env = pe.external_env(REAL)
    if not pe.probe(env)["info"].get("pdb2pqr"):
        pytest.skip("pdb2pqr is not installed in VPX_TEST_PYTHON")
    src = os.path.join(os.path.dirname(__file__), "data", "vpx_test_peptide.pdb")
    r = pq.run_pdb2pqr(env, src, str(tmp_path / "o"), "p", 3.0)
    assert "HIS B3: HIS → HIP" in r["renamed"] and r["no_pka"] == ["HIS B3: HIS → HIP"], (r["renamed"], r["no_pka"])
    assert "HIS A3: HIS → HIP" in r["renamed"] and "ASP A2: ASP → ASH" in r["renamed"]
    assert r["flips"] == [] or all(f[0].startswith("ASN") for f in r["flips"])
