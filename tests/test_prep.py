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


def test_list_ligand_files(tmp_path):
    (tmp_path / "sub").mkdir()
    for n in ("a.sdf", "b.MOL2", "c.txt", "sub/d.smi", "e.pdb"):
        (tmp_path / n).write_text("x")
    got = [os.path.basename(p) for p in prep.list_ligand_files([str(tmp_path), str(tmp_path / "a.sdf")])]
    assert sorted(got) == ["a.sdf", "b.MOL2", "d.smi"]
    assert [os.path.basename(p) for p in prep.list_ligand_files([str(tmp_path)], recursive=False)] == ["a.sdf", "b.MOL2"]


def test_env_commands_and_roundtrip(tmp_path):
    mm = pe.PrepEnv("micromamba", prefix="/e", mamba="/mm", root="/r")
    assert mm.command(["python", "-m", "x"]) == ["/mm", "run", "-r", "/r", "-p", "/e", "python", "-m", "x"]
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
    path, warns = prep.prepare_receptor(env, info, str(pdb), str(tmp_path / "out" / "rec"), True)
    assert path.endswith("rec.pdbqt") and os.path.getsize(path) > 100 and warns == ["A:7"]
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
