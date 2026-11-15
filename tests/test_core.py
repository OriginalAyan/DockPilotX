import os, sys, stat, shutil, importlib
sys.path.insert(0, os.path.dirname(__file__))
import conftest, helpers
from ldsrc import vina_manager as vm, pdbqt_tools as pt, docking as dk
FAKE = os.path.join(conftest.HERE, "fake_vina.py")
import pytest

def write(p, t): open(p, "w").write(t); return p

# ---------- vina_manager
def test_asset_names():
    assert vm.asset_name("1.2.7", "win", "x86_64") == "vina_1.2.7_win.exe"
    assert vm.asset_name("v1.2.7", "mac", "aarch64") == "vina_1.2.7_mac_aarch64"
    assert vm.asset_name("1.2.7", "linux", "x86_64") == "vina_1.2.7_linux_x86_64"

def test_parse_version():
    assert vm.parse_version("AutoDock Vina v1.2.5") == "1.2.5"
    assert vm.parse_version("nothing") is None

def test_check_and_resolve(tmp_path):
    exe = tmp_path / "vina_1.2.7_linux_x86_64"; shutil.copy(FAKE, exe); exe.chmod(0o755)
    split = tmp_path / "vina_split_1.2.7_linux_x86_64"; split.write_text("x")
    assert vm.resolve_vina_path(str(tmp_path)) == str(exe)          # folder -> exe, ignores vina_split
    ok, v = vm.check_vina(str(tmp_path)); assert ok and v == "1.2.7"
    ok, msg = vm.check_vina(str(tmp_path / "nope")); assert not ok
    junk = tmp_path / "notvina"; junk.write_text("hello"); junk.chmod(0o755)
    ok, msg = vm.check_vina(str(junk)); assert not ok

def test_download_with_local_server(tmp_path):
    import http.server, threading, functools
    srv_dir = tmp_path / "srv"; (srv_dir / "dl").mkdir(parents=True)
    name = vm.asset_name("9.9.9", "linux", "x86_64")
    data = open(FAKE, "rb").read() + b"#" * 120000   # > 100 kB sanity limit
    (srv_dir / "dl" / name).write_bytes(data)
    h = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(srv_dir))
    h.log_message = lambda *a, **k: None
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), h); port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    seen = []
    path, ver = vm.download_vina(str(tmp_path / "out"), version="9.9.9", system="linux", arch="x86_64",
                                 url_template="http://127.0.0.1:%d/dl/{asset}" % port,
                                 progress=lambda d, t: seen.append(d))
    httpd.shutdown()
    assert os.access(path, os.X_OK) and ver == "1.2.7" and seen and os.path.exists(path + ".sha256")
    assert not os.path.exists(path + ".part")

def test_download_failure_cleans_part(tmp_path):
    with pytest.raises(vm.VinaError):
        vm.download_vina(str(tmp_path), version="1.0.0", system="linux", arch="x86_64",
                         url_template="http://127.0.0.1:9/{asset}", timeout=2)
    assert not any(f.endswith(".part") for r, d, fs in os.walk(tmp_path) for f in fs)

# ---------- pdbqt_tools
def test_receptor_checks(tmp_path):
    good = write(tmp_path / "r.pdbqt", helpers.rec_text())
    assert not [m for m in pt.check_receptor_pdbqt(str(good)) if m[0] == "error"]
    bad = write(tmp_path / "b.pdbqt", helpers.rec_text().replace("OA", "Qq"))
    assert any(m[0] == "error" for m in pt.check_receptor_pdbqt(str(bad)))
    lig = write(tmp_path / "l.pdbqt", helpers.lig_text())
    assert any("ROOT" in m[1] for m in pt.check_receptor_pdbqt(str(lig)))
    assert pt.check_receptor_pdbqt(str(tmp_path / "missing.pdbqt"))[0][0] == "error"

def test_library_split_and_names(tmp_path):
    lib = write(tmp_path / "lib.pdbqt", helpers.lig_text("aspirin") + helpers.lig_text("a/b c") +
                helpers.lig_text("aspirin"))
    single = write(tmp_path / "one.pdbqt", helpers.lig_text("zzz"))
    assert pt.count_ligands(str(lib)) == 3
    items = pt.stage_ligands([str(lib), str(single)], str(tmp_path / "stage"))
    names = sorted(i.name for i in items)
    assert names == ["a_b_c", "aspirin", "aspirin_2", "one"]       # unique, safe names; single keeps file stem
    assert all(not i.problems for i in items)
    assert open(items[0].path).read().count("TORSDOF") == 1

def test_model_wrapped_and_bad(tmp_path):
    wrapped = "MODEL 1\n" + helpers.lig_text("m1") + "ENDMDL\nMODEL 2\n" + helpers.lig_text("m2") + "ENDMDL\n"
    f = write(tmp_path / "w.pdbqt", wrapped)
    items = pt.stage_ligands([str(f)], str(tmp_path / "s"))
    assert [i.name for i in items] == ["m1", "m2"]
    assert "MODEL" not in open(items[0].path).read()
    bad = write(tmp_path / "bad.pdbqt", "ATOM      1  C   LIG A   1       0.0 0.0 0.0\n")
    it = pt.stage_ligands([str(bad)], str(tmp_path / "s2"))
    assert it[0].problems

def test_folder_input_and_missing(tmp_path):
    d = tmp_path / "d"; d.mkdir()
    write(d / "x.pdbqt", helpers.lig_text("x")); write(d / "y.pdbqt", helpers.lig_text("y")); write(d / "z.txt", "no")
    assert [os.path.basename(p) for p in pt.list_pdbqt_files([str(d)])] == ["x.pdbqt", "y.pdbqt"]
    with pytest.raises(FileNotFoundError): pt.list_pdbqt_files([str(d / "nope.pdbqt")])

# ---------- docking
def make_runner(tmp_path, ligs, **kw):
    exe = tmp_path / "vina"; shutil.copy(FAKE, exe); exe.chmod(0o755)
    rec = write(tmp_path / "my receptor.pdbqt", helpers.rec_text())      # space in name on purpose
    params = kw.pop("params", dk.DockingParams())
    return dk.DockingRunner(str(exe), str(rec), (1, 2, 3), (20, 20, 20), params,
                            str(tmp_path / "work dir"), ligs, **kw), exe

def test_config_text(tmp_path):
    p = tmp_path / "c.txt"
    dk.write_vina_config(str(p), (1, 2, 3), (20, 21, 22.5), dk.DockingParams(exhaustiveness=16, seed=7))
    t = p.read_text()
    for s in ("center_x = 1.0000", "size_z = 22.5000", "exhaustiveness = 16", "scoring = vina"):
        assert s in t
    assert "receptor" not in t and "seed" not in t      # seed is per-ligand and goes on the command line

def test_params_validation():
    assert dk.DockingParams().validate() == []
    assert len(dk.DockingParams(exhaustiveness=0, scoring="x", parallel_jobs=0).validate()) == 3

def test_single_and_batch_run_with_failure_and_resume(tmp_path):
    lib = write(tmp_path / "lib.pdbqt", helpers.lig_text("good1") + helpers.lig_text("badone", atom_type="C", extra="REMARK BADLIG\n") +
                helpers.lig_text("good2"))
    events = []
    r, exe = make_runner(tmp_path, [str(lib)], params=dk.DockingParams(parallel_jobs=2, num_modes=2),
                         on_job=lambda j: events.append((j.name, j.status)),
                         on_log=lambda m: None)
    r.prepare(); r.run()
    st = {j.name: j for j in r.jobs}
    assert st["good1"].status == dk.STATUS_DONE and st["good1"].n_poses == 2
    assert st["badone"].status == dk.STATUS_FAILED and "Unknown AutoDock atom type" in st["badone"].error
    assert st["good2"].best is not None
    assert os.path.exists(os.path.join(r.workdir, "summary.csv"))
    import json; m = json.load(open(os.path.join(r.workdir, "manifest.json")))
    assert m["complete"] and m["vina_version"] == "1.2.7" and m["status_counts"]["failed"] == 1
    assert "--config" in m["example_command"] and '"' in m["example_command"]   # space-quoted for display
    assert any(s == dk.STATUS_RUNNING for _, s in events)
    # resume: second run reuses finished jobs, retries the failed one
    r2, _ = make_runner(tmp_path, [str(lib)])
    r2.prepare(); r2.run()
    st2 = {j.name: j for j in r2.jobs}
    assert st2["good1"].status == dk.STATUS_EXISTING and st2["badone"].status == dk.STATUS_FAILED
    # results folder reload
    loaded = dk.load_results_dir(r.workdir)
    assert sorted(j.name for j in loaded) == ["good1", "good2"]

def test_invalid_ligand_skipped(tmp_path):
    bad = write(tmp_path / "bad.pdbqt", "garbage\n")
    r, _ = make_runner(tmp_path, [str(bad)])
    r.prepare(); r.run()
    assert r.jobs[0].status == dk.STATUS_SKIPPED

def test_cancel(tmp_path):
    ligs = []
    for i in range(6): ligs.append(write(tmp_path / ("l%d.pdbqt" % i), helpers.lig_text("l%d" % i)))
    done = []
    r, exe = make_runner(tmp_path, ligs, on_finished=lambda rr: done.append(1))
    r.prepare(); r.start(); r.cancel(); r.join(30)
    assert done and all(j.status in (dk.STATUS_CANCELLED, dk.STATUS_DONE) for j in r.jobs)
    assert any(j.status == dk.STATUS_CANCELLED for j in r.jobs)

def test_preflight(tmp_path):
    exe = tmp_path / "vina"; shutil.copy(FAKE, exe); exe.chmod(0o755)
    rec = write(tmp_path / "r.pdbqt", helpers.rec_text())
    ok = dk.preflight(str(exe), str(rec), (0, 0, 0), (20, 20, 20), dk.DockingParams(), str(tmp_path))
    assert not [m for m in ok if m[0] == "error"]
    big = dk.preflight(str(exe), str(rec), (0, 0, 0), (40, 40, 40), dk.DockingParams(), str(tmp_path))
    assert any("volume" in m[1] for m in big)
    none = dk.preflight("", "", None, None, dk.DockingParams(exhaustiveness=0), "")
    assert len([m for m in none if m[0] == "error"]) >= 5


# ---------- seeds
def test_seed_is_read_from_vina_output(tmp_path):
    import json, csv
    assert dk.parse_seed("Performing docking (random seed: -1805814086) ...") == -1805814086
    assert dk.parse_seed("nothing here") is None
    ligs = [write(tmp_path / ("s%d.pdbqt" % i), helpers.lig_text("s%d" % i)) for i in range(4)]
    r, _ = make_runner(tmp_path, ligs, params=dk.DockingParams(parallel_jobs=2))
    r.prepare(); r.run()
    seeds = [j.seed for j in r.jobs]
    assert all(isinstance(s, int) for s in seeds) and len(set(seeds)) == 4
    for j in r.jobs:                                   # nothing was forced on the command line
        assert "--seed" not in open(j.log_path).readline() and dk.seed_from_log(j.log_path) == j.seed
    rows = list(csv.DictReader(open(os.path.join(r.workdir, "best_per_ligand.csv"))))
    assert sorted(int(x["seed"]) for x in rows) == sorted(seeds)
    m = json.load(open(os.path.join(r.workdir, "manifest.json")))
    assert m["seeds"] == {j.name: j.seed for j in r.jobs} and "Vina chooses" in m["seed_policy"]
    r2, _ = make_runner(tmp_path, ligs); r2.prepare(); r2.run()                     # resume keeps seeds
    assert sorted(j.seed for j in r2.jobs) == sorted(seeds)
    assert sorted(j.seed for j in dk.load_results_dir(r.workdir)) == sorted(seeds)
    fresh = tmp_path / "fixed"; fresh.mkdir()                                       # user-fixed seed
    r3, _ = make_runner(fresh, ligs[:2], params=dk.DockingParams(seed=12345)); r3.prepare(); r3.run()
    assert [j.seed for j in r3.jobs] == [12345, 12345]
    assert "--seed 12345" in open(r3.jobs[0].log_path).readline()
    assert "fixed seed 12345" in json.load(open(os.path.join(r3.workdir, "manifest.json")))["seed_policy"]

# ---------- ligand properties
def make_lig(coords, types=None, torsdof=0, name="p"):
    types = types or ["C"] * len(coords)
    out = ["REMARK  Name = %s" % name, "ROOT"]
    for i, ((x, y, z), t) in enumerate(zip(coords, types), 1):
        out.append("ATOM  %5d %-4s UNL A%4d    %8.3f%8.3f%8.3f  1.00  0.00    %+6.3f %-2s" % (i, "C%d" % i, 1, x, y, z, 0.1, t))
    out += ["ENDROOT", "TORSDOF %d" % torsdof]
    return "\n".join(out) + "\n"

def test_ligand_dimensions_are_orientation_independent(tmp_path):
    import numpy as np
    pts = [(0, 0, 0), (10, 0, 0), (20, 0, 0), (20, 5, 0), (0, 5, 0), (10, 5, 3)]
    a = pt.ligand_properties(make_lig(pts, torsdof=3).splitlines())
    # rotate the same molecule 40 degrees about z and 30 about x: principal-axis dims must not change
    th, ph = np.radians(40), np.radians(30)
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(ph), -np.sin(ph)], [0, np.sin(ph), np.cos(ph)]])
    rot = [tuple(Rx @ Rz @ np.array(p)) for p in pts]
    b = pt.ligand_properties(make_lig(rot, torsdof=3).splitlines())
    assert np.allclose(a["dims"], b["dims"], atol=0.05) and abs(a["longest"] - b["longest"]) < 1e-3
    assert a["dims"][0] >= a["dims"][1] >= a["dims"][2] and a["torsions"] == 3 and a["n_atoms"] == 6
    assert abs(a["longest"] - (400 + 25) ** 0.5) < 0.01        # (0,0,0)-(20,5,0)
    assert abs(a["charge"] - 0.6) < 1e-6

def test_hbond_counts_and_hydrogens():
    p = pt.ligand_properties(make_lig([(0, 0, 0), (1.2, 0, 0), (2, 1, 0), (2.5, 1.5, 0)], ["C", "OA", "N", "HD"]).splitlines())
    assert p["n_atoms"] == 3 and p["n_acceptors"] == 1 and p["n_polar_h"] == 1

def test_scan_records_and_selective_staging(tmp_path):
    lib = write(tmp_path / "lib.pdbqt", make_lig([(0, 0, 0), (5, 0, 0)], name="small") + make_lig([(0, 0, 0), (30, 0, 0)], name="big") + make_lig([(0, 0, 0), (9, 0, 0)], name="mid"))
    recs = pt.scan_ligands([str(lib)])
    assert [r.name for r in recs] == ["small", "big", "mid"] and max(recs, key=lambda r: r.props["longest"]).name == "big"
    assert pt.scan_ligands([str(lib)], skip_keys=[r.key for r in recs]) == []      # adding same file twice is a no-op
    items = pt.stage_records([recs[1], recs[2]], str(tmp_path / "st"))             # user removed 'small' from the table
    assert [i.name for i in items] == ["big", "mid"] and not os.path.exists(tmp_path / "st" / "small.pdbqt")
    # records go straight into the runner
    r, _ = make_runner(tmp_path, [recs[0], recs[2]]); r.prepare(); r.run()
    assert sorted(j.name for j in r.jobs) == ["mid", "small"] and all(j.status == dk.STATUS_DONE for j in r.jobs)


def test_ligand_csv_roundtrip(tmp_path):
    import csv
    lib = write(tmp_path / "lib.pdbqt", make_lig([(0, 0, 0), (5, 0, 0)], name="a") + make_lig([(0, 0, 0), (30, 0, 0)], name="b", torsdof=2))
    bad = write(tmp_path / "bad.pdbqt", "garbage\n")
    recs = pt.scan_ligands([str(lib), str(bad)])
    out = tmp_path / "t.csv"; pt.write_ligand_csv(recs, str(out))
    rows = list(csv.DictReader(open(out)))
    assert [r["ligand"] for r in rows] == ["bad", "a", "b"] and rows[0]["valid"] == "no" and rows[0]["problems"]   # files are read in sorted order
    assert rows[2]["torsions"] == "2" and rows[2]["longest_A"] == "30.00" and rows[2]["index_in_file"] == "2"


# ---------- poses, CSV, extraction, analysis store, run log
def _run_two(tmp_path):
    ligs = [write(tmp_path / ("p%d.pdbqt" % i), helpers.lig_text("p%d" % i)) for i in range(2)]
    r, _ = make_runner(tmp_path, ligs, params=dk.DockingParams(num_modes=3)); r.prepare(); r.run()
    return r

def test_all_poses_with_rmsd_in_summary_csv(tmp_path):
    import csv
    r = _run_two(tmp_path)
    j = r.jobs[0]
    assert len(j.poses) == 3 and j.poses[1][1:] == (1.1, 2.2) and j.best == min(p[0] for p in j.poses)
    rows = list(csv.DictReader(open(os.path.join(r.workdir, "summary.csv"))))
    assert len(rows) == 6 and [x["pose"] for x in rows[:3]] == ["1", "2", "3"]
    assert rows[1]["rmsd_lb"] == "1.100" and rows[1]["rmsd_ub"] == "2.200" and "time" not in " ".join(rows[0].keys()).lower()
    assert len(list(csv.DictReader(open(os.path.join(r.workdir, "best_per_ligand.csv"))))) == 2

def test_failed_ligand_row_and_filter(tmp_path):
    import csv
    bad = write(tmp_path / "b.pdbqt", helpers.lig_text("b", extra="REMARK BADLIG\n"))
    r, _ = make_runner(tmp_path, [str(bad)]); r.prepare(); r.run()
    rows = list(csv.DictReader(open(os.path.join(r.workdir, "summary.csv"))))
    assert len(rows) == 1 and rows[0]["status"] == "failed" and rows[0]["pose"] == "" and rows[0]["error"]

def test_extract_pose(tmp_path):
    r = _run_two(tmp_path)
    j = r.jobs[0]
    dst = dk.extract_pose(j.out_path, 2, str(tmp_path / "cache" / "x.pdbqt"))
    txt = open(dst).read()
    assert "MODEL" not in txt and txt.count("REMARK VINA RESULT") == 1 and "1.100" in txt and "   1.500" in txt
    with pytest.raises(ValueError): dk.extract_pose(j.out_path, 9, str(tmp_path / "cache" / "y.pdbqt"))

def test_analysis_store_and_csv(tmp_path):
    import csv
    from ldsrc.analysis_store import AnalysisStore
    r = _run_two(tmp_path); a = AnalysisStore(r.workdir)
    n = r.jobs[0].name
    a.set(n, 2, rating=5, hbonds=3, clashes=0, hbond_partners=["LIG O1 - ASP A45 OD1"])
    a.set(n, 3, rating=2)
    b = AnalysisStore(r.workdir)                                        # reloaded from disk
    assert b.get(n, 2)["hbonds"] == 3 and [x[:2] for x in b.rated(4)] == [(n, 2)] and len(b.rated(1)) == 2
    b.set(n, 3, rating=0); assert len(b.rated(1)) == 1                  # clearing a rating
    out = tmp_path / "rated.csv"
    dk.write_pose_csv(r.jobs, str(out), b, only=lambda name, pose, rec: int(rec.get("rating") or 0) >= 4)
    rows = list(csv.DictReader(open(out)))
    assert len(rows) == 1 and rows[0]["pose"] == "2" and rows[0]["hbonds"] == "3" and rows[0]["clashes"] == "0" and "ASP" in rows[0]["hbond_partners"]

def test_run_log_file(tmp_path):
    r = _run_two(tmp_path)
    text = open(os.path.join(r.workdir, "run_log.txt")).read()
    assert "Run started: Vina 1.2.7" in text and "best" in text and text.count("Docking p") >= 2


def test_interaction_cache_roundtrip(tmp_path):
    from ldsrc.analysis_store import InteractionCache
    c = InteractionCache(str(tmp_path))
    assert c.load("lig", 2) == {}
    assert c.save("lig", 2, hbonds=[[0, 5, 1, 3.1]], n_ligand_atoms=4, n_receptor_atoms=9)
    assert c.save("lig", 2, clashes=[[1, 2, 2.4]])
    d = InteractionCache(str(tmp_path)).load("lig", 2)
    assert d["hbonds"] == [[0, 5, 1, 3.1]] and d["clashes"] == [[1, 2, 2.4]] and d["n_receptor_atoms"] == 9
