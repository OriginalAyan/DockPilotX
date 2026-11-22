import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import conftest

def test_layout_and_project(tmp_path):
    from ldsrc import layout
    from ldsrc.project import Project
    d = str(tmp_path)
    p = layout.make(d, "prepared")
    assert os.path.basename(p).startswith("vinax_") and os.path.isfile(os.path.join(p, layout.NOTE_NAME))
    pr = Project(d)
    r1 = pr.new_revision("prot", "prot.pdb")
    r2 = pr.new_revision("prot", "prot.pdb")
    assert (r1["id"], r2["id"]) == ("r001", "r002") and pr.active() is r2 and os.path.isdir(pr.rev_dir(r2))
    f = pr.stage_file(r2, "clean"); open(f, "w").write("x")
    pr.set_stage(r2, "clean", f); pr.set_stage(r2, "hydrogens", f); pr.set_stage(r2, "pdbqt", f)
    assert pr.latest_stage(r2)[0] == "hydrogens"
    pr.set_stage(r2, "clean", f)                       # redoing clean invalidates all later stages
    assert list(r2["stages"]) == ["clean"]
    lig = os.path.join(d, "a.sdf"); pr.set_inputs([lig], box=([0, 1, 2], [10, 10, 10]), params={"exhaustiveness": 8})
    assert pr.data["ligands"]["inputs"] == ["a.sdf"]
    assert pr.save("0.9.0") and os.path.isfile(os.path.join(d, "project.vinapx"))
    q = Project.open(d)
    assert q.active()["id"] == "r002" and q.ligand_paths() == [lig] and q.data["box"]["size"] == [10, 10, 10]
    assert Project.open(os.path.join(d, "project.vinapx")).workdir == os.path.abspath(d)
    assert json.load(open(q.path))["format"] == "vinapx"

def test_copy_seqres(tmp_path):
    from ldsrc.pdb2pqr_tools import copy_seqres
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    a.write_text("SEQRES   1 A    3  ALA GLY SER\nATOM      1  N   ALA A   1\n")
    b.write_text("ATOM      1  N   ALA A   1\n")
    assert copy_seqres(str(a), str(b)) == 1 and b.read_text().startswith("SEQRES")
    assert copy_seqres(str(a), str(b)) == 0            # not twice
