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


def test_approval_and_flips(tmp_path):
    from ldsrc.project import Project
    from ldsrc import pdb2pqr_tools as pq
    pr = Project(str(tmp_path)); r = pr.new_revision("p", "p.pdb")
    q = pr.rev_dir(r) + "/p.pdbqt"; open(q, "w").write("x")
    pr.set_stage(r, "pdbqt", q)
    assert pr.approved() == [] and pr.approve(r) and pr.approved() == [(r, q)]
    pr.set_stage(r, "pdbqt", q); assert not r["approved"] and pr.approved() == []
    os.remove(q); r["approved"] = True; assert pr.approved() == []          # a missing file is never offered

    def atom(n, res, name, num, xyz, chain="A"):
        return "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n" % (n, name, res, chain, num, *xyz)
    a, b = tmp_path / "a.pdb", tmp_path / "b.pdb"
    a.write_text(atom(1, "ASN", "OD1", 6, (9, 0, 0)) + atom(2, "ASN", "ND2", 6, (7, 2, 0)) + atom(3, "GLN", "OE1", 7, (1, 1, 1)) +
                 atom(4, "GLN", "NE2", 7, (2, 2, 2)) + atom(5, "HIS", "ND1", 8, (0, 0, 0)) + atom(6, "GLU", "OE1", 9, (0, 0, 0)))
    b.write_text(atom(1, "ASN", "OD1", 6, (7, 2, 0)) + atom(2, "ASN", "ND2", 6, (9, 0, 0)) + atom(3, "GLN", "OE1", 7, (1, 1, 1.1)) +
                 atom(4, "GLN", "NE2", 7, (2, 2, 2.1)) + atom(5, "HID", "ND1", 8, (0, 0, 0.2)) + atom(6, "GLH", "OE2", 9, (3, 3, 3)))
    fl = pq.find_flips(str(a), str(b))
    assert [f[0] for f in fl] == ["ASN A6"] and 2.8 < fl[0][1] < 2.9, fl      # Gln/His shifts under 1 A, Glu/GLH name swap ignored
    t = pq.build_table([{"group": "N+", "number": 6, "chain": "A", "pka": 8.0, "model_pka": 8.0}], {("A", 6): "HIS"}, {("A", 6): "HID"}, 7.4)
    assert t[0]["renamed"] == ""          # a chain-end row never reports the residue's own rename
