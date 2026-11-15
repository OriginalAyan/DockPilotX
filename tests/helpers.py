import os
def lig_text(name="LIG1", atom_type="C", extra=""):
    return f"""REMARK  Name = {name}
REMARK SMILES CCO
ROOT
ATOM      1  C1  UNL     1       0.000   0.000   0.000  1.00  0.00    +0.000 {atom_type:<2}
ATOM      2  O1  UNL     1       1.200   0.000   0.000  1.00  0.00    -0.300 OA
ENDROOT
{extra}TORSDOF 0
"""
def rec_text():
    return """ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00    -0.350 N 
ATOM      2  HN  ALA A   1       0.500   0.500   0.000  1.00  0.00    +0.160 HD
ATOM      3  CA  ALA A   1       1.000   1.000   0.000  1.00  0.00    +0.180 C 
ATOM      4  O   ALA A   1       2.000   1.000   0.000  1.00  0.00    -0.270 OA
TER
"""
