import os, sys, types
SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src")
pkg = types.ModuleType("ldsrc"); pkg.__path__ = [SRC]; sys.modules["ldsrc"] = pkg
HERE = os.path.dirname(os.path.abspath(__file__))
