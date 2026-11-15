import conftest  # noqa: F401  (sets up the ldsrc package)
import importlib
fp = importlib.import_module("ldsrc.format_preference")


class FakeSettings:
    def __init__(self, mapping):
        self.suffix_to_format_name = dict(mapping); self.saved = 0
    def save(self): self.saved += 1

class FakeFormat:
    def __init__(self, name, suffixes, default_for=(), nicknames=()):
        self.name, self.suffixes, self.default_for, self.nicknames = name, list(suffixes), list(default_for), list(nicknames)

class FakeFormats:
    def __init__(self, mapping, include_pdbt=True):
        self.settings = FakeSettings(mapping)
        self.native = FakeFormat("AutoDock PDBQT", [".pdbqt"], (), ["pdbqt"])
        self.pdbt = FakeFormat("PDBT", [".pdbt", ".pdbqt"], [".pdbt", ".pdbqt"], ["pdbt"])
        self.formats = [self.native] + ([self.pdbt] if include_pdbt else [])
        self._suffix_to_formats = {".pdbqt": list(self.formats)}
        self._user_response_cache = {".pdbqt": self.pdbt}

def setup_function(_):
    fp._SAVED = None

def test_native_becomes_default_then_restore_roundtrip():
    f = FakeFormats({".pdbqt": "PDBT"})
    changed, previous = fp.prefer_native_pdbqt(f)
    assert changed and previous == "PDBT" and fp.override_active()
    assert ".pdbqt" in f.native.default_for and ".pdbqt" not in f.pdbt.default_for
    assert f.settings.suffix_to_format_name[".pdbqt"] == "AutoDock PDBQT"
    # applying twice must not lose the true original
    fp.prefer_native_pdbqt(f)
    assert fp.restore_pdbqt_preference(f) is True
    assert f.native.default_for == [] and f.pdbt.default_for == [".pdbt", ".pdbqt"]
    assert f.settings.suffix_to_format_name[".pdbqt"] == "PDBT"
    assert f._user_response_cache[".pdbqt"] is f.pdbt
    assert not fp.override_active() and fp.restore_pdbqt_preference(f) is False

def test_restore_removes_key_that_did_not_exist():
    f = FakeFormats({}, include_pdbt=False)
    fp.prefer_native_pdbqt(f)
    assert f.settings.suffix_to_format_name[".pdbqt"] == "AutoDock PDBQT"
    fp.restore_pdbqt_preference(f)
    assert ".pdbqt" not in f.settings.suffix_to_format_name and f.native.default_for == []

def test_missing_native_reported():
    f = FakeFormats({".pdbqt": "PDBT"}); f._suffix_to_formats[".pdbqt"] = [f.pdbt]
    try:
        fp.prefer_native_pdbqt(f)
    except LookupError:
        pass
    else:
        raise AssertionError
    assert not fp.override_active()
