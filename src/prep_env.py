"""Set up and use a private Python environment that contains Meeko (and RDKit).

Pure Python (no ChimeraX imports). Three kinds of environment are supported:

* ``micromamba``: a micromamba program is downloaded, then a conda-forge environment is created.
* ``venv``: a standard virtual environment plus ``pip install`` (fallback when micromamba is unavailable).
* ``external``: an existing Python (or conda/venv folder) that already has Meeko; nothing is installed.

Every Meeko command runs in a separate process, so nothing is installed into ChimeraX.
"""

import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

from .vina_manager import no_window_flags

MICROMAMBA_URLS = [
    "https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-{plat}{ext}",
]
MICROMAMBA_API = "https://micro.mamba.pm/api/micromamba/{plat}/latest"
CONDA_PACKAGES = ["python=3.11", "meeko", "rdkit", "scipy", "numpy", "gemmi", "prody"]
PIP_CORE = ["meeko", "rdkit", "scipy", "numpy", "gemmi"]
PIP_OPTIONAL = ["prody"]


class PrepError(RuntimeError):
    pass


# ---------------------------------------------------------------------------------------------
def mamba_platform():
    machine = platform.machine().lower()
    arm = machine in ("arm64", "aarch64")
    if sys.platform.startswith("win"):
        return "win-64", ".exe"
    if sys.platform == "darwin":
        return ("osx-arm64" if arm else "osx-64"), ""
    return ("linux-aarch64" if arm else "linux-64"), ""


def python_in(prefix):
    if sys.platform.startswith("win"):
        for rel in ("python.exe", os.path.join("Scripts", "python.exe")):
            p = os.path.join(prefix, rel)
            if os.path.isfile(p):
                return p
        return os.path.join(prefix, "python.exe")
    for name in ("python", "python3"):
        p = os.path.join(prefix, "bin", name)
        if os.path.isfile(p):
            return p
    return os.path.join(prefix, "bin", "python")


class PrepEnv:
    """Description of an environment; ``command()`` builds an argv that runs a program inside it."""

    def __init__(self, kind, prefix=None, mamba=None, root=None, python=None):
        self.kind = kind            # 'micromamba' | 'venv' | 'external'
        self.prefix = prefix        # env folder (micromamba / venv / external conda or venv folder)
        self.mamba = mamba          # micromamba executable (micromamba kind)
        self.root = root            # micromamba root prefix
        self.python = python or (python_in(prefix) if kind == "micromamba" and prefix else None)

    def to_dict(self):
        return {"kind": self.kind, "prefix": self.prefix, "mamba": self.mamba, "root": self.root,
                "python": self.python}

    @classmethod
    def from_dict(cls, d):
        if not isinstance(d, dict) or not d.get("kind"):
            return None
        return cls(d["kind"], d.get("prefix"), d.get("mamba"), d.get("root"), d.get("python"))

    def command(self, argv):
        """argv (list) -> full command that runs it inside the environment.

        Python is started directly (with ``environ()``): ``micromamba run`` goes through a batch file on
        Windows, which breaks on multi-line or special-character arguments."""
        if argv and argv[0] == "python":
            return [self.python or python_in(self.prefix)] + list(argv[1:])
        if self.kind == "micromamba":
            return [self.mamba, "run", "-r", self.root, "-p", self.prefix] + list(argv)
        return list(argv)

    def environ(self):
        """Environment variables for running programs of this environment (None = inherit)."""
        if self.kind == "external" or not self.prefix:
            return None
        env = dict(os.environ)
        p = self.prefix
        if sys.platform.startswith("win"):
            extra = [p, os.path.join(p, "Library", "mingw-w64", "bin"), os.path.join(p, "Library", "usr", "bin"),
                     os.path.join(p, "Library", "bin"), os.path.join(p, "Scripts"), os.path.join(p, "bin")]
        else:
            extra = [os.path.join(p, "bin")]
        env["PATH"] = os.pathsep.join(extra + [env.get("PATH", "")])
        env["PYTHONNOUSERSITE"] = "1"
        if self.kind == "micromamba":
            env["CONDA_PREFIX"] = p
        env.pop("PYTHONHOME", None)
        env.pop("PYTHONPATH", None)
        return env

    def python_command(self, args):
        return self.command(["python"] + list(args))

    def exists(self):
        if self.kind == "micromamba":
            return bool(self.mamba and os.path.isfile(self.mamba) and self.prefix and os.path.isdir(self.prefix))
        return bool(self.python and os.path.isfile(self.python))


# ---------------------------------------------------------------------------------------------
def run_stream(cmd, log=None, cancelled=None, cwd=None, env=None, timeout=None, idle_timeout=None):
    """Run *cmd*, sending each output line to *log*. Returns the exit code (-1 if cancelled).

    *timeout*: kill after this many seconds in total; *idle_timeout*: kill after this many seconds
    without any output. Both raise PrepError.
    """
    import threading
    import time
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=cwd, env=env,
                                creationflags=no_window_flags())
    except OSError as err:
        raise PrepError("Could not run %s: %s" % (cmd[0], err))
    state = {"last": time.time(), "why": None}
    stop = threading.Event()

    def watchdog():
        start = time.time()
        while not stop.wait(2.0):
            now = time.time()
            if timeout and now - start > timeout:
                state["why"] = "took longer than %d s" % timeout
            elif idle_timeout and now - state["last"] > idle_timeout:
                state["why"] = "produced no output for %d s" % idle_timeout
            if state["why"]:
                proc.kill()
                return

    if timeout or idle_timeout:
        threading.Thread(target=watchdog, daemon=True).start()
    try:
        for raw in iter(proc.stdout.readline, b""):
            state["last"] = time.time()
            line = raw.decode("utf-8", "replace").rstrip()
            if line and log:
                log(line)
            if cancelled is not None and cancelled():
                proc.kill()
                proc.wait()
                return -1
        code = proc.wait()
        if state["why"]:
            raise PrepError("%s %s and was stopped." % (os.path.basename(cmd[0]), state["why"]))
        return code
    finally:
        stop.set()


_PROBE = r'''
import json, sys, importlib.util as u
d = {"python": sys.version.split()[0]}
try:
    import importlib.metadata as m
    d["meeko"] = m.version("meeko")
except Exception as e:
    d["meeko"] = None; d["meeko_error"] = str(e)
for mod in ("rdkit", "scipy", "numpy", "gemmi", "prody"):
    d[mod] = u.find_spec(mod) is not None
for tool in ("mk_prepare_ligand", "mk_prepare_receptor"):
    try:
        d["cli_" + tool] = u.find_spec("meeko.cli." + tool) is not None
    except Exception:
        d["cli_" + tool] = False
print("VPXPROBE" + json.dumps(d))
'''


def probe(env, timeout=180):
    """Return {'ok': bool, 'message': str, 'info': dict}. ok means Meeko and RDKit can be imported."""
    if env is None or not env.exists():
        return {"ok": False, "message": "No preparation environment is set up.", "info": {}}
    lines = []
    script = None
    try:
        fd, script = tempfile.mkstemp(suffix=".py", prefix="vpx_probe_")
        with os.fdopen(fd, "w") as f:
            f.write(_PROBE)
        code = run_stream(env.python_command([script]), log=lines.append, timeout=timeout, env=env.environ())
    except (PrepError, OSError) as err:
        return {"ok": False, "message": str(err), "info": {}}
    finally:
        if script:
            try:
                os.remove(script)
            except OSError:
                pass
    info = {}
    for line in lines:
        if line.startswith("VPXPROBE"):
            try:
                info = json.loads(line[len("VPXPROBE"):])
            except ValueError:
                pass
    if not info:
        return {"ok": False, "message": "The environment did not respond (exit code %s): %s" % (
            code, " | ".join(lines[-3:])[:300]), "info": {}}
    if not info.get("meeko"):
        return {"ok": False, "message": "Meeko is not installed in this environment.", "info": info}
    if not info.get("rdkit"):
        return {"ok": False, "message": "RDKit is missing from this environment.", "info": info}
    return {"ok": True, "message": "Meeko %s ready (Python %s%s)." % (
        info["meeko"], info["python"], "" if info.get("prody") else ", no ProDy"), "info": info}


def tool_argv(env, info, tool):
    """argv prefix for 'mk_prepare_ligand' / 'mk_prepare_receptor' in *env*."""
    if info.get("cli_" + tool):
        return env.command(["python", "-m", "meeko.cli." + tool])
    return env.command([tool + ".py"])


def tool_help(env, info, tool, timeout=120):
    lines = []
    try:
        run_stream(tool_argv(env, info, tool) + ["--help"], log=lines.append, timeout=timeout, env=env.environ())
    except PrepError:
        return ""
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------
def download_micromamba(dest_dir, log=None, cancelled=None, timeout=60):
    """Download the micromamba executable into *dest_dir*; returns its path."""
    plat, ext = mamba_platform()
    os.makedirs(dest_dir, exist_ok=True)
    final = os.path.join(dest_dir, "micromamba" + (".exe" if ext else ""))
    errors = []
    for tmpl in MICROMAMBA_URLS:
        url = tmpl.format(plat=plat, ext=ext)
        try:
            if log:
                log("Downloading micromamba: %s" % url)
            _fetch(url, final + ".part", cancelled, timeout)
            os.replace(final + ".part", final)
            _make_executable(final)
            return final
        except PrepError:
            raise
        except Exception as err:
            errors.append("%s: %s" % (url, err))
    url = MICROMAMBA_API.format(plat=plat)
    try:
        if log:
            log("Downloading micromamba: %s" % url)
        arch = os.path.join(dest_dir, "micromamba.tar.bz2")
        _fetch(url, arch, cancelled, timeout)
        with tarfile.open(arch, "r:bz2") as tf:
            member = next((m for m in tf.getmembers()
                           if os.path.basename(m.name) in ("micromamba", "micromamba.exe") and m.isfile()), None)
            if member is None:
                raise PrepError("micromamba was not found in the downloaded archive.")
            with tf.extractfile(member) as src, open(final, "wb") as dst:
                shutil.copyfileobj(src, dst)
        os.remove(arch)
        _make_executable(final)
        return final
    except PrepError:
        raise
    except Exception as err:
        errors.append("%s: %s" % (url, err))
    raise PrepError("Could not download micromamba.\n" + "\n".join(errors))


def _fetch(url, path, cancelled, timeout):
    req = urllib.request.Request(url, headers={"User-Agent": "ChimeraX-VinaPilotX"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(path, "wb") as f:
        while True:
            if cancelled is not None and cancelled():
                f.close()
                os.remove(path)
                raise PrepError("Cancelled.")
            block = resp.read(1 << 16)
            if not block:
                break
            f.write(block)
    if os.path.getsize(path) < 100_000:
        os.remove(path)
        raise PrepError("The downloaded file is too small to be micromamba.")


def _make_executable(path):
    try:
        os.chmod(path, os.stat(path).st_mode | 0o111)
    except OSError:
        pass


def _quiet(log):
    """micromamba redraws progress bars; keep one line per distinct message."""
    if log is None:
        return None
    last = {"key": None}

    def out(line):
        key = line.split("  ")[0].strip()[:60]
        if "Starting" in line and key == last["key"]:
            return
        last["key"] = key
        log(line)
    return out


def short_base_dir(default_dir):
    """Folder for the conda environment. Windows limits paths to 260 characters and conda packages nest
    deeply, so use a short folder (drive root, else %LOCALAPPDATA%) there; elsewhere *default_dir*."""
    if not sys.platform.startswith("win"):
        return default_dir
    local = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
    for root in (os.path.splitdrive(local)[0] + "\\vpx", os.path.join(local, "vpx")):
        try:
            os.makedirs(root, exist_ok=True)
            test = os.path.join(root, ".write_test")
            with open(test, "w") as f:
                f.write("x")
            os.remove(test)
            return os.path.join(root, "mm")
        except OSError:
            continue
    return default_dir


def create_micromamba_env(base_dir, log=None, cancelled=None, packages=None):
    """Download micromamba and create the conda-forge environment under *base_dir*."""
    mamba = os.path.join(base_dir, "micromamba" + (".exe" if sys.platform.startswith("win") else ""))
    if not os.path.isfile(mamba):
        mamba = download_micromamba(base_dir, log, cancelled)
    root = os.path.join(base_dir, "mamba_root")
    prefix = os.path.join(base_dir, "env")
    if os.path.isdir(prefix):
        shutil.rmtree(prefix, ignore_errors=True)
    cmd = [mamba, "create", "-y", "--no-rc", "-r", root, "-p", prefix, "-c", "conda-forge",
           "--override-channels"] + list(packages or CONDA_PACKAGES)
    if log:
        log("Creating the environment (this can take several minutes) ...")
    try:
        code = run_stream(cmd, _quiet(log), cancelled, idle_timeout=900)
    except PrepError:
        shutil.rmtree(base_dir, ignore_errors=True)
        raise
    if code == -1:
        shutil.rmtree(base_dir, ignore_errors=True)
        raise PrepError("Cancelled.")
    if code != 0:
        if log:
            log("Removing the partly created environment ...")
        shutil.rmtree(base_dir, ignore_errors=True)
        raise PrepError("micromamba could not create the environment (exit code %d)." % code)
    return PrepEnv("micromamba", prefix=prefix, mamba=mamba, root=root)


def _python_works(exe):
    """True if *exe* is a Python (>=3.9) that can create virtual environments with pip."""
    try:
        out = subprocess.run([exe, "-c", "import sys, venv, ensurepip; print(sys.version_info[:2] >= (3, 9))"],
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=30,
                             creationflags=no_window_flags())
        return out.returncode == 0 and b"True" in out.stdout
    except (OSError, subprocess.SubprocessError):
        return False


def find_python_for_venv():
    """A Python able to create a venv: system Python first, then the one ChimeraX ships."""
    cands = []
    for name in ("python3", "python", "py"):
        p = shutil.which(name)
        if p and "WindowsApps" not in p:
            cands.append(p)
    ver = "%d.%d" % sys.version_info[:2]
    for rel in ("python.exe", "bin/python3", "bin/python%s" % ver, "bin/python"):
        cands.append(os.path.join(sys.exec_prefix, rel))
    cands.append(os.path.join(os.path.dirname(sys.executable), "python.exe"))
    seen = set()
    for c in cands:
        if c in seen or not os.path.isfile(c):
            continue
        seen.add(c)
        if _python_works(c):
            return c
    return None


def create_venv_env(base_dir, log=None, cancelled=None):
    """Fallback: virtual environment + pip install."""
    py = find_python_for_venv()
    if py is None:
        raise PrepError("No usable Python (3.9 or newer, with venv) was found for the fallback environment.")
    prefix = os.path.join(base_dir, "venv")
    if os.path.isdir(prefix):
        shutil.rmtree(prefix, ignore_errors=True)
    os.makedirs(base_dir, exist_ok=True)
    if log:
        log("Creating a virtual environment with %s" % py)
    code = run_stream([py, "-m", "venv", prefix], log, cancelled)
    if code != 0:
        raise PrepError("Could not create the virtual environment (exit code %d)." % code)
    vpy = python_in(prefix)
    for args, optional in (([ "install", "--upgrade", "pip"], True), (["install"] + PIP_CORE, False),
                           (["install"] + PIP_OPTIONAL, True)):
        if log:
            log("pip " + " ".join(args))
        code = run_stream([vpy, "-m", "pip"] + args + ["--disable-pip-version-check"], log, cancelled,
                          idle_timeout=900)
        if code == -1:
            raise PrepError("Cancelled.")
        if code != 0 and not optional:
            raise PrepError("pip could not install %s (exit code %d)." % (", ".join(PIP_CORE), code))
        if code != 0 and log:
            log("(optional step failed; continuing)")
    return PrepEnv("venv", prefix=prefix, python=vpy)


def external_env(path):
    """Environment from a python executable, or from a conda/venv folder containing one."""
    path = os.path.expanduser(str(path).strip().strip('"'))
    if os.path.isfile(path):
        return PrepEnv("external", prefix=None, python=os.path.abspath(path))
    if os.path.isdir(path):
        py = python_in(path)
        if os.path.isfile(py):
            return PrepEnv("external", prefix=os.path.abspath(path), python=py)
    return None
