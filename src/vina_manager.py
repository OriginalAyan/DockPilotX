"""Locate, verify and download the AutoDock Vina command-line program.

Pure Python: nothing here imports ChimeraX, so it can be unit-tested outside
ChimeraX. Network access uses only the standard library (urllib), which in
ChimeraX uses ChimeraX's bundled certificates.
"""

import hashlib
import json
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import urllib.error
import urllib.request

GITHUB_REPO = "ccsb-scripps/AutoDock-Vina"
API_LATEST = "https://api.github.com/repos/%s/releases/latest" % GITHUB_REPO
WEB_LATEST = "https://github.com/%s/releases/latest" % GITHUB_REPO
DOWNLOAD_URL = "https://github.com/%s/releases/download/{tag}/{asset}" % GITHUB_REPO

# Used only if the GitHub API and the redirect lookup are both unreachable.
FALLBACK_VERSION = "1.2.7"

_VERSION_RE = re.compile(r"v?(\d+\.\d+\.\d+)")


class VinaError(RuntimeError):
    """Raised for problems locating, verifying or downloading Vina."""


def no_window_flags():
    """subprocess creationflags that stop a console window flashing on Windows."""
    if sys.platform.startswith("win"):
        return getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)
    return 0


# --------------------------------------------------------------------------
# platform / asset naming
# --------------------------------------------------------------------------
def detect_platform():
    """Return (system, arch) as used in Vina release asset names.

    system: 'win', 'mac' or 'linux'; arch: 'x86_64' or 'aarch64'.
    """
    if sys.platform.startswith("win"):
        system = "win"
    elif sys.platform == "darwin":
        system = "mac"
    else:
        system = "linux"
    machine = platform.machine().lower()
    if machine in ("arm64", "aarch64"):
        arch = "aarch64"
    else:
        arch = "x86_64"
    return system, arch


def asset_name(version, system=None, arch=None):
    """Release asset file name, e.g. vina_1.2.7_win.exe or vina_1.2.7_mac_aarch64.

    The Windows build is x86-64 only (it also runs on Windows-on-ARM through
    emulation).
    """
    if system is None or arch is None:
        d_system, d_arch = detect_platform()
        system = system or d_system
        arch = arch or d_arch
    version = version.lstrip("v")
    if system == "win":
        return "vina_%s_win.exe" % version
    if system in ("mac", "linux"):
        return "vina_%s_%s_%s" % (version, system, arch)
    raise VinaError("Unsupported platform: %s/%s" % (system, arch))


# --------------------------------------------------------------------------
# verifying an executable
# --------------------------------------------------------------------------
def parse_version(text):
    """Extract 'x.y.z' from Vina's --version output, or None."""
    m = _VERSION_RE.search(text or "")
    return m.group(1) if m else None


def version_tuple(version):
    try:
        return tuple(int(p) for p in version.split("."))
    except (AttributeError, ValueError):
        return (0,)


def resolve_vina_path(path):
    """Accept either the executable or the folder that contains it.

    Returns an absolute path to a candidate file, or None.
    """
    if not path:
        return None
    path = os.path.expanduser(str(path).strip().strip('"'))
    if os.path.isfile(path):
        return os.path.abspath(path)
    if os.path.isdir(path):
        candidates = []
        for name in sorted(os.listdir(path)):
            low = name.lower()
            if not low.startswith("vina") or "split" in low:
                continue
            full = os.path.join(path, name)
            if os.path.isfile(full) and not low.endswith((".txt", ".md", ".pdf", ".zip", ".sha256")):
                candidates.append(full)
        if candidates:
            # prefer the newest-looking version, then the plain name
            candidates.sort(key=lambda p: (version_tuple(parse_version(os.path.basename(p)) or "0"),
                                           os.path.basename(p)), reverse=True)
            return os.path.abspath(candidates[0])
    return None


def check_vina(path, timeout=20):
    """Run '<path> --version'. Return (ok, version_or_message)."""
    exe = resolve_vina_path(path)
    if exe is None:
        return False, "No Vina executable found at: %s" % path
    try:
        proc = subprocess.run(
            [exe, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=timeout, creationflags=no_window_flags())
    except PermissionError:
        return False, "Not permitted to run %s (check that the file is executable)." % exe
    except OSError as err:
        return False, "Could not run %s: %s" % (exe, err)
    except subprocess.TimeoutExpired:
        return False, "%s did not respond to --version within %d s." % (exe, timeout)
    text = proc.stdout.decode("utf-8", "replace")
    version = parse_version(text)
    if "vina" not in text.lower() and version is None:
        return False, "%s does not look like AutoDock Vina (output: %s)" % (exe, text.strip()[:200])
    if version is None:
        return True, "unknown"
    return True, version


def find_vina_on_system(extra_dirs=()):
    """Look for an already-installed Vina. Returns a path or None."""
    names = ["vina", "vina.exe", "autodock_vina"]
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    for d in extra_dirs:
        found = resolve_vina_path(d)
        if found:
            return found
    return None


# --------------------------------------------------------------------------
# latest release lookup + download
# --------------------------------------------------------------------------
def _request(url, timeout):
    req = urllib.request.Request(url, headers={"User-Agent": "ChimeraX-VinaPilotX",
                                               "Accept": "application/vnd.github+json"})
    return urllib.request.urlopen(req, timeout=timeout)


def latest_release(timeout=15):
    """Return {'version': 'x.y.z', 'tag': 'vx.y.z', 'digests': {asset: sha256}}.

    Tries the GitHub API first (gives checksums), then the /releases/latest
    redirect, and finally a built-in fallback version.
    """
    try:
        with _request(API_LATEST, timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        tag = data["tag_name"]
        digests = {}
        for a in data.get("assets", []):
            d = a.get("digest") or ""
            if d.startswith("sha256:"):
                digests[a["name"]] = d.split(":", 1)[1]
        return {"version": tag.lstrip("v"), "tag": tag, "digests": digests, "source": "api"}
    except Exception:
        pass
    try:
        class _NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *a, **k):
                return None
        opener = urllib.request.build_opener(_NoRedirect)
        req = urllib.request.Request(WEB_LATEST, headers={"User-Agent": "ChimeraX-VinaPilotX"})
        try:
            opener.open(req, timeout=timeout)
        except urllib.error.HTTPError as err:
            loc = err.headers.get("Location", "")
            m = re.search(r"/tag/(v?\d+\.\d+\.\d+)", loc)
            if m:
                tag = m.group(1)
                return {"version": tag.lstrip("v"), "tag": tag, "digests": {}, "source": "redirect"}
    except Exception:
        pass
    return {"version": FALLBACK_VERSION, "tag": "v" + FALLBACK_VERSION, "digests": {},
            "source": "fallback"}


def sha256_of(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            block = f.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def download_vina(dest_root, version=None, progress=None, cancelled=None, timeout=30,
                  system=None, arch=None, url_template=DOWNLOAD_URL):
    """Download the Vina executable for this platform.

    dest_root: folder under which vina-<version>/ is created.
    progress(done_bytes, total_bytes_or_None) and cancelled() -> bool are optional.
    Returns (path_to_executable, version). Raises VinaError.
    """
    info = latest_release(timeout) if version is None else {
        "version": version.lstrip("v"), "tag": "v" + version.lstrip("v"), "digests": {}, "source": "given"}
    ver = info["version"]
    name = asset_name(ver, system, arch)
    url = url_template.format(tag=info["tag"], asset=name)
    out_dir = os.path.join(dest_root, "vina-%s" % ver)
    os.makedirs(out_dir, exist_ok=True)
    final = os.path.join(out_dir, name)
    part = final + ".part"
    try:
        with _request(url, timeout) as resp, open(part, "wb") as f:
            total = resp.headers.get("Content-Length")
            total = int(total) if total and total.isdigit() else None
            done = 0
            while True:
                if cancelled is not None and cancelled():
                    raise VinaError("Download cancelled.")
                block = resp.read(1 << 16)
                if not block:
                    break
                f.write(block)
                done += len(block)
                if progress is not None:
                    progress(done, total)
        if os.path.getsize(part) < 100_000:
            raise VinaError("Downloaded file is suspiciously small (%d bytes); not a Vina binary."
                            % os.path.getsize(part))
        digest = sha256_of(part)
        expected = info["digests"].get(name)
        if expected and expected.lower() != digest.lower():
            raise VinaError("Checksum mismatch for %s (expected %s, got %s). File discarded."
                            % (name, expected, digest))
        os.replace(part, final)
    except VinaError:
        _silent_remove(part)
        raise
    except Exception as err:
        _silent_remove(part)
        raise VinaError("Could not download %s: %s" % (url, err))
    mode = os.stat(final).st_mode
    os.chmod(final, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    with open(final + ".sha256", "w") as f:
        f.write("%s  %s\n" % (digest, name))
    ok, msg = check_vina(final)
    if not ok:
        raise VinaError("Downloaded file did not run: %s" % msg)
    return final, msg


def _silent_remove(path):
    try:
        os.remove(path)
    except OSError:
        pass
