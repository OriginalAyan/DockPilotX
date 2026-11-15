"""Session-level controller shared by the GUI tool and the ``vinapilotx`` commands.

Keeps persistent settings, owns the Vina executable, runs docking jobs in a
background thread and relays progress to the GUI on ChimeraX's main thread.
"""

import json
import os
import threading

from . import vina_manager, docking
from .docking import DockingParams, DockingRunner

APP_NAME = "VinaPilotX"


def data_dir():
    """Folder for settings and downloaded programs (per user)."""
    try:
        from chimerax import app_dirs
        base = app_dirs.user_data_dir
    except Exception:
        base = os.path.join(os.path.expanduser("~"), ".local", "share")
    path = os.path.join(base, APP_NAME)
    os.makedirs(path, exist_ok=True)
    return path


class Settings:
    """Tiny JSON settings store (vina path, last folders, last parameters)."""

    def __init__(self, path=None):
        self.path = path or os.path.join(data_dir(), "settings.json")
        self.data = {}
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                self.data = json.load(f)
        except (OSError, ValueError):
            self.data = {}

    def get(self, key, default=None):
        return self.data.get(key, default)

    def set(self, key, value):
        self.data[key] = value
        self.save()

    def save(self):
        try:
            tmp = self.path + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(self.data, f, indent=2)
            os.replace(tmp, self.path)
        except OSError:
            pass


class DockingController:
    """One per ChimeraX session. Listeners receive (event, payload) on the main thread.

    Events: 'vina' (ok, version_or_message), 'download_progress' (done, total),
    'job' (DockingJob), 'progress' (done, total), 'log' (str), 'finished' (runner),
    'busy' (bool).
    """

    def __init__(self, session):
        self.session = session
        self.settings = Settings()
        self.listeners = []
        self.runner = None
        self.vina_exe = None
        self.vina_version = None
        self.workdir = None
        self._download_cancel = threading.Event()
        self._downloading = False
        saved = self.settings.get("vina_path")
        if saved:
            self.set_vina(saved, quiet=True)
        else:
            found = vina_manager.find_vina_on_system([os.path.join(data_dir(), "vina")])
            if not found:
                found = self._find_downloaded()
            if found:
                self.set_vina(found, quiet=True)

    # ---- plumbing ---------------------------------------------------------
    def add_listener(self, fn):
        if fn not in self.listeners:
            self.listeners.append(fn)

    def remove_listener(self, fn):
        if fn in self.listeners:
            self.listeners.remove(fn)

    def _emit(self, event, payload=None):
        for fn in list(self.listeners):
            try:
                fn(event, payload)
            except RuntimeError:
                pass  # a Qt widget that was already deleted

    def _emit_threadsafe(self, event, payload=None):
        """Safe to call from any thread."""
        self.session.ui.thread_safe(self._emit, event, payload)

    def log(self, message, level="info"):
        getattr(self.session.logger, level)(message)
        self._emit("log", message)

    def _log_threadsafe(self, message):
        self.session.ui.thread_safe(self.log, message)

    # ---- Vina -------------------------------------------------------------
    def _find_downloaded(self):
        root = os.path.join(data_dir(), "vina")
        if not os.path.isdir(root):
            return None
        best = None
        for d in os.listdir(root):
            p = vina_manager.resolve_vina_path(os.path.join(root, d))
            if p:
                ver = vina_manager.version_tuple(vina_manager.parse_version(d) or "0")
                if best is None or ver > best[0]:
                    best = (ver, p)
        return best[1] if best else None

    def set_vina(self, path, quiet=False):
        """Verify and remember a Vina executable (or the folder containing it)."""
        ok, info = vina_manager.check_vina(path)
        if ok:
            self.vina_exe = vina_manager.resolve_vina_path(path)
            self.vina_version = info
            self.settings.set("vina_path", self.vina_exe)
            if not quiet:
                self.log("AutoDock Vina %s ready: %s" % (info, self.vina_exe))
        else:
            self.vina_exe = None
            self.vina_version = None
            if not quiet:
                self.log(info, "warning")
        self._emit("vina", (ok, info))
        return ok, info

    @property
    def downloading(self):
        return self._downloading

    def download_vina(self):
        """Download the newest Vina for this platform in a background thread."""
        if self._downloading:
            return
        self._downloading = True
        self._download_cancel.clear()
        self._emit("busy", True)
        dest = os.path.join(data_dir(), "vina")

        def work():
            try:
                self._log_threadsafe("Looking up the latest AutoDock Vina release ...")
                path, ver = vina_manager.download_vina(
                    dest, progress=lambda d, t: self._emit_threadsafe("download_progress", (d, t)),
                    cancelled=self._download_cancel.is_set)
                self.session.ui.thread_safe(self._download_done, path, None)
            except Exception as err:
                self.session.ui.thread_safe(self._download_done, None, str(err))

        threading.Thread(target=work, name="VinaPilotX-download", daemon=True).start()

    def cancel_download(self):
        self._download_cancel.set()

    def _download_done(self, path, error):
        self._downloading = False
        self._emit("busy", False)
        if error:
            self.log("Vina download failed: %s" % error, "warning")
            self._emit("vina", (False, error))
        else:
            self.set_vina(path)

    # ---- native PDBQT reader override ----------------------------------------
    @property
    def native_pdbqt_wanted(self):
        return bool(self.settings.get("native_pdbqt", True))

    def set_native_pdbqt(self, enabled, quiet=False):
        """Turn the 'open .pdbqt with ChimeraX's native reader' override on or off."""
        from . import format_preference as fp
        try:
            if enabled:
                fp.prefer_native_pdbqt(self.session.data_formats, force=True)
            else:
                fp.restore_pdbqt_preference(self.session.data_formats)
        except LookupError as err:
            self.log("Could not enable the native .pdbqt reader: %s" % err, "warning")
            return False
        self.settings.set("native_pdbqt", bool(enabled))
        if not quiet:
            self.log("Native .pdbqt opening %s." % ("ENABLED (plain 'open file.pdbqt' uses ChimeraX's AutoDock "
                                                  "PDBQT reader)" if enabled else
                                                  "disabled; ChimeraX's original .pdbqt behaviour is restored"))
        return True

    # ---- docking ----------------------------------------------------------
    @property
    def running(self):
        return self.runner is not None and self.runner._thread is not None \
            and self.runner._thread.is_alive()

    def start_docking(self, receptor, ligands, workdir, box_center, box_size, params,
                      recursive=False, overwrite=False):
        """Validate, stage and launch docking. Returns (ok, message_or_job_count)."""
        if self.running:
            return False, "A docking run is already in progress."
        msgs = docking.preflight(self.vina_exe, receptor, box_center, box_size, params, workdir)
        if not ligands:
            msgs.append(("error", "No ligand PDBQT files or folders selected."))
        for level, text in msgs:
            self.log(("Warning: " if level == "warning" else "Error: ") + text,
                     "warning")
        if any(level == "error" for level, _ in msgs):
            return False, "; ".join(t for l, t in msgs if l == "error")
        runner = DockingRunner(
            self.vina_exe, receptor, box_center, box_size, params, workdir, ligands,
            recursive=recursive, overwrite=overwrite,
            on_log=self._log_threadsafe,
            on_job=lambda job: self._emit_threadsafe("job", job),
            on_progress=lambda d, t: self._emit_threadsafe("progress", (d, t)),
            on_finished=self._runner_finished)
        try:
            jobs = runner.prepare()
        except Exception as err:
            return False, "Could not prepare the docking project: %s" % err
        if not jobs:
            return False, "No ligands found in the selected inputs."
        self.runner = runner
        self.workdir = os.path.abspath(workdir)
        self.settings.set("last_workdir", self.workdir)
        self.log("Vina %s | receptor %s | %d ligand(s) | box center (%.2f, %.2f, %.2f) size (%.2f, %.2f, %.2f)"
                 % (self.vina_version, os.path.basename(receptor), len(jobs), *box_center, *box_size))
        self._emit("busy", True)
        self._emit("started", jobs)
        runner.start()
        return True, len(jobs)

    def _runner_finished(self, runner):
        self.session.ui.thread_safe(self._finished_main, runner)

    def _finished_main(self, runner):
        counts = {}
        for j in runner.jobs:
            counts[j.status] = counts.get(j.status, 0) + 1
        summary = ", ".join("%d %s" % (n, s) for s, n in sorted(counts.items()))
        best = [j for j in runner.jobs if j.best is not None]
        msg = "Docking %s: %s." % ("cancelled" if runner.cancelled else "finished", summary or "no jobs")
        if best:
            top = min(best, key=lambda j: j.best)
            msg += " Best: %s (%.2f kcal/mol)." % (top.name, top.best)
        msg += " Results in %s" % runner.workdir
        self.log(msg)
        self._emit("busy", False)
        self._emit("finished", runner)

    def cancel_docking(self):
        if self.running:
            self.runner.cancel()
            self.log("Cancelling docking ...")


def get_controller(session):
    ctrl = getattr(session, "_vinapilotx_controller", None)
    if ctrl is None:
        ctrl = DockingController(session)
        session._vinapilotx_controller = ctrl
    return ctrl
