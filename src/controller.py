"""Session-level controller shared by the GUI tool and the ``vinapilotx`` commands.

Keeps persistent settings, owns the Vina executable, runs docking jobs in a
background thread and relays progress to the GUI on ChimeraX's main thread.
"""

import json
import os
import threading

import shutil

from . import vina_manager, docking, prep, prep_env
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
        self._prep_cancel = threading.Event()
        self.prep_env = prep_env.PrepEnv.from_dict(self.settings.get("prep_env"))
        self.prep_info = self.settings.get("prep_info") or {}
        self.prep_message = ""
        self.prep_ok = bool(self.prep_env and self.prep_env.exists() and self.prep_info.get("meeko"))
        self.prep_busy = False
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

    # ---- Meeko preparation environment ---------------------------------------
    @property
    def prep_dir(self):
        return os.path.join(data_dir(), "prep")

    def _prep_log(self, line):
        self.session.ui.thread_safe(self._emit, "prep_log", line)

    def _prep_set_busy(self, busy):
        self.prep_busy = busy
        self._emit("prep_busy", busy)

    def _prep_apply(self, result):
        self.prep_ok = bool(result.get("ok"))
        self.prep_info = result.get("info", {})
        self.prep_message = result.get("message", "")
        self.settings.set("prep_info", self.prep_info if self.prep_ok else {})
        self._emit("prep_env", result)

    def cancel_prep(self):
        self._prep_cancel.set()

    def _prep_thread(self, work, done):
        """Run work() in a thread; done(result, error) runs on the main thread."""
        if self.prep_busy:
            self.log("Another preparation task is still running.", "warning")
            return False
        self._prep_cancel.clear()
        self._prep_set_busy(True)

        def run():
            try:
                res, err = work(), None
            except prep_env.PrepError as e:
                res, err = None, str(e)
            except Exception as e:  # never let a thread die silently
                res, err = None, "%s: %s" % (type(e).__name__, e)

            def finish():
                self._prep_set_busy(False)
                done(res, err)
            self.session.ui.thread_safe(finish)
        threading.Thread(target=run, daemon=True).start()
        return True

    def verify_prep(self):
        """Check the saved environment (runs in the background)."""
        env = self.prep_env

        def work():
            return prep_env.probe(env)

        def done(res, err):
            self._prep_apply(res if res else {"ok": False, "message": err or "failed", "info": {}})
        return self._prep_thread(work, done)

    def setup_prep_env(self, method="auto"):
        """Create the Meeko environment: method 'auto' (micromamba, then venv), 'micromamba' or 'venv'."""
        base_dir = self.prep_dir

        def work():
            log, cancelled = self._prep_log, self._prep_cancel.is_set
            env, errors, mm_base = None, [], None
            if method in ("auto", "micromamba"):
                try:
                    mm_base = prep_env.short_base_dir(os.path.join(base_dir, "micromamba"))
                    env = prep_env.create_micromamba_env(mm_base, log, cancelled)
                    res = prep_env.probe(env)
                    if not res["ok"]:
                        raise prep_env.PrepError(res["message"])
                except prep_env.PrepError as err:
                    if mm_base:
                        shutil.rmtree(mm_base, ignore_errors=True)
                    if str(err) == "Cancelled.":
                        raise
                    errors.append("micromamba: %s" % err)
                    env = None
                    log("micromamba route failed: %s" % err)
                    if method == "auto":
                        log("Trying the venv + pip fallback ...")
            if env is None and method in ("auto", "venv"):
                env = prep_env.create_venv_env(base_dir, log, cancelled)
            if env is None:
                raise prep_env.PrepError("; ".join(errors) or "Setup failed.")
            return env, prep_env.probe(env)

        def done(res, err):
            if err:
                self._prep_apply({"ok": False, "message": err, "info": {}})
                self.log("Meeko environment setup failed: %s" % err, "warning")
                return
            env, result = res
            if result["ok"]:
                self.prep_env = env
                self.settings.set("prep_env", env.to_dict())
            self._prep_apply(result)
        return self._prep_thread(work, done)

    def use_external_python(self, path):
        env = prep_env.external_env(path)
        if env is None:
            self._prep_apply({"ok": False, "message": "No Python found at: %s" % path, "info": {}})
            return False

        def work():
            return env, prep_env.probe(env)

        def done(res, err):
            if err:
                self._prep_apply({"ok": False, "message": err, "info": {}})
                return
            env_, result = res
            if result["ok"]:
                self.prep_env = env_
                self.settings.set("prep_env", env_.to_dict())
            self._prep_apply(result)
        return self._prep_thread(work, done)

    def remove_prep_env(self):
        env = self.prep_env
        target = None
        if env is not None and env.prefix:
            if env.kind == "micromamba" and os.path.basename(os.path.dirname(env.prefix)) in ("mm", "micromamba"):
                target = os.path.dirname(env.prefix)
            elif env.kind == "venv" and os.path.basename(env.prefix) == "venv":
                target = env.prefix
        if target:
            shutil.rmtree(target, ignore_errors=True)
        self.prep_env, self.prep_info = None, {}
        self.settings.set("prep_env", None)
        self.settings.set("prep_info", {})
        self._prep_apply({"ok": False, "message": "No preparation environment is set up.", "info": {}})

    def prepare_receptor_async(self, clean_pdb, out_base, delete_bad=True):
        env, info = self.prep_env, self.prep_info

        def work():
            return prep.prepare_receptor(env, info, clean_pdb, out_base, delete_bad,
                                         log=self._prep_log, cancelled=self._prep_cancel.is_set)

        def done(res, err):
            self._emit("prep_receptor_done", (res, err))
        return self._prep_thread(work, done)

    def prepare_ligands_async(self, inputs, out_dir, stage_dir, add_h=True):
        env, info = self.prep_env, self.prep_info

        def progress(d, t):
            self.session.ui.thread_safe(self._emit, "prep_progress", (d, t))

        def work():
            return prep.prepare_ligands(env, info, inputs, out_dir, stage_dir, add_h, log=self._prep_log,
                                        progress=progress, cancelled=self._prep_cancel.is_set)

        def done(res, err):
            self._emit("prep_ligands_done", (res, err))
        return self._prep_thread(work, done)

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
