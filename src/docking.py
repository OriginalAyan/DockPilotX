"""Run AutoDock Vina for one or many ligands, with resume, cancel and a manifest.

Pure Python (no ChimeraX imports). Callbacks are invoked from worker threads;
the GUI layer is responsible for marshalling them to the main thread.
"""

import csv
import datetime
import json
import os
import platform
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, asdict

from .vina_manager import no_window_flags, check_vina
from .pdbqt_tools import check_receptor_pdbqt, scan_ligands, stage_records, LigandRecord

TOOL_VERSION = "0.6.2"
SCORING_FUNCTIONS = ("vina", "vinardo")
VINA_VOLUME_WARN = 27000.0  # Å^3, Vina's own warning threshold

STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_DONE = "done"
STATUS_EXISTING = "done (existing)"
STATUS_FAILED = "failed"
STATUS_CANCELLED = "cancelled"
STATUS_SKIPPED = "skipped"

_RESULT_RE = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")


@dataclass
class DockingParams:
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    cpu: int = 0            # CPUs per Vina process; 0 = automatic
    parallel_jobs: int = 1  # ligands docked at the same time
    seed: int = 0           # 0 = let Vina choose (its reported seed is recorded); >0 = fixed seed
    scoring: str = "vina"

    def validate(self):
        errs = []
        if not 1 <= int(self.exhaustiveness) <= 1024:
            errs.append("exhaustiveness must be between 1 and 1024")
        if not 1 <= int(self.num_modes) <= 100:
            errs.append("num_modes must be between 1 and 100")
        if not 0 < float(self.energy_range) <= 100:
            errs.append("energy_range must be > 0 and <= 100 kcal/mol")
        if int(self.cpu) < 0:
            errs.append("cpu cannot be negative (0 = automatic)")
        if int(self.parallel_jobs) < 1:
            errs.append("parallel jobs must be at least 1")
        if int(self.seed) < 0:
            errs.append("seed cannot be negative (0 = random)")
        if self.scoring not in SCORING_FUNCTIONS:
            errs.append("scoring must be one of: %s" % ", ".join(SCORING_FUNCTIONS))
        return errs

    def cpu_per_job(self):
        if int(self.cpu) > 0:
            return int(self.cpu)
        if int(self.parallel_jobs) > 1:
            return max(1, (os.cpu_count() or 1) // int(self.parallel_jobs))
        return 0


@dataclass
class DockingJob:
    name: str
    ligand_path: str
    out_path: str
    log_path: str
    status: str = STATUS_PENDING
    returncode: int = None
    poses: list = field(default_factory=list)  # [(affinity, rmsd_lb, rmsd_ub), ...]; pose number = index + 1
    seconds: float = 0.0
    seed: int = None
    error: str = ""
    progress: float = 0.0  # 0..100 for the running job

    @property
    def best(self):
        return min(p[0] for p in self.poses) if self.poses else None

    @property
    def n_poses(self):
        return len(self.poses)


# --------------------------------------------------------------------------
def parse_vina_out(path):
    """Return a list of (affinity, rmsd_lb, rmsd_ub) for each pose in a Vina output file."""
    poses = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                m = _RESULT_RE.match(line)
                if m:
                    poses.append(tuple(float(x) for x in m.groups()))
    except OSError:
        pass
    return poses


def box_volume(size):
    return float(size[0]) * float(size[1]) * float(size[2])


def write_vina_config(path, center, size, params):
    """Write the shared Vina configuration file (box + search parameters).

    File paths are deliberately NOT written here: they are passed on the command
    line as separate arguments, which is safe for paths containing spaces.
    """
    lines = [
        "# AutoDock Vina configuration written by VinaPilotX",
        "center_x = %.4f" % center[0],
        "center_y = %.4f" % center[1],
        "center_z = %.4f" % center[2],
        "size_x = %.4f" % size[0],
        "size_y = %.4f" % size[1],
        "size_z = %.4f" % size[2],
        "exhaustiveness = %d" % int(params.exhaustiveness),
        "num_modes = %d" % int(params.num_modes),
        "energy_range = %s" % ("%g" % float(params.energy_range)),
        "scoring = %s" % params.scoring,
    ]
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")


def build_command(vina_exe, config_path, receptor, ligand, out_path, cpu=0, seed=None):
    cmd = [vina_exe, "--config", config_path, "--receptor", receptor,
           "--ligand", ligand, "--out", out_path]
    if seed is not None and int(seed) > 0:
        cmd += ["--seed", str(int(seed))]
    if cpu and int(cpu) > 0:
        cmd += ["--cpu", str(int(cpu))]
    return cmd


def preflight(vina_exe, receptor, box_center, box_size, params, workdir):
    """Check everything that can be checked before docking starts.

    Returns a list of (level, message) with level in {'error', 'warning'}.
    """
    msgs = []
    if not vina_exe:
        msgs.append(("error", "AutoDock Vina is not set up. Choose vina.exe or download it on the Settings tab."))
    else:
        ok, info = check_vina(vina_exe)
        if not ok:
            msgs.append(("error", info))
    if not receptor:
        msgs.append(("error", "No receptor PDBQT selected."))
    else:
        msgs += check_receptor_pdbqt(receptor)
    if box_center is None or box_size is None:
        msgs.append(("error", "No search box defined. Create one on the Box tab."))
    else:
        if any(float(s) <= 0 for s in box_size):
            msgs.append(("error", "Box sizes must all be > 0 Å."))
        elif box_volume(box_size) > VINA_VOLUME_WARN:
            msgs.append(("warning", "Search box volume is %.0f Å^3 (> %.0f). Vina warns that very large "
                                    "boxes need higher exhaustiveness to stay reliable."
                         % (box_volume(box_size), VINA_VOLUME_WARN)))
    for e in params.validate():
        msgs.append(("error", e))
    if not workdir:
        msgs.append(("error", "No working directory selected."))
    if sys.platform.startswith("win"):
        for label, p in (("Vina", vina_exe), ("receptor", receptor), ("working directory", workdir)):
            if p and not _is_ascii(p):
                msgs.append(("warning", "The %s path contains non-ASCII characters; Vina on Windows "
                                        "may fail to open it. Use a plain path such as C:/docking." % label))
    return msgs


def _is_ascii(text):
    try:
        str(text).encode("ascii")
        return True
    except UnicodeEncodeError:
        return False


# --------------------------------------------------------------------------
class DockingRunner:
    """Docks a list of ligands against one receptor with a shared box.

    Layout created inside *workdir*::

        receptor/<name>.pdbqt      copy of the receptor used
        ligands/<name>.pdbqt       one file per ligand (libraries are split)
        results/<name>_out.pdbqt   Vina poses
        logs/<name>.log            Vina console output
        vina_config.txt            box + search parameters
        manifest.json              versions, parameters, commands, timings
        summary.csv                one row per ligand
    """

    def __init__(self, vina_exe, receptor, box_center, box_size, params, workdir,
                 ligands, recursive=False, overwrite=False,
                 on_log=None, on_job=None, on_progress=None, on_finished=None):
        self.vina_exe = vina_exe
        self.receptor_src = os.path.abspath(receptor)
        self.box_center = tuple(float(v) for v in box_center)
        self.box_size = tuple(float(v) for v in box_size)
        self.params = params
        self.workdir = os.path.abspath(workdir)
        self.ligands = list(ligands)  # paths (files/folders) or LigandRecord objects
        self.recursive = recursive
        self.overwrite = overwrite
        self._user_on_log = on_log or (lambda msg: None)
        self._log_lock = threading.Lock()
        self.on_job = on_job or (lambda job: None)
        self.on_progress = on_progress or (lambda done, total: None)
        self.on_finished = on_finished or (lambda runner: None)

        self.jobs = []
        self.receptor = None
        self.config_path = None
        self.vina_version = None
        self.started = None
        self.ended = None
        self._cancel = threading.Event()
        self._lock = threading.Lock()
        self._procs = set()
        self._thread = None
        self._done_count = 0

    def log(self, message):
        """Show a message (GUI/log window) and append it to run_log.txt in the working folder."""
        try:
            os.makedirs(self.workdir, exist_ok=True)
            with self._log_lock, open(os.path.join(self.workdir, "run_log.txt"), "a", encoding="utf-8") as f:
                f.write("%s  %s\n" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message))
        except OSError:
            pass
        self._user_on_log(message)

    # ---- paths -----------------------------------------------------------
    @property
    def results_dir(self):
        return os.path.join(self.workdir, "results")

    @property
    def logs_dir(self):
        return os.path.join(self.workdir, "logs")

    @property
    def ligands_dir(self):
        return os.path.join(self.workdir, "ligands")

    # ---- preparation -----------------------------------------------------
    def prepare(self):
        """Create folders, stage receptor + ligands, write config. Returns job list."""
        import shutil
        for d in ("receptor", "ligands", "results", "logs"):
            os.makedirs(os.path.join(self.workdir, d), exist_ok=True)
        rec_dst = os.path.join(self.workdir, "receptor", os.path.basename(self.receptor_src))
        if os.path.abspath(rec_dst) != self.receptor_src:
            shutil.copyfile(self.receptor_src, rec_dst)
        self.receptor = rec_dst
        self.config_path = os.path.join(self.workdir, "vina_config.txt")
        write_vina_config(self.config_path, self.box_center, self.box_size, self.params)
        ok, ver = check_vina(self.vina_exe)
        self.vina_version = ver if ok else None
        if self.ligands and all(isinstance(l, LigandRecord) for l in self.ligands):
            records = self.ligands
        else:
            records = scan_ligands(self.ligands, self.recursive)
        items = stage_records(records, self.ligands_dir)
        self.jobs = []
        for it in items:
            job = DockingJob(
                name=it.name, ligand_path=it.path,
                out_path=os.path.join(self.results_dir, it.name + "_out.pdbqt"),
                log_path=os.path.join(self.logs_dir, it.name + ".log"))
            if it.problems:
                job.status = STATUS_SKIPPED
                job.error = "Invalid ligand PDBQT: " + "; ".join(it.problems)
            self.jobs.append(job)
        return self.jobs

    # ---- execution -------------------------------------------------------
    def start(self):
        """Run in a background thread; returns immediately."""
        self._thread = threading.Thread(target=self.run, name="VinaPilotX-runner", daemon=True)
        self._thread.start()

    def join(self, timeout=None):
        if self._thread is not None:
            self._thread.join(timeout)

    def cancel(self):
        self._cancel.set()
        with self._lock:
            procs = list(self._procs)
        for p in procs:
            try:
                p.terminate()
            except OSError:
                pass

    @property
    def cancelled(self):
        return self._cancel.is_set()

    def run(self):
        """Blocking run of all jobs (called by start() in its own thread)."""
        self.started = datetime.datetime.now()
        try:
            self.log("Run started: Vina %s | receptor %s | box center (%.2f, %.2f, %.2f) size (%.2f, %.2f, %.2f) | "
                     "exhaustiveness %d, poses %d, scoring %s" % (
                         self.vina_version, os.path.basename(self.receptor_src), *self.box_center, *self.box_size,
                         self.params.exhaustiveness, self.params.num_modes, self.params.scoring))
            if not self.jobs:
                self.prepare()
            total = len(self.jobs)
            self._write_manifest(partial=True)
            runnable = []
            for job in self.jobs:
                if job.status == STATUS_SKIPPED:
                    self.log("Skipping %s: %s" % (job.name, job.error))
                    self._finish_job(job, total)
                elif not self.overwrite and self._load_existing(job):
                    self.log("%s: using existing result (%.2f kcal/mol)" % (job.name, job.best))
                    self._finish_job(job, total)
                else:
                    runnable.append(job)
            workers = max(1, min(int(self.params.parallel_jobs), len(runnable) or 1))
            self.log("Docking %d ligand(s) with %d parallel job(s)." % (len(runnable), workers))
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(self._run_job, job, total) for job in runnable]
                for fut in futures:
                    try:
                        fut.result()
                    except Exception as err:  # never let one job kill the batch
                        self.log("Internal error: %s" % err)
            for job in self.jobs:
                if job.status in (STATUS_PENDING, STATUS_RUNNING):
                    job.status = STATUS_CANCELLED
        finally:
            self.ended = datetime.datetime.now()
            try:
                self.write_summary()
                self._write_manifest(partial=False)
            except OSError as err:
                self.log("Could not write summary files: %s" % err)
            self.on_finished(self)

    def _load_existing(self, job):
        if os.path.isfile(job.out_path) and os.path.isfile(job.log_path):
            poses = parse_vina_out(job.out_path)
            if poses:
                job.poses = list(poses)
                job.seed = seed_from_log(job.log_path)
                job.status = STATUS_EXISTING
                return True
        return False

    def _finish_job(self, job, total):
        with self._lock:
            self._done_count += 1
            done = self._done_count
        self.on_job(job)
        self.on_progress(done, total)

    def _run_job(self, job, total):
        if self._cancel.is_set():
            job.status = STATUS_CANCELLED
            self._finish_job(job, total)
            return
        job.seed = None  # filled in from Vina's own output ("random seed: N")
        cmd = build_command(self.vina_exe, self.config_path, self.receptor, job.ligand_path,
                            job.out_path, self.params.cpu_per_job(),
                            int(self.params.seed) if int(self.params.seed) > 0 else None)
        job.status = STATUS_RUNNING
        job.progress = 0.0
        self.on_job(job)
        self.log("Docking %s ..." % job.name)
        t0 = time.monotonic()
        tail = []
        head = b""
        stars = 0
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    cwd=self.workdir, creationflags=no_window_flags())
        except OSError as err:
            job.status = STATUS_FAILED
            job.error = "Could not start Vina: %s" % err
            self.log("%s: %s" % (job.name, job.error))
            self._finish_job(job, total)
            return
        with self._lock:
            self._procs.add(proc)
        try:
            with open(job.log_path, "wb") as log:
                log.write(("# command: %s\n" % " ".join(_q(c) for c in cmd)).encode("utf-8"))
                buf = b""
                while True:
                    chunk = proc.stdout.read1(4096) if hasattr(proc.stdout, "read1") else proc.stdout.read(4096)
                    if not chunk:
                        break
                    log.write(chunk)
                    if job.seed is None and len(head) < 16384:
                        head += chunk
                        seed = parse_seed(head.decode("utf-8", "replace"))
                        if seed is not None:
                            job.seed = seed
                            self.on_job(job)
                    stars += chunk.count(b"*")
                    buf += chunk
                    parts = re.split(rb"[\r\n]", buf)
                    buf = parts.pop()
                    for p in parts:
                        txt = p.decode("utf-8", "replace").strip()
                        if txt and set(txt) - set("*|-% 0123456789"):
                            tail.append(txt)
                            del tail[:-8]
                    new_progress = min(100.0, stars / 51.0 * 100.0)
                    if new_progress - job.progress >= 2.0:
                        job.progress = new_progress
                        self.on_job(job)
            proc.wait()
        finally:
            with self._lock:
                self._procs.discard(proc)
        job.seconds = time.monotonic() - t0
        job.returncode = proc.returncode
        job.progress = 100.0 if proc.returncode == 0 else job.progress
        if self._cancel.is_set() and proc.returncode != 0:
            job.status = STATUS_CANCELLED
        else:
            poses = parse_vina_out(job.out_path) if proc.returncode == 0 else []
            if proc.returncode == 0 and poses:
                job.poses = list(poses)
                job.status = STATUS_DONE
                self.log("%s: best %.2f kcal/mol (%d poses, %.1f s, Vina seed %s)"
                            % (job.name, job.best, job.n_poses, job.seconds,
                               job.seed if job.seed is not None else "not reported"))
            else:
                job.status = STATUS_FAILED
                job.error = _summarize_error(tail, proc.returncode)
                self.log("%s FAILED: %s" % (job.name, job.error))
        self._finish_job(job, total)

    # ---- output files ----------------------------------------------------
    def write_summary(self):
        """summary.csv: one row per POSE (all poses of every ligand); best_per_ligand.csv: one row per ligand."""
        path = os.path.join(self.workdir, "summary.csv")
        write_pose_csv(self.jobs, path)
        best = os.path.join(self.workdir, "best_per_ligand.csv")
        with open(best, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ligand", "status", "best_affinity_kcal_mol", "n_poses", "seed", "output_file", "error"])
            for j in self.jobs:
                w.writerow([j.name, j.status, "" if j.best is None else "%.3f" % j.best, j.n_poses,
                            "" if j.seed is None else j.seed,
                            os.path.relpath(j.out_path, self.workdir) if os.path.exists(j.out_path) else "",
                            j.error])
        return path

    def _write_manifest(self, partial):
        first = self.jobs[0] if self.jobs else None
        example = ""
        if first is not None and self.receptor:
            example = " ".join(_q(c) for c in build_command(
                self.vina_exe, self.config_path, self.receptor, first.ligand_path,
                first.out_path, self.params.cpu_per_job(),
                int(self.params.seed) if int(self.params.seed) > 0 else None))
        counts = {}
        for j in self.jobs:
            counts[j.status] = counts.get(j.status, 0) + 1
        data = {
            "tool": "VinaPilotX", "tool_version": TOOL_VERSION, "complete": not partial,
            "started": self.started.isoformat(timespec="seconds") if self.started else None,
            "ended": self.ended.isoformat(timespec="seconds") if self.ended else None,
            "vina_executable": self.vina_exe, "vina_version": self.vina_version,
            "receptor_source": self.receptor_src,
            "receptor_used": self.receptor,
            "box": {"center": self.box_center, "size": self.box_size},
            "parameters": asdict(self.params),
            "seed_policy": ("fixed seed %d passed to Vina for every ligand" % self.params.seed)
                           if int(self.params.seed) > 0
                           else "Vina chooses a random seed per ligand; the seed it reported is recorded below",
            "seeds": {j.name: j.seed for j in self.jobs if j.seed is not None},
            "example_command": example,
            "n_ligands": len(self.jobs), "status_counts": counts,
            "python": sys.version.split()[0], "platform": platform.platform(),
        }
        with open(os.path.join(self.workdir, "manifest.json"), "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)


def _q(arg):
    return '"%s"' % arg if (" " in arg or "\t" in arg) else arg


def _summarize_error(tail, returncode):
    useful = [t for t in tail if t.lower().startswith(("error", "parse error", "pdbqt", "atom",
                                                       "unknown", "warning", "terminate"))]
    text = " | ".join((useful or tail)[-3:])
    if not text:
        text = "Vina exited with code %s and produced no pose file." % returncode
    return text[:400]


_VINA_SEED_RE = re.compile(r"random seed:\s*(-?\d+)", re.IGNORECASE)
_CMD_SEED_RE = re.compile(r"--seed\s+(\d+)")


def parse_seed(text):
    """The seed Vina reports, e.g. 'Performing docking (random seed: -1805814086) ...'."""
    m = _VINA_SEED_RE.search(text)
    return int(m.group(1)) if m else None


def seed_from_log(log_path):
    """Seed Vina reported in a saved log (falls back to '--seed N' in our command header)."""
    try:
        with open(log_path, "r", errors="replace") as f:
            head = f.read(16384)
    except OSError:
        return None
    seed = parse_seed(head)
    if seed is None:
        m = _CMD_SEED_RE.search(head.split("\n", 1)[0])
        seed = int(m.group(1)) if m else None
    return seed


def load_results_dir(workdir):
    """Rebuild DockingJob objects from an existing working directory (for the Results tab)."""
    jobs = []
    rdir = os.path.join(workdir, "results")
    if not os.path.isdir(rdir):
        return jobs
    for fn in sorted(os.listdir(rdir), key=str.lower):
        if not fn.endswith("_out.pdbqt"):
            continue
        name = fn[:-len("_out.pdbqt")]
        out = os.path.join(rdir, fn)
        poses = parse_vina_out(out)
        job = DockingJob(name=name, ligand_path=os.path.join(workdir, "ligands", name + ".pdbqt"),
                         out_path=out, log_path=os.path.join(workdir, "logs", name + ".log"))
        if poses:
            job.poses = list(poses)
            job.seed = seed_from_log(job.log_path)
            job.status = STATUS_EXISTING
            job.progress = 100.0
        else:
            job.status = STATUS_FAILED
            job.error = "No poses found in output file."
        jobs.append(job)
    return jobs


POSE_CSV_HEADER = ["ligand", "pose", "affinity_kcal_mol", "rmsd_lb", "rmsd_ub", "rating", "hbonds", "hbonds_relaxed", "clashes",
                   "hbond_partners", "seed", "status", "error"]


def write_pose_csv(jobs, path, analysis=None, only=None):
    """Write one row per pose. *analysis* (AnalysisStore) adds rating / H-bonds / clashes.

    *only*: optional callable(job_name, pose_number, record_dict) -> bool to filter rows.
    Ligands without poses (failed/skipped) get a single row carrying status and error.
    """
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(POSE_CSV_HEADER)
        for j in jobs:
            seed = "" if j.seed is None else j.seed
            if not j.poses:
                if only is None:
                    w.writerow([j.name, "", "", "", "", "", "", "", "", "", seed, j.status, j.error])
                continue
            for n, (aff, lb, ub) in enumerate(j.poses, start=1):
                rec = analysis.get(j.name, n) if analysis is not None else {}
                if only is not None and not only(j.name, n, rec):
                    continue
                w.writerow([j.name, n, "%.3f" % aff, "%.3f" % lb, "%.3f" % ub,
                            rec.get("rating") or "", "" if rec.get("hbonds") is None else rec["hbonds"],
                            "" if rec.get("hbonds_relaxed") is None else rec["hbonds_relaxed"],
                            "" if rec.get("clashes") is None else rec["clashes"],
                            "; ".join(rec.get("hbond_partners", [])), seed, j.status, j.error])


def extract_pose(out_path, pose_number, dest_path):
    """Write pose *pose_number* (1-based) of a Vina output file as a stand-alone single-model PDBQT."""
    current, lines, found = 0, [], False
    with open(out_path, "r", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            rec = line.split(None, 1)[0] if line.strip() else ""
            if rec == "MODEL":
                current += 1
                lines = []
                continue
            if rec == "ENDMDL":
                if current == pose_number:
                    found = True
                    break
                continue
            if current == pose_number:
                lines.append(line)
    if not found:
        raise ValueError("Pose %d not found in %s" % (pose_number, out_path))
    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    with open(dest_path, "w", newline="\n") as out:
        out.write("\n".join(lines) + "\n")
    return dest_path
