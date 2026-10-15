"""Minimal stand-in for ecalautoctrl's HandlerBase / HCALCondorHandler.

State lives in <outdir>/status.json:

    {"<task>": "<status>",
     "commit": ..., "cluster": ...,
     "jobs": {"<id>": {"cluster": ..., "status": ..., "exit": ..., "wall": ...}}}

with the Automator status words new, processing, done, reprocess, incomplete.

One HTCondor job is submitted per work item. The code is shipped with the job
as a `git archive` tarball of HEAD, so a run always corresponds to a commit.
Each job writes <outdir>/jobs/<id>.json when it ends (the stand-in for
`ecalautomation.py jobctrl --done/--failed`); the task is done when every job
reports done.

Runs with the lxplus system python3 (3.9): HTCondor is driven through
condor_submit / condor_q, since the htcondor bindings are not installed there.
"""
import glob
import json
import os
import re
import shutil
import subprocess
import tarfile
import time
from datetime import datetime, timezone

NEW = "new"
PROCESSING = "processing"
DONE = "done"
REPROCESS = "reprocess"
INCOMPLETE = "incomplete"

HELD = 5  # HTCondor JobStatus


def _now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _run(cmd, **kwargs):
    return subprocess.run(cmd, check=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, universal_newlines=True,
                          **kwargs).stdout


class HandlerBase:
    task = None         # task name, set by subclasses
    script = None       # job script, relative to the repository root

    def __init__(self, outdir, condor_dir, repo_dir, dry_run=False,
                 flavour="workday", memory=4000):
        self.outdir = outdir
        self.condor_dir = condor_dir
        self.repo_dir = repo_dir
        self.dry_run = dry_run
        self.flavour = flavour
        self.memory = memory
        self.status_file = os.path.join(outdir, "status.json")

    # --- to be provided by subclasses --------------------------------------

    def getWorkItems(self):
        """Return (items, extra): items is a list of {"arguments": str}."""
        raise NotImplementedError

    def inputFiles(self):
        """Extra files transferred with every job (besides the code tarball)."""
        return []

    def onDoneAction(self, st):
        """Called once when all jobs are done: copy the condor logs to EOS."""
        logdir = os.path.join(self.outdir, "logs")
        os.makedirs(logdir, exist_ok=True)
        cluster = st["cluster"]
        for path in glob.glob(os.path.join(self.condor_dir, "*.%s.*" % cluster)):
            shutil.copy2(path, logdir)

    # --- state --------------------------------------------------------------

    def load_status(self):
        if os.path.exists(self.status_file):
            with open(self.status_file) as f:
                return json.load(f)
        return {self.task: NEW, "jobs": {}}

    def save_status(self, st):
        os.makedirs(self.outdir, exist_ok=True)
        tmp = self.status_file + ".tmp"
        with open(tmp, "w") as f:
            json.dump(st, f, indent=2, sort_keys=True)
        os.replace(tmp, self.status_file)

    # --- code delivery ------------------------------------------------------

    def make_tarball(self):
        """Archive the code into condor_dir; return (commit, tarball path).

        Refuses if tracked files are modified, so every run is reproducible
        from its commit. A dry run accepts a modified tree and archives the
        working tree instead (tracked + untracked, not ignored), labelled
        "<commit>-dirty".
        """
        git = ["git", "-C", self.repo_dir]
        head = _run(git + ["rev-parse", "HEAD"]).strip()
        modified = _run(git + ["status", "--porcelain", "--untracked-files=no"]).strip()
        os.makedirs(self.condor_dir, exist_ok=True)
        tarball = os.path.join(self.condor_dir, "isotrack.tar.gz")
        if not modified:
            _run(git + ["archive", "--format=tar.gz", "-o", tarball, "HEAD"])
            return head, tarball
        if not self.dry_run:
            raise RuntimeError("tracked files are modified in %s; commit them first:\n%s"
                               % (self.repo_dir, modified))
        print("WARNING: dry run with modified tracked files; archiving the working tree")
        files = _run(git + ["ls-files", "-c", "-o", "--exclude-standard", "-z"]).split("\0")
        with tarfile.open(tarball, "w:gz") as tar:
            for name in sorted(f for f in files if f):
                tar.add(os.path.join(self.repo_dir, name), arcname=name)
        return head + "-dirty", tarball

    # --- HTCondor -----------------------------------------------------------

    def write_submit(self, items, commit, tarball):
        script = os.path.join(self.condor_dir, os.path.basename(self.script))
        shutil.copy2(os.path.join(self.repo_dir, self.script), script)
        inputs = [tarball] + [os.path.abspath(f) for f in self.inputFiles()]
        lines = [
            "# %s, written %s" % (self.task, _now()),
            "universe                = vanilla",
            "executable              = %s" % script,
            "initialdir              = %s" % self.condor_dir,
            "should_transfer_files   = YES",
            "when_to_transfer_output = ON_EXIT",
            "transfer_input_files    = %s" % ", ".join(inputs),
            "transfer_output_files   = job_$(ProcId).json",
            'environment             = "ISOTRACK_COMMIT=%s"' % commit,
            "output                  = job_$(ProcId).$(ClusterId).out",
            "error                   = job_$(ProcId).$(ClusterId).err",
            "log                     = job.$(ClusterId).log",
            '+JobFlavour             = "%s"' % self.flavour,
            "request_memory          = %d" % self.memory,
            "",
        ]
        for item in items:
            # Automator argument order: TASK JOB_ID EOS_DIR <item arguments>
            lines.append('arguments = "%s $(ProcId) %s %s"'
                         % (self.task, self.outdir, item["arguments"]))
            lines.append("queue 1")
        path = os.path.join(self.condor_dir, "submit.sub")
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return path

    def condor_submit(self, subfile):
        out = _run(["condor_submit", os.path.basename(subfile)], cwd=self.condor_dir)
        m = re.search(r"submitted to cluster (\d+)", out)
        if not m:
            raise RuntimeError("could not parse condor_submit output:\n" + out)
        return int(m.group(1))

    def queue_state(self, cluster):
        """Jobs of the cluster still in the queue: {proc_id: (status, hold reason)}."""
        out = _run(["condor_q", str(cluster), "-af:t", "ProcId", "JobStatus", "HoldReason"])
        jobs = {}
        for line in out.splitlines():
            fields = line.split("\t")
            if len(fields) >= 2:
                jobs[fields[0]] = (int(fields[1]), fields[2] if len(fields) > 2 else "")
        return jobs

    # --- workflow -----------------------------------------------------------

    def submit(self, wait=True, poll=300, st=None):
        st = st or self.load_status()
        if st[self.task] not in (NEW, REPROCESS):
            print("%s is '%s'; nothing to submit" % (self.task, st[self.task]))
            return st
        items, _ = self.getWorkItems()
        commit, tarball = self.make_tarball()
        subfile = self.write_submit(items, commit, tarball)
        if self.dry_run:
            print("Dry run: wrote %s and %s (commit %s); not submitted.\n"
                  % (tarball, subfile, commit))
            with open(subfile) as f:
                print(f.read())
            return st
        self._archive_job_results(st)
        cluster = self.condor_submit(subfile)
        print("Submitted %d job(s) to cluster %d" % (len(items), cluster))
        st[self.task] = PROCESSING
        st.update(commit=commit, cluster=cluster, submitted=_now())
        st["jobs"] = {str(i): {"cluster": cluster, "status": PROCESSING,
                               "exit": None, "wall": None}
                      for i in range(len(items))}
        self.save_status(st)
        return self.wait(poll) if wait else st

    def _archive_job_results(self, st):
        """Move job results of an earlier attempt aside, so they are not reread."""
        jobdir = os.path.join(self.outdir, "jobs")
        old = glob.glob(os.path.join(jobdir, "*.json"))
        if old:
            dest = os.path.join(jobdir, "old", str(st.get("cluster", "unknown")))
            os.makedirs(dest, exist_ok=True)
            for path in old:
                shutil.move(path, dest)

    def check(self):
        """Advance a processing task once its jobs have left the queue."""
        st = self.load_status()
        if st[self.task] != PROCESSING:
            return st
        queued = self.queue_state(st["cluster"])
        if queued:
            for jid, (status, reason) in queued.items():
                if jid in st["jobs"]:
                    st["jobs"][jid]["status"] = "held" if status == HELD else PROCESSING
                if status == HELD:
                    print("job %s.%s is held: %s" % (st["cluster"], jid, reason))
            self.save_status(st)
            return st
        for jid, job in st["jobs"].items():
            path = os.path.join(self.outdir, "jobs", "%s.json" % jid)
            if os.path.exists(path):
                with open(path) as f:
                    res = json.load(f)
                job.update(status=res["status"], exit=res["exit"], wall=res["wall"])
            else:
                job.update(status="failed", exit=None)  # job left no result
        if all(job["status"] == DONE for job in st["jobs"].values()):
            self.onDoneAction(st)
            st[self.task] = DONE
        else:
            st[self.task] = REPROCESS
        st["finished"] = _now()
        self.save_status(st)
        return st

    def wait(self, poll=300):
        while True:
            st = self.check()
            if st[self.task] != PROCESSING:
                return st
            if any(job["status"] == "held" for job in st["jobs"].values()):
                print("Stopped waiting: a job is held. Fix it, then run `status` again.")
                return st
            time.sleep(poll)

    def resubmit(self, wait=True, poll=300, force=False):
        st = self.load_status()
        if st[self.task] == DONE and force:
            st[self.task] = REPROCESS  # saved only once the jobs are submitted
        if st[self.task] != REPROCESS:
            print("%s is '%s'; resubmit needs 'reprocess' (or 'done' with --force)"
                  % (self.task, st[self.task]))
            return st
        return self.submit(wait=wait, poll=poll, st=st)
