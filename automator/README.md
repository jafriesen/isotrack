# isotrack automator (MVP)

Runs `mk_IsoTrack.C` for one era as a single HTCondor job, with the code taken
from a git commit, and keeps the task state in `<outdir>/status.json`.
A minimal stand-in for the HCAL Automation handlers (`ecalautoctrl`).

Run on lxplus, with the system `python3`, from a clean checkout (commit first:
the job ships `git archive HEAD`).

```bash
# Write the tarball and submit file only
python3 automator/run-isotrack.py submit --era EG25C1 --file-list EG25C1EA.txt --tag mvp-test --dry-run

# Submit and return; then follow up with status
python3 automator/run-isotrack.py submit --era EG25C1 --file-list EG25C1EA.txt --tag mvp-test --no-wait
python3 automator/run-isotrack.py status --era EG25C1 --tag mvp-test

# After a failure (status 'reprocess'), or with --force to redo a 'done' task
python3 automator/run-isotrack.py resubmit --era EG25C1 --file-list EG25C1EA.txt --tag mvp-test
```

Without `--no-wait`, `submit` polls every `--poll` seconds until the job ends.

- Outputs: `--outdir`, default `/eos/user/j/jfriesen/isotrack-dev/runs/<era>/<tag>/`,
  with `status.json`, `jobs/<id>.json`, the ROOT files, `pdf/`, `settings.txt`,
  `compile.log`, `README.txt`, and `logs/` (condor logs, copied when done).
- Condor files: `~/isotrack-dev/condor/runs/<era>_<tag>/` on AFS (condor cannot
  write its logs to EOS).
- Status words: `new`, `processing`, `done`, `reprocess`.

Check a run against the EG25C1 reference (must be identical):

```bash
source automator/env.sh
python3 automator/compare_to_reference.py --dir /eos/user/j/jfriesen/isotrack-dev/runs/EG25C1/mvp-test --tag mvp-test --era EG25C1
```
