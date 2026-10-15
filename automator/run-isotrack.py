#!/usr/bin/env python3
"""Submit and track IsoTrack calibration jobs on HTCondor (one job per era).

Run on lxplus with the system python3, from anywhere in the isotrack checkout.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from isotrack_handler import IsoTrackCalibHandler  # noqa: E402

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EOS_RUNS = "/eos/user/j/jfriesen/isotrack-dev/runs"
CONDOR_RUNS = os.path.expanduser("~/isotrack-dev/condor/runs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("submit", "status", "resubmit"):
        p = sub.add_parser(name)
        p.add_argument("--era", required=True)
        p.add_argument("--tag", required=True)
        p.add_argument("--outdir", help="default: %s/<era>/<tag>" % EOS_RUNS)
        p.add_argument("--condor-dir", help="default: %s/<era>_<tag>" % CONDOR_RUNS)
        if name != "status":
            p.add_argument("--file-list", required=True,
                           help="list of input ROOT files (transferred with the job)")
            p.add_argument("--dry-run", action="store_true",
                           help="write the tarball and submit file, but do not submit")
            p.add_argument("--no-wait", action="store_true",
                           help="return after submitting; use `status` to follow up")
            p.add_argument("--poll", type=int, default=300, help="seconds between checks")
        if name == "resubmit":
            p.add_argument("--force", action="store_true", help="also resubmit a 'done' task")
    args = parser.parse_args()

    handler = IsoTrackCalibHandler(
        era=args.era, tag=args.tag,
        file_list=getattr(args, "file_list", None),
        outdir=args.outdir or os.path.join(EOS_RUNS, args.era, args.tag),
        condor_dir=args.condor_dir or os.path.join(CONDOR_RUNS, "%s_%s" % (args.era, args.tag)),
        repo_dir=REPO_DIR,
        dry_run=getattr(args, "dry_run", False))

    if args.command == "submit":
        st = handler.submit(wait=not args.no_wait, poll=args.poll)
    elif args.command == "resubmit":
        st = handler.resubmit(wait=not args.no_wait, poll=args.poll, force=args.force)
    else:
        st = handler.check()
    print(json.dumps(st, indent=2, sort_keys=True))
    return 0 if st[handler.task] in ("new", "processing", "done") else 1


if __name__ == "__main__":
    sys.exit(main())
