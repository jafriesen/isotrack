#!/usr/bin/env python3
"""Compare a run's correction factors with a reference run, bin by bin.

Checks every hf_dd_* histogram (bin contents and errors, including under- and
overflow) in CorrFact_<tag>_<era>.root and its _even/_odd variants. Prints the
maximum absolute difference per histogram and exits 0 only if all are identical.
Needs pyROOT: run under automator/env.sh.
"""
import argparse
import os
import sys

import ROOT

VARIANTS = ("", "_even", "_odd")


def histograms(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("cannot open %s" % path)
    hists = {}
    for key in f.GetListOfKeys():
        if key.GetName().startswith("hf_dd_"):
            h = key.ReadObj()
            h.SetDirectory(0)
            hists[key.GetName()] = h
    f.Close()
    return hists


def compare(new, ref):
    """Return (max |content diff|, max |error diff|), or None if binnings differ."""
    if new.GetNcells() != ref.GetNcells():
        return None
    dc = de = 0.0
    for i in range(new.GetNcells()):
        dc = max(dc, abs(new.GetBinContent(i) - ref.GetBinContent(i)))
        de = max(de, abs(new.GetBinError(i) - ref.GetBinError(i)))
    return dc, de


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir", required=True, help="run output directory")
    p.add_argument("--tag", required=True, help="version tag in the run's file names")
    p.add_argument("--era", required=True)
    p.add_argument("--ref-dir", default="/eos/user/j/jfriesen/isotrack-dev/reference/EG25C1")
    p.add_argument("--ref-tag", default="lxplus_ref")
    p.add_argument("--ref-era", help="default: same as --era")
    args = p.parse_args()
    ref_era = args.ref_era or args.era

    identical = True
    for v in VARIANTS:
        new_path = os.path.join(args.dir, "CorrFact%s_%s_%s.root" % (v, args.tag, args.era))
        ref_path = os.path.join(args.ref_dir, "CorrFact%s_%s_%s.root" % (v, args.ref_tag, ref_era))
        print("%s\n  vs %s" % (new_path, ref_path))
        new, ref = histograms(new_path), histograms(ref_path)
        for name in sorted(set(new) | set(ref)):
            if name not in new or name not in ref:
                print("  %-10s only in %s" % (name, "run" if name in new else "reference"))
                identical = False
                continue
            res = compare(new[name], ref[name])
            if res is None:
                print("  %-10s different binning" % name)
                identical = False
                continue
            print("  %-10s max |d content| = %.3g, max |d error| = %.3g" % (name, res[0], res[1]))
            identical = identical and res == (0.0, 0.0)
    print("IDENTICAL" if identical else "DIFFERENT")
    return 0 if identical else 1


if __name__ == "__main__":
    sys.exit(main())
