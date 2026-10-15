# Pinned ROOT environment for the isotrack automator jobs.
# Same LCG view as the EG25C1 reference run (ROOT 6.40.04, gcc13).
export LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_110a/x86_64-el9-gcc13-opt

# The LCG setup script uses unset variables, so relax `set -u` while sourcing it.
case $- in
  *u*) set +u; source "$LCG_VIEW/setup.sh"; set -u ;;
  *)   source "$LCG_VIEW/setup.sh" ;;
esac
