#!/bin/bash
# HTCondor job: build isotrack from the shipped tarball, run mk_IsoTrack.C for
# one era and copy the outputs to EOS.
# Arguments (Automator order): TASK JOB_ID EOS_DIR ERA FILELIST TAG
set -euo pipefail

TASK=$1
JOB_ID=$2
EOS_DIR=$3
ERA=$4
FILELIST=$5
TAG=$6

SCRATCH=$PWD
COMMIT=${ISOTRACK_COMMIT:-unknown}
START=$(date +%s)
RC=1  # stays non-zero unless the run and the copy to EOS both succeed

# Always leave a job result, also when a step fails: the stand-in for
# `ecalautomation.py jobctrl --done/--failed`.
finish() {
  local status=failed
  [ "$RC" -eq 0 ] && status=done
  local json="{\"task\": \"$TASK\", \"job\": \"$JOB_ID\", \"status\": \"$status\", \"exit\": $RC, \"wall\": $(( $(date +%s) - START )), \"commit\": \"$COMMIT\", \"host\": \"$(hostname)\"}"
  echo "$json" > "$SCRATCH/job_${JOB_ID}.json"
  { mkdir -p "$EOS_DIR/jobs" && echo "$json" > "$EOS_DIR/jobs/${JOB_ID}.json"; } \
    || echo "WARNING: could not write $EOS_DIR/jobs/${JOB_ID}.json"
  echo "job result: $json"
}
trap finish EXIT

# 1. Unpack the code; the macro needs these folders before it creates them.
mkdir src
tar -xzf isotrack.tar.gz -C src
cd src
mkdir -p rootfiles pdf logs
source automator/env.sh

echo "task:   $TASK, job $JOB_ID"
echo "host:   $(hostname)"
echo "date:   $(date -u +%FT%TZ)"
echo "commit: $COMMIT"
klist || true
echo "root:   $(root-config --version)"

# 2. Compile and check that every library was built.
root -l -b -q mk_compile.C > logs/compile.log 2>&1
for lib in IsoTrack drawIsoTrack compareCorrFact drawPileup rebinProfiles \
           drawCovariance hybridCorrFact IsoTrackV2; do
  if [ ! -f ${lib}_C.so ]; then
    echo "ERROR: ${lib}_C.so was not built"
    cat logs/compile.log
    exit 2
  fi
done

# 3. Run.
RUN_START=$(date +%s)
set +e
root -l -b -q "mk_IsoTrack.C(\"$ERA\",\"$TAG\",\"$SCRATCH/$FILELIST\")"
ROOT_RC=$?
set -e
WALL=$(( $(date +%s) - RUN_START ))
echo "root exit code: $ROOT_RC, wall time: ${WALL}s"

# 4. Record the analysis switches used.
sed -n '/^\/\/ Override new developments/,/^\/\/ Merge depths/p' IsoTrack.C > settings.txt

# 5. Copy the outputs (this run's pdf/ only: the scratch dir starts empty).
#    After a failed run copy what exists; after a good run every copy must work.
copy() { if [ "$ROOT_RC" -eq 0 ]; then cp -rp "$@"; else cp -rp "$@" 2>/dev/null || true; fi; }
mkdir -p "$EOS_DIR"
copy IsoTrack_${TAG}_${ERA}.root "$EOS_DIR/"
copy rootfiles/*${TAG}_${ERA}* "$EOS_DIR/"
copy pdf "$EOS_DIR/"
copy logs/compile.log settings.txt "$EOS_DIR/"

# 6. README, as for the reference run.
cat > "$EOS_DIR/README.txt" <<README
IsoTrack run: era=$ERA version=$TAG
task:       $TASK, job $JOB_ID
date:       $(date -u +%FT%TZ)
host:       $(hostname)
commit:     $COMMIT
LCG view:   $LCG_VIEW
ROOT:       $(root-config --version)
wall time:  ${WALL}s
root exit:  $ROOT_RC
inputs ($FILELIST):
$(cat "$SCRATCH/$FILELIST")
README

# 7./8. The EXIT trap writes the job result; exit with root's code.
RC=$ROOT_RC
exit $RC
