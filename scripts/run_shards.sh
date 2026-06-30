#!/usr/bin/env bash
# Run the shard partition for G(N+1,k). For each shard i in 0 .. NPARTS-1, generate that
# block of the unlabeled N-point posets and filter it into moment residues:
#
#     genposetg N t m i NPARTS | poset_moment_filter > OUTDIR/shard_i.txt
#
# Re-running is idempotent (complete shards are skipped), so it resumes after interruption.
# The full N=15 run used NPARTS=307200. The shards are independent, so run them in
# whatever order and concurrency you like: set JOBS to how many to run at once here, or drive
# the per-shard command with whatever scheduler you have.
#
# Environment (all optional):
#   N         points per parent           (default 15)
#   NPARTS    number of shards             (default 307200)
#   OUTDIR    output directory             (default ./shards)
#   JOBS      parallel workers             (default 4)
#   GENPOSETG path to genposetg            (default genposetg, from nauty/gtools)
#   FILTER    path to the built kernel      (default ./poset_moment_filter)
set -euo pipefail

N=${N:-15}
NPARTS=${NPARTS:-307200}
OUTDIR=${OUTDIR:-shards}
JOBS=${JOBS:-4}
GENPOSETG=${GENPOSETG:-genposetg}
FILTER=${FILTER:-./poset_moment_filter}

mkdir -p "$OUTDIR"
export N NPARTS OUTDIR GENPOSETG FILTER

seq 0 $((NPARTS - 1)) | xargs -P "$JOBS" -I{} bash -c '
    i="{}"
    out="$OUTDIR/shard_${i}.txt"
    if [ -s "$out" ] && grep -q "^mom 4" "$out"; then exit 0; fi   # already complete
    "$GENPOSETG" "$N" t m "$i" "$NPARTS" | "$FILTER" > "$out.tmp" && mv "$out.tmp" "$out"
'
