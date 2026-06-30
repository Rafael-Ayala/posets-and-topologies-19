#!/usr/bin/env python3
"""Build a per-shard residue table (shards.csv) from the poset_moment_filter shard outputs.

Each shard i was produced by
    genposetg N t m i NPARTS | poset_moment_filter > shard_i.txt
This reads every shard_*.txt in a results directory and writes one CSV row per shard:
    shard, parents, children, maxd, m0_p0, m0_p1, ..., m4_p3
i.e. the input (shard index) and the output (the 20 moment residues, k = 0..4 over the
four 61-bit primes). Lines other than the counts and residues are skipped, so the residue table
contains only inputs and outputs.

Usage: make_shards_csv.py <results_dir> <out.csv> [--nparts N]
"""
import sys, os, re, glob

PRIMES = [2305843009213693951, 2305843009213693921, 2305843009213693907, 2305843009213693723]
PIDX = {p: i for i, p in enumerate(PRIMES)}


def main():
    rdir, out = sys.argv[1], sys.argv[2]
    nparts = int(sys.argv[sys.argv.index('--nparts') + 1]) if '--nparts' in sys.argv else None

    files = sorted(glob.glob(os.path.join(rdir, 'shard_*.txt')))
    cols = ['shard', 'parents', 'children', 'maxd'] + [f'm{k}_p{j}' for k in range(5) for j in range(4)]
    seen, nbad = set(), 0
    with open(out, 'w') as fo:
        fo.write(','.join(cols) + '\n')
        for f in files:
            m = re.search(r'shard_(\d+)\.txt$', f)
            if not m:
                continue
            sid = int(m.group(1))
            parents = children = maxd = None
            res = {}
            with open(f) as fi:
                for line in fi:
                    mo = re.match(r'parents (\d+)', line)
                    if mo: parents = int(mo.group(1)); continue
                    mo = re.match(r'children (\d+)', line)
                    if mo: children = int(mo.group(1)); continue
                    mo = re.match(r'maxd (\d+)', line)
                    if mo: maxd = int(mo.group(1)); continue
                    mo = re.match(r'mom (\d+) prime (\d+) res (\d+)', line)
                    if mo and int(mo.group(2)) in PIDX:
                        res[(int(mo.group(1)), PIDX[int(mo.group(2))])] = int(mo.group(3))
            if parents is None or len(res) != 20:
                nbad += 1
                continue
            seen.add(sid)
            row = [sid, parents, children, maxd] + [res[(k, j)] for k in range(5) for j in range(4)]
            fo.write(','.join(str(x) for x in row) + '\n')
    print(f"wrote {len(seen)} shard rows to {out} ({nbad} incomplete/skipped)")
    if nparts:
        missing = sorted(set(range(nparts)) - seen)
        print(f"missing {len(missing)} of {nparts} shards" + (f" e.g. {missing[:10]}" if missing else ""))


if __name__ == '__main__':
    main()
