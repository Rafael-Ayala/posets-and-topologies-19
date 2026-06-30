#!/usr/bin/env python3
"""Harvest poset_moment_filter shard outputs with FULL VALIDATION (audit-hardened).
Usage: harvest.py <results_dir> [--nparts N] [--n15]
Globs shard files internally (no ARG_MAX issues), validates each shard's completeness,
detects duplicates/gaps, enforces every global anchor, then CRT-reconstructs moments."""
import sys, os, re, glob
from functools import reduce

PRIMES = [2305843009213693951, 2305843009213693921, 2305843009213693907, 2305843009213693723]
A000112_15 = 68275077901156
P15 = 77567171020440688353049939
P16 = 83480529785490157813844256579
ES_G16_1 = 28441643117705315333254490986318
ES_G16_2 = 11344858065618251316427764256980898

def crt(residues, primes):
    M = reduce(lambda a, b: a * b, primes)
    x = 0
    for r, p in zip(residues, primes):
        Mi = M // p
        x = (x + r * Mi * pow(Mi % p, p - 2, p)) % M
    return x, M

def main():
    rdir = sys.argv[1]
    nparts = None
    for i, a in enumerate(sys.argv):
        if a == '--nparts': nparts = int(sys.argv[i + 1])
    check_n15 = '--n15' in sys.argv

    files = sorted(glob.glob(os.path.join(rdir, 'shard_*.txt')))
    print(f"found {len(files)} shard files in {rdir}")

    parents = 0; children = 0; maxd = 0
    mom = {(k, p): 0 for k in range(5) for p in PRIMES}
    sumw = {p: 0 for p in PRIMES}
    target_n = None
    seen_ids = set(); bad = []
    for fn in files:
        sid = int(re.search(r'shard_(\d+)\.txt', fn).group(1))
        if sid in seen_ids:
            bad.append((sid, 'DUPLICATE ID')); continue
        seen_ids.add(sid)
        got_moms = 0; got_parents = False
        local = []
        try:
            for line in open(fn):
                m = re.match(r'poset_moment_filter results \(target G\((\d+),k\)\)', line)
                if m:
                    tn = int(m.group(1))
                    if target_n is None: target_n = tn
                    elif tn != target_n: bad.append((sid, f'target mismatch {tn}')); break
                m = re.match(r'parents (\d+)', line)
                if m: parents += int(m.group(1)); got_parents = True
                m = re.match(r'children (\d+)', line)
                if m: children += int(m.group(1))
                m = re.match(r'maxd (\d+)', line)
                if m: maxd = max(maxd, int(m.group(1)))
                m = re.match(r'sumw prime (\d+) res (\d+)', line)
                if m: local.append(('w', int(m.group(1)), int(m.group(2))))
                m = re.match(r'mom (\d+) prime (\d+) res (\d+)', line)
                if m:
                    local.append((int(m.group(1)), int(m.group(2)), int(m.group(3))))
                    got_moms += 1
        except OSError as e:
            bad.append((sid, f'READ ERROR {e}')); continue
        if not got_parents or got_moms != 20:
            bad.append((sid, f'INCOMPLETE (parents={got_parents}, mom_lines={got_moms}/20)'))
            continue
        for (k, p, r) in local:
            if k == 'w': sumw[p] = (sumw[p] + r) % p
            else: mom[(k, p)] = (mom[(k, p)] + r) % p

    print(f"\n=== shard validation ===")
    if bad:
        print(f"  {len(bad)} BAD SHARDS (excluded from sums, MUST be re-run):")
        for sid, why in bad[:20]: print(f"    shard {sid}: {why}")
        if len(bad) > 20: print(f"    ... and {len(bad)-20} more")
    if nparts:
        missing = set(range(nparts)) - seen_ids
        print(f"  missing shard IDs: {len(missing)}" + (f" e.g. {sorted(missing)[:10]}" if missing else ""))
    print(f"  usable shards: {len(seen_ids) - len(bad)}")

    print(f"\n=== global anchors ===")
    ok_all = True
    def chk(name, got, exp):
        nonlocal ok_all
        good = (got == exp); ok_all &= good
        print(f"  {name}: {'OK' if good else f'FAIL (got {got}, expect {exp})'}")
    if check_n15:
        chk("parents == A000112(15)", parents, A000112_15)
        for p in PRIMES: chk(f"sum-w mod {p} == p(15)", sumw[p], P15 % p)
    else:
        print(f"  parents = {parents}  (pass --n15 to enforce = A000112(15))")
        print(f"  sum-w residues: {[sumw[p] for p in PRIMES]}")
    print(f"  children = {children}, maxd = {maxd}")

    print(f"\n=== CRT-reconstructed moments ===")
    vals = []
    for k in range(5):
        val, M = crt([mom[(k, p)] for p in PRIMES], PRIMES)
        vals.append(val)
        print(f"G({target_n},{k}) = {val}")
    if check_n15 and target_n == 16:
        chk("G(16,0) == p(16)", vals[0], P16)
        chk("G(16,1) == ES-derived", vals[1], ES_G16_1)
        chk("G(16,2) == ES-derived", vals[2], ES_G16_2)
        tlo, thi = 5316669 * 10**30, 5368268 * 10**30
        chk2 = tlo <= vals[3] <= thi
        print(f"  G(16,3) in magnitude window [5.317e36, 5.369e36]: {'OK' if chk2 else 'OUTSIDE'}")
        ok_all &= chk2
    print(f"\nHARVEST GATE: {'PASSED' if ok_all and not bad else 'FAILED / ATTENTION'}")
    return 0 if ok_all and not bad else 1

if __name__ == '__main__':
    sys.exit(main())
