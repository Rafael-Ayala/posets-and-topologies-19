#!/usr/bin/env python3
"""Verify the per-shard residue table: sum the residues prime by prime, CRT-reconstruct
G(16,k), and check against the published values and A000112(15).

Usage: verify_residue_table.py <shards.csv|shards.csv.gz>
"""
import sys, csv, gzip
from functools import reduce

PRIMES = [2305843009213693951, 2305843009213693921, 2305843009213693907, 2305843009213693723]
A000112_15 = 68275077901156
KNOWN = {
    0: 83480529785490157813844256579,            # G(16,0) = P(16)
    1: 28441643117705315333254490986318,          # G(16,1)  (Erne-Stege)
    2: 11344858065618251316427764256980898,       # G(16,2)  (Erne-Stege)
    3: 5322963351172775869497071016032650486,     # G(16,3)  (the new moment)
    4: 2954997625790351969485154266039478036626,  # G(16,4)
}


def crt(rs, ps):
    M = reduce(lambda a, b: a * b, ps)
    x = 0
    for r, p in zip(rs, ps):
        Mi = M // p
        x = (x + r * Mi * pow(Mi % p, p - 2, p)) % M
    return x


def main():
    path = sys.argv[1]
    opener = gzip.open if path.endswith('.gz') else open
    S = {(k, j): 0 for k in range(5) for j in range(4)}
    parents = rows = bad = 0
    with opener(path, 'rt') as f:
        rd = csv.reader(f)
        next(rd)  # header
        for row in rd:
            if len(row) != 24:
                bad += 1
                continue
            rows += 1
            parents += int(row[1])
            for k in range(5):
                for j in range(4):
                    S[(k, j)] = (S[(k, j)] + int(row[4 + k * 4 + j])) % PRIMES[j]
    G = [crt([S[(k, j)] for j in range(4)], PRIMES) for k in range(5)]

    print(f"rows = {rows}   malformed = {bad}")
    ok = (parents == A000112_15) and bad == 0
    print(f"  parents sum == A000112(15): {'OK' if parents == A000112_15 else f'MISMATCH ({parents})'}")
    for k in range(5):
        good = G[k] == KNOWN[k]
        ok &= good
        print(f"  G(16,{k}): {'OK' if good else 'MISMATCH'}   {G[k]}")
    print("VERIFY:", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
