#!/usr/bin/env python3
"""
Reconstruct the lower-diagonal inputs to the Erne-Stege reduction

    P(19) = A - 11628*G(15,4) + 969*G(16,3)

directly from the antichain-count histograms of the posets on <= 15 points
(data/d_histograms.txt). Only the 16-point frontier moment G(16,3) needs the big
insertion harvest (harvest.py / verify_residue_table.py); everything below it -- the
constant A and the supporting moment G(15,4) -- is a direct sum over enumerated posets,
which this script carries out and checks against the published values.

For each m, the histogram gives the number of labeled m-point posets with each value of
d(Q) = number of antichains (= number of order ideals). Hence

    G(m,k) = sum_d count[d] * d^k,

so each diagonal moment G(m,19-m) is just a weighted power sum over the enumerated posets.
The constant is

    A = C(20,2)*P(18) - sum_{m=0}^{14} (-1)^(19-m) C(18-m,2) C(19,m) G(m,19-m),

with C(20,2) = 190.

Checks (all must pass):
  - every histogram totals A001035(m)  (a complete enumeration);
  - the histogram moments reproduce the known P(2..16) via the Brinkmann-McKay identity
    P(N) = sum_{k>=1} (-1)^(k-1) C(N,k) G(N-k,k);
  - the reconstructed G(15,4) and A equal the published values, and A - 11628*G(15,4)
    equals B_CONST as used by assemble_p19.py.

Usage: assemble_constant.py [data/d_histograms.txt]
"""
import sys, re, os
from math import comb

A001035 = {0:1,1:1,2:3,3:19,4:219,5:4231,6:130023,7:6129859,8:431723379,9:44511042511,
    10:6611065248783,11:1396281677105899,12:414864951055853499,13:171850728381587059351,
    14:98484324257128207032183,15:77567171020440688353049939,
    16:83480529785490157813844256579,17:122152541250295322862941281269151,
    18:241939392597201176602897820148085023}
A_PUBLISHED    = 5325468436052842213619347019464238237629
G154_PUBLISHED = 846002793378179474085677125510787278
B_CONST        = A_PUBLISHED - 11628 * G154_PUBLISHED   # as used in assemble_p19.py


def load(fn):
    H = {}
    cur = None
    for line in open(fn):
        if line.lstrip().startswith('#'):
            continue
        m = re.match(r'd-histogram for n=(\d+):', line)
        if m:
            cur = int(m.group(1)); H.setdefault(cur, {}); continue
        m = re.match(r'\s*d=(\d+):\s+(\d+)', line)
        if m and cur is not None:
            H[cur][int(m.group(1))] = int(m.group(2))
    return H


def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(os.path.dirname(__file__), '..', 'data', 'd_histograms.txt')
    H = load(fn)
    ok = True

    print("=== histogram completeness (total == A001035(m)) ===")
    for n in sorted(H):
        tot = sum(H[n].values()); good = (tot == A001035.get(n)); ok &= good
        print(f"  m={n:2d}: total {tot}  {'OK' if good else 'FAIL'}")

    def G(m, k):
        return 1 if m == 0 else sum(c * d**k for d, c in H[m].items())

    print("\n=== Brinkmann-McKay identity P(N)=sum (-1)^(k-1) C(N,k) G(N-k,k) ===")
    for N in range(2, 17):
        val = sum((-1)**(k-1) * comb(N, k) * G(N-k, k) for k in range(1, N+1))
        good = (val == A001035[N]); ok &= good
        print(f"  P({N:2d}): {'OK' if good else 'FAIL got '+str(val)}")

    g154 = G(15, 4)
    A = 190 * A001035[18] - sum((-1)**(19-m) * comb(18-m, 2) * comb(19, m) * G(m, 19-m)
                                for m in range(0, 15))
    print("\n=== reconstructed lower-diagonal inputs ===")
    print(f"  G(15,4) = {g154}")
    print(f"            published? {'OK' if g154 == G154_PUBLISHED else 'FAIL'}")
    print(f"  A       = {A}")
    print(f"            published? {'OK' if A == A_PUBLISHED else 'FAIL'}")
    print(f"  A - 11628*G(15,4) == B_CONST (assemble_p19.py)? "
          f"{'OK' if A - 11628*g154 == B_CONST else 'FAIL'}")
    ok &= (g154 == G154_PUBLISHED) and (A == A_PUBLISHED) and (A - 11628*g154 == B_CONST)

    print("\nCONSTANT-ASSEMBLY:", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
