#!/usr/bin/env python3
"""
Validate the FAST-ZETA reformulation of the fast-zeta kernel against the proven O(d^2) reference.

Claims to validate:
 (Z1) down-zeta on the ideal lattice via restricted Yates: process elements in a linear
      extension order; for each ideal B with x maximal in B: f[B] += f[B \ {x}].
      Result: f_hat[B] = sum_{D in L, D subseteq B} f[D].
 (Z2) up-zeta dual (reverse order, transposed update) gives g_hat[D] = sum_{I superseteq D} g[I].
 (Z3) children collapse: S_k = sum_J sum_{D subseteq B*(J)} (c_sub[J]+c_sup[D])^k
                            = sum_J sum_j C(k,j) c_sub[J]^(k-j) * M_j(B*(J)),
      with M_j = down-zeta of (c_sup ** j), M_0 = c_sub. Must reproduce F_k(Q) exactly.
Total kernel cost: 1 up-zeta + 5 down-zetas + O(d) recombination  ~ O(d * avg#max(B)).
"""
from itertools import combinations, product
from math import comb
import sys
sys.path.insert(0, '.')
from reference_moments import all_labeled_posets, ideals, F_moments

def linext(pred, n):
    placed = 0; order = []
    while len(order) < n:
        for x in range(n):
            if not (placed >> x) & 1 and (pred[x] & ~placed) == 0:
                order.append(x); placed |= 1 << x
    return order

def zeta_down(L, idx, order, f):
    """f_hat[B] = sum over D in L, D subseteq B of f[D], via restricted Yates. O(sum #max)."""
    g = list(f)
    for x in order:                      # linear-extension order
        bit = 1 << x
        for bi, B in enumerate(L):
            if B & bit:
                Bm = B & ~bit
                if Bm in idx and is_max(B, x, pred_g):   # x maximal in B
                    g[bi] += g[idx[Bm]]
    return g

# helper: x maximal in B iff no y in B with x < y  (pred_g[y] contains x)
def is_max(B, x, pred):
    m = B
    while m:
        y = (m & -m).bit_length() - 1
        if y != x and (pred[y] >> x) & 1:
            return False
        m &= m - 1
    return True

def zeta_up(L, idx, order, f, pred):
    """g_hat[D] = sum over I in L, I superseteq D of f[I]. Dual: reverse order, transposed."""
    g = list(f)
    for x in reversed(order):
        bit = 1 << x
        for bi, B in enumerate(L):
            if B & bit and is_max(B, x, pred):
                g[idx[B & ~bit]] += g[bi]
    return g

def F_moments_zeta(pred, n, kmax=4):
    global pred_g
    pred_g = pred
    L = ideals(pred, n); d = len(L)
    idx = {B: i for i, B in enumerate(L)}
    order = linext(pred, n)
    full = (1 << n) - 1
    # c_sup via up-zeta of ones
    c_sup = zeta_up(L, idx, order, [1]*d, pred)
    # M_j via down-zeta of c_sup^j (j=0 gives c_sub)
    M = [zeta_down(L, idx, order, [c_sup[i]**j for i in range(d)]) for j in range(kmax+1)]
    c_sub = M[0]
    S = [0]*(kmax+1)
    for jdx, J in enumerate(L):
        U = full & ~J
        B = J; u = U
        while u:
            x = (u & -u).bit_length() - 1
            B &= pred[x]; u &= u - 1
        bidx = idx[B]                      # B*(J) is an ideal -> in L
        cj = c_sub[jdx]
        for k in range(kmax+1):
            S[k] += sum(comb(k, j) * cj**(k-j) * M[j][bidx] for j in range(k+1))
    return S

if __name__ == "__main__":
    pred_g = None
    for n in (3, 4, 5):
        nfail = 0; cnt = 0
        for pred in all_labeled_posets(n):
            a = F_moments(pred, n, 4)          # proven O(d^2) reference
            b = F_moments_zeta(pred, n, 4)     # fast-zeta reformulation
            cnt += 1
            if a != b:
                nfail += 1
                if nfail <= 3: print(f"  MISMATCH n={n} pred={pred}: ref={a} zeta={b}")
        print(f"n={n}: {cnt} posets, {'ALL MATCH' if nfail==0 else f'{nfail} MISMATCHES'}")
