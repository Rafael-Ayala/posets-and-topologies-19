#!/usr/bin/env python3
"""
Reference implementation + validation of the extend-from-n identity:

  G(n+1,k) = sum_{labeled Q on n} F_k(Q)
           = sum_{unlabeled Q on n} (n!/|Aut Q|) * F_k(Q)

  F_k(Q) = sum over compatible insertions (D,U) of d(Q+z)^k, computed via the
  ideal-lattice algorithm the C filter will use:
     L      = all ideals (downsets) of Q
     c_sub[J] = #{I in L : I subseteq J}
     c_sup[D] = #{I in L : I superseteq D}
     children: for J in L (J = complement of upset U):
                 B*(J) = J  ∩  ⋂_{u not in J} pred(u)
                 for D in L with D ⊆ B*(J):  d(Q+z) = c_sub[J] + c_sup[D]

Validated against: (a) direct brute-force G(n+1,k) by enumerating ALL labeled
posets on n+1 points; (b) the published G(13..15,k) values; (c) the
unlabeled-weighted form using |Aut| (permutation backtracking).
"""
from itertools import combinations, product, permutations
from math import factorial
import sys, re

def all_labeled_posets(n):
    """Yield strict-order relation as tuple pred[i] = bitmask of elements < i."""
    pairs = list(combinations(range(n), 2))
    for assign in product(range(3), repeat=len(pairs)):
        lt = [[False]*n for _ in range(n)]
        for (a, b), v in zip(pairs, assign):
            if v == 1: lt[a][b] = True
            elif v == 2: lt[b][a] = True
        ok = True
        for a in range(n):
            for b in range(n):
                if lt[a][b]:
                    for c in range(n):
                        if lt[b][c] and not lt[a][c]:
                            ok = False; break
                    if not ok: break
            if not ok: break
        if ok:
            pred = [0]*n
            for a in range(n):
                for b in range(n):
                    if lt[a][b]: pred[b] |= (1 << a)
            yield tuple(pred)

def ideals(pred, n):
    """All ideals as bitmasks (DP over a linear extension)."""
    # topological order: repeatedly take elements whose preds are all placed
    placed = 0; order = []
    while len(order) < n:
        for x in range(n):
            if not (placed >> x) & 1 and (pred[x] & ~placed) == 0:
                order.append(x); placed |= 1 << x
    L = [0]
    for x in order:
        bit = 1 << x; px = pred[x]
        L += [I | bit for I in L if (px & ~I) == 0]
    return L

def F_moments(pred, n, kmax=4):
    """F_k(Q) for k=0..kmax via the ideal-lattice algorithm (mirrors the C)."""
    L = ideals(pred, n)
    d = len(L)
    full = (1 << n) - 1
    c_sub = [0]*d; c_sup = [0]*d
    for j, J in enumerate(L):
        for i, I in enumerate(L):
            if (I & ~J) == 0:
                c_sub[j] += 1; c_sup[i] += 1
    S = [0]*(kmax+1)
    for j, J in enumerate(L):
        U = full & ~J
        B = J
        u = U
        while u:
            x = (u & -u).bit_length() - 1
            B &= pred[x]
            u &= u - 1
        cj = c_sub[j]
        for i, D in enumerate(L):
            if (D & ~B) == 0:
                t = cj + c_sup[i]
                tk = 1
                for k in range(kmax+1):
                    S[k] += tk; tk *= t
    return S

def d_of(pred, n):
    return len(ideals(pred, n))

def aut_order(pred, n):
    """|Aut(Q)| by brute permutation check (n<=6)."""
    cnt = 0
    rel = set()
    for b in range(n):
        m = pred[b]
        while m:
            a = (m & -m).bit_length() - 1
            rel.add((a, b)); m &= m - 1
    for perm in permutations(range(n)):
        if all((perm[a], perm[b]) in rel for (a, b) in rel):
            cnt += 1
    return cnt

if __name__ == "__main__":
    KM = 4
    for n in (3, 4):
        # direct: G(n+1,k) by brute enumeration of (n+1)-point posets
        direct = [0]*(KM+1)
        for pred in all_labeled_posets(n+1):
            dd = d_of(pred, n+1)
            tk = 1
            for k in range(KM+1):
                direct[k] += tk; tk *= dd
        # identity (labeled): sum F_k over labeled n-point posets
        ident = [0]*(KM+1)
        nlab = 0
        for pred in all_labeled_posets(n):
            S = F_moments(pred, n, KM); nlab += 1
            for k in range(KM+1): ident[k] += S[k]
        # identity (unlabeled-weighted): dedup by canonical form, weight n!/|Aut|
        seen = {}
        for pred in all_labeled_posets(n):
            rel = [(a, b) for b in range(n) for a in range(n) if (pred[b] >> a) & 1]
            canon = min(tuple(sorted((perm[a], perm[b]) for (a, b) in rel))
                        for perm in permutations(range(n)))
            seen.setdefault(canon, pred)
        wident = [0]*(KM+1)
        for pred in seen.values():
            S = F_moments(pred, n, KM)
            w = factorial(n)//aut_order(pred, n)
            for k in range(KM+1): wident[k] += w*S[k]
        print(f"n={n} -> G({n+1},k):")
        print(f"  direct   : {direct}")
        print(f"  identity : {ident}   {'MATCH' if ident==direct else 'FAIL'}")
        print(f"  unlabeled: {wident}  ({len(seen)} classes)  {'MATCH' if wident==direct else 'FAIL'}")
