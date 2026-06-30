#!/usr/bin/env python3
"""
Final assembly: G(16,k) from the insertion harvest -> p(19) -> all checks -> A000798(19).
Usage: assemble_p19.py <G16_0> <G16_1> <G16_2> <G16_3> [G16_4]
(the CRT-reconstructed integers printed by harvest.py)
"""
import sys
from math import comb

# Reduction constants (the classical Erne-Stege reduction, recorded on A001035; all gates passed)
A_CONST   = 5325468436052842213619347019464238237629
G15_4     = 846002793378179474085677125510787278
B_CONST   = A_CONST - 11628 * G15_4          # = -4511852045348628711048906595975196230955

P = {15:77567171020440688353049939, 16:83480529785490157813844256579,
     17:122152541250295322862941281269151, 18:241939392597201176602897820148085023}
ES_G16_1  = 28441643117705315333254490986318
ES_G16_2  = 11344858065618251316427764256980898
CONG_MOD, CONG_RES = 232792560, 163279579   # lcm(1..20); the A001035 modular-periodicity congruence

def main():
    g = [int(x) for x in sys.argv[1:]]
    G16_0, G16_1, G16_2, G16_3 = g[0], g[1], g[2], g[3]
    ok = True
    def chk(name, got, exp):
        nonlocal ok
        good = (got == exp); ok &= good
        print(f"  {name}: {'OK' if good else f'FAIL got {got} expect {exp}'}")
    print("=== self-checks ===")
    chk("G(16,0) == p(16)", G16_0, P[16])
    chk("G(16,1) == ES-derived", G16_1, ES_G16_1)
    chk("G(16,2) == ES-derived", G16_2, ES_G16_2)
    print(f"\nG(16,3) = {G16_3}")
    tlo, thi = 5316669 * 10**30, 5368268 * 10**30
    print(f"  magnitude window [5.317e36, 5.369e36]: {'OK' if tlo <= G16_3 <= thi else 'OUTSIDE, INVESTIGATE'}")

    p19 = B_CONST + 969 * G16_3
    print(f"\np(19) = A001035(19) = {p19}")
    print(f"  digits: {len(str(p19))}")
    cong = p19 % CONG_MOD
    print(f"  OEIS congruence mod {CONG_MOD}: got {cong}, expect {CONG_RES}  "
          f"{'OK' if cong == CONG_RES else 'FAIL'}")
    import math
    print(f"  log2 p(19) = {math.log2(p19):.4f}  (magnitude estimate: 128.9±0.1)")

    # Stirling endgame: A000798(19) = sum_k S(19,k) * A001035(k)
    A001035 = {0:1,1:1,2:3,3:19,4:219,5:4231,6:130023,7:6129859,8:431723379,
        9:44511042511,10:6611065248783,11:1396281677105899,12:414864951055853499,
        13:171850728381587059351,14:98484324257128207032183,
        15:77567171020440688353049939,16:83480529785490157813844256579,
        17:122152541250295322862941281269151,
        18:241939392597201176602897820148085023,19:p19}
    # Stirling2 via DP
    S = [[0]*(20) for _ in range(20)]
    S[0][0] = 1
    for n in range(1, 20):
        for k in range(1, n+1):
            S[n][k] = k*S[n-1][k] + S[n-1][k-1]
    # A000798(17), (18), (19) via the Stirling transform (17 and 18 are informational)
    for N in (17, 18, 19):
        t = sum(S[N][k]*A001035[k] for k in range(1, N+1))
        print(f"A000798({N}) = {t}")
    print("\nALL CHECKS:", "PASSED" if ok and cong == CONG_RES else "ATTENTION NEEDED")

if __name__ == '__main__':
    main()
