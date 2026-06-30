/*
 * poset_moment_filter.c: fast-zeta filter for G(n+1,k), k=0..4.
 *
 * Reads unlabeled m-point posets (m = n <= 15) from stdin in digraph6 form (transitively-closed
 * digraphs, the output of genposetg ... t) and emits the residues of the moments sum_z d(Q+z)^k
 * modulo four 61-bit primes.
 *
 *  - The per-parent moment is read off the ideal lattice by zeta transforms
 *    (restricted Yates along a linear extension; BHKKNP O(d*n)):
 *      c_sup = up-zeta of ones; M_j = down-zeta of c_sup^j (M_0 = c_sub);
 *      children collapse binomially: sum_{D<=B*(J)} (c_sub[J]+c_sup[D])^k
 *        = sum_j C(k,j) c_sub[J]^{k-j} M_j(B*(J)),  B*(J) = J & G[J],
 *      G[J] = G[J u {x}] & pred[x] by one stored covering arc per ideal.
 *    Validated bit-exactly against the O(d^2) reference (zeta_kernel_ref.py, all posets n<=5,
 *    audit prototype n<=16 incl. antichains).
 *  - Per-parent self-asserts: M0[full] == d and c_sup[empty] == d (transform integrity),
 *    n! % |Aut| == 0, d <= 32768.
 *  - Parent-level anchors: panchor k = sum w * d(Q)^k mod each prime (k=0..4). Summed over a
 *    complete run these equal G(15,k), five 244-bit pipeline checks.
 *  - Deterministic sample emission (1 per 2^20 parents) for offline independent recompute.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAXN 32
#define WORDSIZE 64
#include "nauty.h"
#include "nautinv.h"
#include "gtools.h"

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;
typedef uint16_t u16;
typedef uint8_t u8;

#define MAX_D 32769          /* ideals at n<=15: <= 2^15 */
#define MAX_E (MAX_D * 15)

static u16 L[MAX_D], mmask[MAX_D], Gm[MAX_D];
static u32 IDX[1 << 15];
static u32 asrc[MAX_E], adst[MAX_E];
static u32 aoff[17], acnt[16];
static u8  up_elem[MAX_D];
static u32 up_src[MAX_D];
static u32 csup[MAX_D];
static u64 M0[MAX_D], M1[MAX_D], M2[MAX_D], M3[MAX_D];
static u128 M4[MAX_D];

static const u64 PRIMES[4] = {
    2305843009213693951ULL, 2305843009213693921ULL,
    2305843009213693907ULL, 2305843009213693723ULL };

static u64 acc[5][4], accw[4], panchor[5][4];
static u128 total_children = 0;
static long long total_parents = 0;
static int max_d_seen = 0;

static u64 mulmod(u64 a, u64 b, u64 p) { return (u64)(((u128)a * b) % p); }
static u64 u128mod(u128 x, u64 p) { return (u64)(x % p); }
static void print_u128(FILE *f, u128 x) {
    char buf[44]; int i = 43; buf[i] = 0;
    if (x == 0) { fputc('0', f); return; }
    while (x > 0) { buf[--i] = '0' + (int)(x % 10); x /= 10; }
    fputs(buf + i, f);
}

/* build ideal list + max-masks along a linear extension; fills IDX; returns d */
static int build_ideals(int n, const u32 *pred, int *order) {
    u32 placed = 0; int no = 0;
    while (no < n) {
        for (int x = 0; x < n; x++)
            if (!((placed >> x) & 1) && (pred[x] & ~placed) == 0) {
                order[no++] = x; placed |= (1u << x);
            }
    }
    L[0] = 0; mmask[0] = 0; IDX[0] = 0; int d = 1;
    for (int oi = 0; oi < n; oi++) {
        int x = order[oi]; u16 bit = (u16)(1u << x); u16 px = (u16)pred[x];
        int dprev = d;
        for (int i = 0; i < dprev; i++)
            if ((L[i] & px) == px) {
                if (d >= MAX_D) { fprintf(stderr, "FATAL ideal overflow\n"); exit(1); }
                L[d] = L[i] | bit;
                mmask[d] = (u16)((mmask[i] & ~px) | bit);
                IDX[L[d]] = (u32)d;
                d++;
            }
    }
    return d;
}

int main(void) {
    graph g[1024];
    int lab[MAXN], ptn[MAXN], orbits[MAXN];
    statsblk stats;
    static DEFAULTOPTIONS_DIGRAPH(options);
    options.getcanon = FALSE;

    int n = 0, m_nauty;
    int order[MAXN], pos[MAXN];
    memset(acc, 0, sizeof(acc)); memset(accw, 0, sizeof(accw));
    memset(panchor, 0, sizeof(panchor));

    while (1) {
        boolean digflag;
        graph *gp = readgg(stdin, g, 0, &m_nauty, &n, &digflag);
        if (gp == NULL) break;
        if (n > 15) { fprintf(stderr, "FATAL n=%d > 15\n", n); exit(1); }
        total_parents++;

        /* adjacency -> closure -> pred masks */
        u32 adj[MAXN]; memset(adj, 0, n * sizeof(u32));
        for (int i = 0; i < n; i++) {
            set *gi = GRAPHROW(g, i, m_nauty);
            for (int j = 0; j < n; j++)
                if (ISELEMENT(gi, j)) adj[i] |= (1u << j);
        }
        u32 reach[MAXN]; memcpy(reach, adj, n * sizeof(u32));
        for (int k = 0; k < n; k++)
            for (int i = 0; i < n; i++)
                if (reach[i] & (1u << k)) reach[i] |= reach[k];
        u32 pred[MAXN];
        for (int i = 0; i < n; i++) pred[i] = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (reach[i] & (1u << j)) pred[j] |= (1u << i);

        int d = build_ideals(n, pred, order);
        if (d > max_d_seen) max_d_seen = d;
        for (int oi = 0; oi < n; oi++) pos[order[oi]] = oi;

        /* ---- arc list bucketed by stage (linext position of the removed max element) ---- */
        memset(acnt, 0, n * sizeof(u32));
        for (int b = 0; b < d; b++) {
            u16 t = mmask[b];
            while (t) { int x = __builtin_ctz(t); acnt[pos[x]]++; t &= (u16)(t - 1); }
        }
        aoff[0] = 0;
        for (int s = 0; s < n; s++) aoff[s + 1] = aoff[s] + acnt[s];
        u32 fill[16]; memcpy(fill, aoff, n * sizeof(u32));
        for (int b = 0; b < d; b++) {
            u16 t = mmask[b];
            while (t) {
                int x = __builtin_ctz(t); t &= (u16)(t - 1);
                u32 dst = IDX[(u16)(L[b] & ~(1u << x))];
                u32 a = fill[pos[x]]++;
                asrc[a] = (u32)b; adst[a] = dst;
                up_elem[dst] = (u8)x; up_src[dst] = (u32)b;
            }
        }
        u32 E = aoff[n];

        /* ---- up-zeta: c_sup[D] = #{I in L : I >= D} ---- */
        for (int i = 0; i < d; i++) csup[i] = 1;
        for (int s = n - 1; s >= 0; s--)
            for (u32 a = aoff[s]; a < aoff[s + 1]; a++)
                csup[adst[a]] += csup[asrc[a]];
        if (csup[0] != (u32)d) { fprintf(stderr, "FATAL csup self-check parent %lld\n", total_parents); exit(1); }

        /* ---- fused down-zeta: M_j = zeta(c_sup^j), M_0 = c_sub ---- */
        for (int i = 0; i < d; i++) {
            u64 c = csup[i], c2 = c * c;
            M0[i] = 1; M1[i] = c; M2[i] = c2; M3[i] = c2 * c; M4[i] = (u128)c2 * c2;
        }
        for (int s = 0; s < n; s++)
            for (u32 a = aoff[s]; a < aoff[s + 1]; a++) {
                u32 src = asrc[a], dst = adst[a];
                M0[src] += M0[dst]; M1[src] += M1[dst]; M2[src] += M2[dst];
                M3[src] += M3[dst]; M4[src] += M4[dst];
            }
        u16 full = (u16)((1u << n) - 1);
        if (M0[IDX[full]] != (u64)d) { fprintf(stderr, "FATAL csub self-check parent %lld\n", total_parents); exit(1); }

        /* ---- G masks via stored covering arcs, reverse creation order ---- */
        for (int i = d - 1; i >= 0; i--)
            Gm[i] = (L[i] == full) ? full : (u16)(Gm[up_src[i]] & pred[up_elem[i]]);

        /* ---- combine: binomial collapse ---- */
        u64 s0 = 0; u128 s1 = 0, s2 = 0, s3 = 0, s4 = 0;
        for (int j = 0; j < d; j++) {
            u32 bi = IDX[(u16)(L[j] & Gm[j])];
            u64 cj = M0[j], cj2 = cj * cj;
            u64 m0 = M0[bi], m1 = M1[bi], m2 = M2[bi], m3 = M3[bi];
            s0 += m0;
            s1 += (u128)cj * m0 + m1;
            s2 += (u128)cj2 * m0 + 2 * (u128)cj * m1 + m2;
            s3 += (u128)(cj2 * cj) * m0 + 3 * (u128)cj2 * m1 + 3 * (u128)cj * m2 + m3;
            s4 += (u128)(cj2 * cj2) * m0 + 4 * (u128)(cj2 * cj) * m1
                + 6 * (u128)cj2 * m2 + 4 * (u128)cj * m3 + M4[bi];
        }
        total_children += s0;

        /* ---- |Aut| ---- */
        densenauty(g, lab, ptn, orbits, &options, &stats, m_nauty, n, NULL);
        double aut_d = stats.grpsize1 * pow(10.0, stats.grpsize2);
        u64 aut = (u64)(aut_d + 0.5);
        u64 nfact = 1; for (int i = 2; i <= n; i++) nfact *= (u64)i;
        if (nfact % aut != 0) { fprintf(stderr, "FATAL |Aut|=%llu !| %d! parent %lld\n",
                                        (unsigned long long)aut, n, total_parents); exit(1); }
        u64 w = nfact / aut;

        /* ---- fold ---- */
        for (int p = 0; p < 4; p++) {
            u64 P = PRIMES[p], wm = w % P;
            accw[p] = (accw[p] + wm) % P;
            acc[0][p] = (acc[0][p] + mulmod(wm, (u64)(s0 % P), P)) % P;
            acc[1][p] = (acc[1][p] + mulmod(wm, u128mod(s1, P), P)) % P;
            acc[2][p] = (acc[2][p] + mulmod(wm, u128mod(s2, P), P)) % P;
            acc[3][p] = (acc[3][p] + mulmod(wm, u128mod(s3, P), P)) % P;
            acc[4][p] = (acc[4][p] + mulmod(wm, u128mod(s4, P), P)) % P;
            u64 dm = 1;
            for (int k = 0; k <= 4; k++) {
                panchor[k][p] = (panchor[k][p] + mulmod(wm, dm, P)) % P;
                dm = mulmod(dm, (u64)d % P, P);
            }
        }
        (void)E;

        if ((total_parents & 1048575) == 1) {
            fprintf(stdout, "sample %lld d %d s0 %llu s1m %llu s3m %llu aut %llu\n",
                    total_parents, d, (unsigned long long)s0,
                    (unsigned long long)u128mod(s1, PRIMES[0]),
                    (unsigned long long)u128mod(s3, PRIMES[0]),
                    (unsigned long long)aut);
        }
        if (total_parents % 2000000 == 0)
            fprintf(stderr, "\r  %lldM parents, max_d=%d  ", total_parents / 1000000, max_d_seen);
    }

    fprintf(stderr, "\r  %lld parents done, max_d=%d\n", total_parents, max_d_seen);

    printf("poset_moment_filter results (target G(%d,k))\n", n + 1);
    printf("parents %lld\n", total_parents);
    printf("children "); print_u128(stdout, total_children); printf("\n");
    printf("maxd %d\n", max_d_seen);
    for (int p = 0; p < 4; p++)
        printf("sumw prime %llu res %llu\n",
               (unsigned long long)PRIMES[p], (unsigned long long)accw[p]);
    for (int k = 0; k <= 4; k++)
        for (int p = 0; p < 4; p++)
            printf("mom %d prime %llu res %llu\n", k,
                   (unsigned long long)PRIMES[p], (unsigned long long)acc[k][p]);
    for (int k = 0; k <= 4; k++)
        for (int p = 0; p < 4; p++)
            printf("panchor %d prime %llu res %llu\n", k,
                   (unsigned long long)PRIMES[p], (unsigned long long)panchor[k][p]);
    return 0;
}
