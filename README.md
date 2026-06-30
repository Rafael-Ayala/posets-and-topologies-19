# Labeled partial orders and topologies on 19 points

Code and data for the computation of

- **P(19) = A001035(19)** = the number of labeled partial orders (equivalently labeled
  T0 topologies) on 19 points, and
- **T(19) = A000798(19)** = the number of labeled topologies on 19 points,

obtained by harvesting the 16-point moment **G(16,3)** from the unlabeled 15-point posets by
one-point insertion, using a single-pass moment-transfer kernel. See the accompanying paper
(in preparation) for the method; this repository contains the executable pipeline, the validation
references, and the per-shard residue table.

## Results

```
P(19) = A001035(19) = 646099441937791106493755218560442089979
T(19) = A000798(19) = 689054943207246404281592791142107048261
G(16,3) = 5322963351172775869497071016032650486
```

`T(19)` is obtained from `P(19)` by the Stirling transform (`assemble_p19.py`); `G(16,3)` is the
new 16-point moment the computation produces.

## Method

For each unlabeled 15-point poset `Q`, the kernel reads the entire 16-point moment
`sum_z d(Q+z)^k` off `Q`'s ideal lattice in one pass (up/down zeta transforms plus a binomial
collapse), accumulating it modulo four 61-bit primes. Summed over all
A000112(15) = 68,275,077,901,156 parents and reconstructed by the Chinese Remainder Theorem,
this gives `G(16,k)` for `k = 0..4`. `P(19)` then follows from the Erné–Stege reduction

```
P(19) = A - 11628*G(15,4) + 969*G(16,3)
```

where `A` is a fixed constant absorbing `P(18)` and the lower-diagonal moments `G(m,19-m)` for
`m <= 14`, and `G(15,4) = 846002793378179474085677125510787278` is a moment over the 15-point posets.
Both are summed directly from the antichain-count histograms in `data/d_histograms.txt` — a complete
enumeration of the posets on at most 15 points — and reconstructed and checked against these published
values by `scripts/assemble_constant.py`. Only `G(16,3)` requires the 16-point harvest above.

This route combines the classical Erné–Stege moment reduction (1991) with the isomorphism-free harvest
of Heitzig and Reinhold (2000), carried one diagonal beyond the published record; the contribution is the
computed values, not a new method.

## Repository layout

```
src/ 
  poset_moment_filter.c   production kernel (fast-zeta, C): reads digraph6 m-point posets from
                          stdin, emits the moment residues sum_z d(Q+z)^k mod four 61-bit primes
  reference_moments.py    O(d^2) ideal-lattice moment reference (the recurrence the C kernel
                          uses), cross-checked against brute-force enumeration of labeled posets
  zeta_kernel_ref.py      checks the kernel's fast-zeta math bit-exactly against reference_moments
scripts/
  run_shards.sh           runs the genposetg | poset_moment_filter pipeline over the shards;
                          idempotent, so it resumes after interruption
  harvest.py              validates the shards, CRT-reconstructs G(16,k), and checks it against
                          the anchors (the known values listed in step 2 below)
  assemble_p19.py         G(16,k) -> P(19) -> congruence/magnitude checks -> A000798(19)
  assemble_constant.py    reconstructs the constant A and G(15,4) from data/d_histograms.txt (the
                          lower-diagonal moments), checked against the published values
  make_shards_csv.py      builds the per-shard residue table (data/shards.csv) from shard outputs
  verify_residue_table.py sums data/shards.csv.gz prime by prime, CRT-reconstructs G(16,k),
                          and checks it against the published values and A000112(15)
data/
  expected_values.json    G(13,k), G(14,k), G(15,k), k=0..4 (validation at lower diagonals)
  moments.csv             G(m,k) = sum over labeled posets on m points of d^k, for m=0..16,
                          k=0..4 (G_m_0 = A001035(m)); the moment table behind the reduction
  d_histograms.txt        antichain-count (d) histograms of the labeled posets on n=1..15 points;
                          assemble_constant.py sums the lower-diagonal moments and A from these
  shards.csv.gz           per-shard residue table (gzip): one row per shard -- its index, then the
                          20 moment residues that shard contributed (plus 3 diagnostic counts)
```

## Dependencies

- [`nauty`/`gtools`](https://pallini.di.uniroma1.it/) — provides `genposetg` (poset generation)
  and the headers (`nauty.h`, `gtools.h`) the kernel links against.
- A C compiler with `__int128`, and Python 3 for the harvest/assembly scripts (standard
  library only).

Build the kernel against your `nauty` tree, e.g.

```
cc -O3 -march=native -I<nauty_dir> src/poset_moment_filter.c <nauty_dir>/nauty.a -o poset_moment_filter
```

## Reproducing the computation

### Verify the published values without rerunning

The shipped residue table reconstructs every result in seconds:

```
python3 scripts/verify_residue_table.py data/shards.csv.gz
```

This sums `data/shards.csv.gz` prime by prime, CRT-reconstructs `G(16,0..4)`, and checks them
against the published values and `parents == A000112(15)`, printing `VERIFY: PASSED`.

The lower half of the reduction needs no rerun either: `python3 scripts/assemble_constant.py`
reconstructs the constant `A` and `G(15,4)` from the antichain-count histograms in
`data/d_histograms.txt`, checks every histogram total against `A001035`, reproduces `P(2..16)` via
the Brinkmann–McKay identity, and confirms `A` and `G(15,4)` against the published values.

The moment arithmetic can also be checked with no kernel build: `python3 src/reference_moments.py`
validates the O(d^2) ideal-lattice algorithm against brute-force enumeration of labeled posets,
and `python3 src/zeta_kernel_ref.py` validates the fast-zeta reformulation against that reference
(both standard-library only).

Steps 1-3 below reproduce the residue table itself from scratch, a large computation.

1. **Generate + filter** (the shard partition). For each shard `i` in `0 .. NPARTS-1`:

   ```
   genposetg 15 t m i NPARTS | ./poset_moment_filter > shard_i.txt
   ```

   The production run used `NPARTS = 307200`. The shards are independent and may be run in any
   order and at any concurrency. `run_shards.sh` wraps the loop (`JOBS` sets the number run
   concurrently); alternatively, the per-shard command can be driven directly under any
   scheduler. It is idempotent — complete shards are skipped, so it resumes after interruption.

2. **Harvest**: validate all shards and reconstruct the moments.

   ```
   python3 scripts/harvest.py <results_dir> --nparts 307200 --n15
   ```

   `<results_dir>` is wherever the shard files were written — `run_shards.sh`'s `OUTDIR`
   (default `./shards`), or, if the per-shard command was run directly, the directory holding
   the `shard_i.txt` files.

   The *anchors* are quantities already known independently that a correct reconstruction must
   reproduce exactly; `harvest.py` recomputes each and aborts (`HARVEST GATE: FAILED`) on any
   mismatch, so a corrupted shard or CRT slip is caught before it can reach `P(19)`. With
   `--n15` the anchors are:

   - `parents == A000112(15)` — every one of the 68,275,077,901,156 unlabeled 15-point parents
     was summed exactly once (no missing or duplicated shards);
   - `sum w == P(15)` modulo each prime — the parent weights `15!/|Aut Q|` sum to the labeled
     15-point total `P(15)`;
   - `G(16,0) == P(16)` — the 0th moment is a count, so it must equal the labeled 16-point
     total `P(16)`;
   - `G(16,1)` and `G(16,2)` — match the closed-form values the Erné–Stege relations predict;
   - `G(16,3)` — the new moment, which has no independent reference value; it is checked
     against the predicted magnitude window `[5.317e36, 5.369e36]`.

3. **Assemble**: turn the reconstructed `G(16,k)` into `P(19)` and run the external checks.

   ```
   python3 scripts/assemble_p19.py <G16_0> <G16_1> <G16_2> <G16_3> [<G16_4>]
   ```

   This reports `P(19)`, the OEIS congruence check, the magnitude, and `A000798(17..19)`.
   The arguments are the `G(16,k)` integers `harvest.py` printed in step 2, in order
   (`G(16,4)` is optional).

The same kernel run one diagonal lower (over the posets on at most 14 points) produces
`G(13,k)`, `G(14,k)`, `G(15,k)`, recorded in `data/expected_values.json`, which independently
reproduce the known `P(17)` and `P(18)` — the validation ladder for the 16-point run.

## Residue table

`data/shards.csv.gz` (gzip-compressed) has one row per shard (24 columns):

- `shard` — the shard index `i`, the **only** input: shard `i` was produced by
  `genposetg 15 t m i 307200 | poset_moment_filter`, so `i` alone fixes which block of the
  15-point posets that row covers.
- `m{k}_p{j}` (the last 20 columns, `k = 0..4`, `j = 0..3`) — the **output**: moment `k`
  reduced modulo the `j`-th of the four 61-bit primes. Summing these prime by prime over all
  307,200 shards and applying CRT reproduces `G(16,k)`.
- `parents`, `children`, `maxd` — diagnostic counts (the number of 15-point parents and
  16-point children the shard processed, and the largest ideal-lattice size `d`, i.e. the
  number of order ideals or downsets, seen among the shard's parents); not part of the
  residue sums, kept only for sanity checks.

Most readers consume the shipped `data/shards.csv.gz` directly. To regenerate it from a fresh
shard rerun, build the table from the shard outputs and gzip it:

```
python3 scripts/make_shards_csv.py <results_dir> shards.csv --nparts 307200 && gzip shards.csv
```

`scripts/verify_residue_table.py` then consumes the resulting `shards.csv.gz`.

## References

- M. Erné and K. Stege, *Counting finite posets and topologies*, Order **8** (1991), 247–265.
- J. Heitzig and J. Reinhold, *The number of unlabeled orders on fourteen elements*, Order **17** (2000), 333–341.
- G. Brinkmann and B. D. McKay, *Posets on up to 16 points*, Order **19** (2002), 147–179.
- A. Björklund, T. Husfeldt, P. Kaski, M. Koivisto, J. Nederlof, P. Parviainen, *Fast zeta transforms for lattices with few irreducibles*, SODA 2012, 1436–1444.
- B. D. McKay and A. Piperno, *Practical graph isomorphism, II*, J. Symbolic Comput. **60** (2014), 94–112 (nauty/gtools).
- OEIS: [A001035](https://oeis.org/A001035), [A000798](https://oeis.org/A000798), [A000112](https://oeis.org/A000112).

## License

The source code (`src/`) and scripts (`scripts/`) are released under the MIT License (see `LICENSE`).
The data files (`data/`) are released under CC-BY-4.0 (see `data/LICENSE.txt`).

## Citation

If you use these values or this code, please cite the accompanying paper (in preparation) and the OEIS
entries [A001035](https://oeis.org/A001035) and [A000798](https://oeis.org/A000798).
