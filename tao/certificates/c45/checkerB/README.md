# checker B — independent exact verifier for the C45 upper bound

Written from the mathematical statement only (no code or data of checker A, Griego or
Yoo was read). Two programs:

| file | what it does |
|---|---|
| `histB.cpp` | C++17 (pthreads, no other deps). Exact base histogram `nu delta` for a set of base primes. |
| `evalB.py` | Python 3 stdlib only. Reads a JSON spec, verifies the structure, computes S exactly. |
| `run3.sh` | Runs `histB` on validation sets (a), (b), (c) and prints sha256 / time / memory. |
| `specs/*.json` | Specs for targets (a) and (b). |
| `out/` | Histograms and logs produced on the Mac. |

## Build and run

```
g++ -O3 -march=native -std=c++17 -pthread histB.cpp -o histB      # (clang++ on the Mac)
./histB --selftest                                 # brute-force comparisons, 324 runs
./histB 3 5 7 11 13 17 19 29 31 37 41 61 73 > griego.txt     # progress/checks on stderr
python3 evalB.py specs/griego.json
```

`histB` options: `--threads=N` (default: all cores), `--order=auto|given|q,q,...`,
`--tail=K` (default 3), `--canon=tree|affine|none` (default tree), `--maxstates=N`
(exit 3 if a stored layer gets larger; useful for probing an order), `--plan` (print
the order and cut moduli, then exit), `--brute` (also run the direct enumeration over
a mod M0 and compare), `--bruteonly`, `--quiet`. Primes may be given in any order and
separated by spaces or commas. The output has only the nonzero `delta` entries, in
ascending `nu`. Exit status is nonzero if any internal check fails.

## Algorithm of histB

**Reformulation.** For each base prime q (order d = ord_q(2)), a mod q is either not a
power of 2 (q − d residues, which remove nothing) or 2^s for exactly one s mod d (this
removes the class r ≡ s mod d). By CRT, delta is the weighted distribution of
|Z_T0 minus the union of the chosen classes| over independent choices, one per prime.

**Cut DP.** The primes are processed in a fixed order. After k primes, let
Lp = lcm(processed orders), Lr = lcm(remaining orders) and g_k = gcd(Lp, Lr). The state
is the vector c on Z_{g_k}, where c[x] = #{r in Z_T0 : r ≡ x mod g_k, r survives the
processed primes}. This state is complete. Take r uniform in Z_T0 and fix r mod g_k.
Then (r mod Lp, r mod Lr) runs uniformly over the compatible pairs, so survival of the
past and survival of the future are independent given r mod g_k.

Transition with a prime of order d (d | Lr): let G = lcm(g, d), expand
e[X] = c[X mod g]/(G/g) on Z_G (the division is exact and is checked), zero the class
X ≡ s mod d, and project to g_{k+1} (which divides G) by summing fibres. The weight is
multiplied by q − d for "no class" and by 1 for each class. At the end g = 1 and
c[0] = nu. State entries are at most T0/g, so they are stored as 1 or 2 bytes. Weights
are fixed-width 128/192/256-bit integers, chosen from M0, with overflow checks on
every operation.

**Symmetry.** Let D* = {gcd(d, g) : d remaining}. Any permutation of Z_g that preserves
every partition "x mod e" for e in D* lifts to a permutation of Z_Lr that preserves
the classes of every remaining prime. The lift keeps the higher CRT digits fixed (the
proof is in the header comment of `histB.cpp`). So it leaves the future
nu-distribution unchanged, and states are merged up to this group:
* If D* is a divisibility chain, the group is the iterated wreath product, i.e. the
  automorphisms of the tree of nested blocks. The canonical form is a recursive sort of
  sibling subtrees, and needs no enumeration. Example: when only primes of order 36
  remain, the canonical form is just the sorted vector.
* Otherwise a componentwise subgroup is used. Z_g = prod Z_{l^a}, and each component
  gets a tree group whose levels are the l-adic valuations that occur. The largest
  component is canonicalised by sorting; the others are enumerated (switching to their
  affine subgroups if the enumeration exceeds 3000 elements). If the group equals
  Aff(Z_g), a lexicographic-minimum search over x → u x + t is used, because it is
  faster.

Only orbit-invariance matters for correctness, and every mode is checked against the
brute force.

**Order (auto).** First the primes whose order has 2-adic valuation ≥ 3 (highest
first). Then, for each prime l ≥ 5 dividing T0 (larger l first), the primes whose
order is divisible by l, as one contiguous block. Then the remaining orders not
divisible by 9. The orders divisible by 9 come last. This keeps g_k ≤ 84 and makes the
long middle stretch g = 12. For the 20-prime set the order is
97,257,241,17,127,43,29,113,31,11,41,61,13,5,7,3,73,19,37,109 with
g = 1,48,24,24,12,84,84,84,12,60,60,60,12,12,12,12,12,36,36,36,1.

**Tail.** The last K primes (default 3) are not merged. For every stored state they
are enumerated depth-first with small-integer multipliers, collected in a local
`uint64` histogram (the per-state mass is the product of the tail primes, checked to
be < 2^63), and then flushed as W × count. The last two primes are fused in closed
form: on Z_Q with Q = lcm(dA, dB),
nu = tot − mA[sA] − mB[sB] + [sA ≡ sB mod gcd] · e[CRT(sA, sB)].

**Parallelism and memory.** Each layer is a hash set split into 4096 shards. Each shard
has its own spin lock, an open-addressing index (32-bit tag + 32-bit record index) and
a block arena of records [weight | packed key]. Worker threads take parent blocks from
an atomic counter, and free each parent block once it has been expanded. Cost is
roughly (key bytes + weight bytes + ~16) bytes per state, e.g. about 70 B per state at
g = 36 with 128-bit weights.

**Built-in checks, run on every invocation.**
* sum delta = M0
* sum delta·nu = T0·φ(M0)
* sum delta·nu² = T0 · Σ_{x∈Z_T0} Π_q (q − 2 + [d_q | x])
* sum delta·nu³ = T0 · Σ_{x,y} Π_q (q − #{0, x, y mod d_q})

These are independent closed forms obtained by counting tuples of r's.

**Self-test.** `--selftest` compares 18 small sets, e.g. {3,5,7,13}, {3,5,7,13,17},
{3,7,5,13,17,11,43} and sets built from the "hard" primes 73, 19, 37, 109, 257, 97,
241, 127, 43, 29, 113. Each set runs in 18 configurations (order auto/given × tail
1, 2, 3 × canon tree/affine/none) against a direct enumeration over all a mod M0. That
enumeration uses the literal definition with no use of the order structure: for each a
it takes the union over q of {r in Z_T0 : 2^r ≡ a mod q}. Result: 324/324 match.

## Algorithm of evalB.py

The script checks:
* primality of all primes;
* the orders, computed from the factorisation of q − 1, certified by 2^k ≡ 1 and
  2^{k/r} ≢ 1, and also confirmed by direct iteration for q < 2·10^7;
* T0 = lcm of the base orders;
* each cluster has one common prime order p_j, coprime to T0, with the p_j distinct;
* T = T0·Π p_j equals the lcm of all orders;
* the histogram is ascending, and sum = M0, first moment = T0·φ(M0), second moment as
  above.

**Cluster distributions.** These use the occupancy formula
δ_C(p − k) = Σ_m e_m · C(p,k) · k! · S2(m,k), where e_m are the coefficients of
Π (x + q − p). They are cross-checked by:
* an independent DP over the number of occupied classes;
* a "grouped literal" enumeration, where each residue's removed class is computed
  from 2^r mod q (used when (p+1)^m ≤ 2·10^6);
* literal enumeration of b mod M_C (used when M_C ≤ `--brute-limit`, default 3·10^6;
  with `--brute-limit=12000000` it also covers [47,178481], [431,9719] and
  [2351,4513]).

Each cluster distribution must sum to M_C, and its first moment must be p·φ(M_C).

**Combination.** The multiplier distribution D(H), H = Π h_j, is a dict convolution. It
is sorted, with prefix sums of D and D·H. For each nu0 the threshold
2 M nu0 H λ_den < T φ(M) λ_num is found by bisection, so S = A/(T φ λ) + B/(2M) with
exact integers A, B. The global first moment Σ_a nu(a) = T·φ(M) is checked.

**λ.** λ = Σ_{k=1}^{200} 1/(k 2^k). Since log 2 = −log(1 − 1/2) = Σ_{k≥1} 1/(k 2^k) and
all terms are positive, λ < log 2. The error is below 1/(201·2^200) < 4·10^-63. This is
not the atanh series. For information the script also prints S(λ) − S(λ + err) ≥ 0.

**Output.** The exact fraction, 40 decimals (truncated), and S rounded UP to 15
decimals.

## Validation results (Mac, Apple M-series, 14 cores, 24 GB)

| target | sha256 of histogram | entries | match | wall | peak RSS | max stored layer |
|---|---|---|---|---|---|---|
| (a) Griego base, 13 primes | f35f2f3d0a14c04ecc0a8f65a0b2e8b469f358a344b3e0038bee83e1e65c7d08 | 2104 | yes | 0.05 s | 0.10 GB | 14,892 |
| (b) Yoo base, 14 primes | 4ac1d42106f662d810a47876ff60fbc93c9e7b2d0f96c4b7c4fa8c3351d0aeee | 4217 | yes | 0.18 s | 0.13 GB | 57,524 |
| (c) Yoo ∪ {43,127,241}, 17 primes | 6a9246781ee53f732f9117abb420513f5045fe88debc3914a318cbbc84b4479a | 4407 | yes | 3.5 s | 0.43 GB | 2,541,964 |
| 20-prime set (list in the task) | 4a7d768934902d6d3e8ce246c0a1665ef678f5135e70260869316aa112b7377a | 4556 | no reference | 441 s | 9.0 GB | 97,499,214 |

(With `--tail=2`: (a) 0.05 s, (b) 0.19 s, (c) 3.7 s and 0.83 GB.)

`evalB.py`:
* (a) S = 0.4902494078111542556868594292410419912779… and rounded up
  **0.490249407811155**. Both match. 0.8 s.
* (b) S = 0.4901800632900604571004505740261343540790… and rounded up
  **0.490180063290061**. Both match. 13 s, with 2,635,776 distinct multiplier
  products.

Robustness matrix (all produce the reference sha256), in `proto/matrix.log`:
* Griego with the given order × tail 1/2/3 × canon tree/affine, plus canon=none;
* Yoo × tail 1/2/3 × canon tree/affine, plus two other explicit orders and
  `--threads=1`;
* (c) × tail 1/2/3 × canon tree/affine.

### 20-prime set {3,5,7,11,13,17,19,29,31,37,41,61,73,257,43,127,241,113,109,97}

It fits on the Mac. The layer sizes (auto order) were:

| step | prime | layer (g) | states |
|---|---|---|---|
| 11 | 61 | Z_12 | 2,767,359 |
| 12 | 13 | Z_12 | 24,154,492 |
| 13 | 5 | Z_12 | 30,642,856 |
| 14 | 7 | Z_12 | 43,310,503 |
| 15 | 3 | Z_12 | 39,069,905 |
| 16 | 73 | Z_36 | 97,499,214 |

The tail then enumerates 19, 37, 109. Results: 4556 entries, all four moment checks
OK, sha256 4a7d7689…377a. The final binary took 441 s wall (5117 s CPU), with peak
RSS 9.0 GB and a peak footprint of 9.1 GB: about 48 s for the stored DP and the rest in
the tail.

The same file was produced by:
* the earlier unfused-tail implementation (894 s);
* the fused-tail build (633 s);
* a different processing order,
  `--order=257,97,17,241,29,113,43,127,61,41,11,31,3,7,13,5,73,19,109,37`
  (429 s, 8.65 GB). Its intermediate layers differ, but the layers after 3 and after
  73 have identical state counts, as they must.

With `--tail=2` the layer after 19 would be stored as well. That run was stopped after
more than 25 minutes on the Mac (memory pressure), so `--tail=3` is the default.

## Notes for the 64-core box

* Memory is roughly 70–90 B per state in the largest layer, plus the parent layer
  while it is being consumed (parent blocks are freed progressively).
* Use `--plan` to inspect the cut moduli. Use `--maxstates=N` to probe an order
  cheaply.
* The order matters a lot. Primes such as 433 (order 72 = 8·9) or 181 (order
  180 = 4·9·5) touch two "blocks" at once and will raise g in the middle stretch.
  Try several explicit `--order=` lists with `--maxstates` before a long run.
* T0 must be ≤ 65535 and M0 < 2^256 (both checked).
