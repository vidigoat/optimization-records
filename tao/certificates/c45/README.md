# C45 (Romanoff's constant): certified upper bound 0.489596115074071

Let A = {n odd : n = p + 2^k, p prime, k >= 0}. This certificate proves

    upper density of A  <  0.489596115074071

(exact value of the certified rational sum: S = 0.48959611507407010687293152717386177650...).
That is 5.84e-4 below the current record 0.490180063290061 [Y2026, Yoo] and 6.53e-4 below
0.490249407811155 [G2026, Griego].

## Method (same obstruction inequality as on the optimization-constants page)

P = finite set of odd primes, M = prod P, T = lcm ord_q(2). For a mod M,
nu(a) = #{r mod T : gcd(a - 2^r, M) = 1}, delta(nu) = #{a mod M : nu(a) = nu}, and

    dbar(A) <= S = sum_nu delta(nu) * min( 1/(2M), nu / (T phi(M) lambda) ),   lambda < log 2 rational.

P = base primes (orders dividing T0 = 5040) ∪ clusters (all primes of one prime order p, p coprime to T0
and to the other cluster orders). By CRT, nu(a) = nu0(a mod M0) * prod_clusters h_c, so S is computed exactly
from the base histogram delta0 and the cluster h-distributions.

### Prime set (39 primes)

* Base, 21 primes (T0 = 5040):
  3, 5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 61, 73, 257 (Yoo's base) **+ 43, 127, 241, 113, 109, 97, 151** (new;
  orders 14, 7, 24, 28, 36, 48, 15).
* Clusters (order: primes): 11: 23, 89 | 13: 8191 | 23: 47, 178481 | 29: 233, 1103, 2089 | 37: 223 |
  43: 431, 9719 | 47: 2351, 4513 | 53: 6361 | 73: 439 | 79: 2687 | 83: 167 | 113: 3391 | 131: 263 |
  179: 359, 1433 | 191: 383 | 233: 1399 | 239: 479, 1913, 5737, 176383 | 251: 503
  (union of Griego's and Yoo's clusters).
* T = 5040 * 11*13*23*29*37*43*47*53*73*79*83*113*131*179*191*233*239*251.

The gain comes almost entirely from enlarging the base: single-prime marginal gains measured on Yoo's base were
43: 1.74e-4, 109: 9.6e-5, 241: 7.2e-5, 113: 7.2e-5, 97: 7.1e-5, 127: 6.8e-5, 151: 4.7e-5, then 331: 2.2e-5,
71: 1.4e-5, 433/673/337/181: ~1e-5 each. Adding every known Mersenne factor < 1e7 of 2^p-1 (p <= 1200) as
extra clusters changes S by only ~1e-6, so clusters were not extended further.

## Files

* `data/base_final.txt` — exact 21-prime base histogram, lines `nu delta0(nu)` (4559 entries,
  sha256 238174d0239c9cec3718238416d62d201405471c2a63127ba2397ebba1fcc41d, sum = M0).
* `data/base_20.txt` — 20-prime base (without 151), certifies 0.489640182196384 (intermediate result).
* `data/base_griego13.txt`, `data/base_yoo14.txt` — Griego's and Yoo's seed histograms (recomputed).
* `specs/A_*.json`, `specs/B_*.json` — certificate specs (base primes, histogram file, clusters).
* `checkerA/` — histogram generator `basedp3.cpp` + exact evaluator `exact_eval_a.py` (+ `brute.py`).
* `checkerB/` — independently written generator `histB.cpp` + evaluator `evalB.py` (written from the
  mathematical statement only, without access to checker A's code; see `checkerB/README.md`).
* `run_checkerA.sh`, `run_checkerB.sh` — rebuild everything from first principles.
* `logs/` — run logs.

## Verification performed

| | checker A | checker B |
|---|---|---|
| histogram algorithm | DP over primes; state = counts mod L = lcm(remaining orders); states merged under affine maps x->ux+t of Z_L; last 3-4 primes enumerated as patterns with the final prime applied in closed form | DP with state on Z_g, g = gcd(lcm processed, lcm remaining); symmetry = automorphisms of nested-block trees (wreath products) or affine group; different processing order (9-block last); built-in 0th-3rd moment checks |
| evaluation | Python Fractions; lambda = 2 sum_{j<60} 1/((2j+1)3^(2j+1)) | Python stdlib integers/Fractions; lambda = sum_{k<=200} 1/(k 2^k); clusters via occupancy (Stirling) formula cross-checked by brute force; checks sum, 1st, 2nd moments of the histogram |
| Griego (13-prime seed, 36 primes) | 0.4902494078111542556868594292... -> 0.490249407811155 | identical |
| Yoo (14-prime seed, 38 primes) | 0.4901800632900604571004505740... -> 0.490180063290061 | identical |
| 20-prime base (38 primes total) | 0.4896401821963830395194824005... -> 0.489640182196384 | identical; histogram sha256 identical (4a7d7689...) |
| **21-prime base (39 primes total)** | **0.4895961150740701068729315271... -> 0.489596115074071** | identical; histogram sha256 identical (238174d0...); 0th-3rd moment checks OK |

Both generators reproduce Griego's and Yoo's seed histograms bit-for-bit (sha256 f35f2f3d..., 4ac1d421...),
and both match direct enumeration over a mod M0 on small prime sets.

Runtimes (checker B histogram for the 21-prime base: 20.2 min, 699.5M stored states, 49 GB, same machine). Runtimes (checker A, c7a.16xlarge, 64 vCPU): 20-prime base 5 min (149M states, 20 GB); 21-prime base 31.7 min
(187.5M states before the final stage, 23 GB). Exact evaluation: ~15 s (A), ~1 min (B) on a Mac
(9,676,800 multiplier entries).
