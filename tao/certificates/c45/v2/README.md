# C45 (Romanoff's constant): upper bound 0.489473818147233 (32-prime base, 57 primes)

    upper density of {n odd : n = p + 2^k}  <  0.489473818147233

Exact certified rational sum: S = 0.489473818147232759198284099287928974701900478242653238729654...
(numerator and denominator in `data/xA3.txt.S.json`). That is 7.93e-6 below the 54-prime bound
0.489481743486002 of R. Salavatov ([teorth/optimizationproblems#290](https://github.com/teorth/optimizationproblems/pull/290)),
which itself extended the 21-prime base of our earlier certificate (`../README.md`, 0.489596115074071).

## Inequality

Same finite obstruction inequality as on the optimization-constants page and in `../README.md`: for a finite set P of odd
primes, M = prod P, T = lcm ord_q(2), lambda < log 2 rational,

    dbar(A) <= S = sum_nu delta(nu) * min( 1/(2M), nu / (T phi(M) lambda) ).

P = base (orders dividing T0) ∪ clusters (orders coprime to T0 and to each other), so nu = nu0 * prod_c h_c and S is
computed exactly from the base histogram and the cluster occupancy laws.

## Prime set (57 primes)

* Base, 32 primes, T0 = 12252240 = 2^4 3^2 5 7 11 13 17: Salavatov's 29-prime base
  3,5,7,11,13,17,19,29,31,37,41,61,73,257, 43,127,241,113,109,97,151, 8191,53,157,23,89,397,683,137
  **+ 953, 1613, 2113** (orders 68, 52, 44; all divide T0, so T0 is unchanged).
* Clusters: those of Salavatov's certificate (the clusters of Griego and Yoo without orders 11 and 13; 25 primes),
  unchanged.

## Generator (`c45x.cpp`)

The base primes are split into groups whose order-lcms L_G satisfy: B groups share only r mod 12, A groups only r mod 4,
and the quotients L_G/12 (B) and L_G/4 (A) are pairwise coprime with product T0/12. Then, for a mod M0,

    nu0(a) = sum_{s mod 4} sigma_s(a) beta_s(a),   sigma_s = sum_{c mod 12, c = s (4)} prod_{G in B} alpha_G(c),
    beta_s = prod_{G in A} alpha_G(s),

where alpha_G(c) counts residues r mod L_G with r = c (mod 12 or 4) not killed by the primes of G, and the alpha_G are
independent over a by CRT. The groups used here:

    B:3,7,5,13   B:37,109,19,73   B:31,11,41,61,151   B:17,257,97,241
    A:29,113,43,127   A:8191,53,157,1613   A:23,89,397,683,2113   A:137,953

L = 12, 36, 60, 48 | 28, 52, 44, 68; quotients 1, 3, 5, 4 | 7, 13, 11, 17. Per-group laws are exact (128-bit counts);
the B side is stored up to the 48 affine maps of Z/12 and then sigma up to the affine maps of Z/4 (128-bit orbit
weights); the A side is the full product law (256-bit counts); the histogram uses 384-bit counts. Only bins
nu0 <= C = 1535541 are written; a pair (sigma, beta) is skipped only when sum(sigma) * min(beta) > C, which implies
nu0 > C. The evaluator checks (C + 1) * H_min * kappa >= 1/2, so every class above C contributes exactly 1/(2M), and
gets their total count as M0 minus the written counts.

## Evaluator (`c45eval.py`)

Python Fractions. lambda = 2 sum_{j<60} 1/((2j+1) 3^(2j+1)) < log 2; cluster laws by the exact occupancy recursion;
checks that every listed order is the true order of 2, that cluster orders do not divide T0, and the cut condition.

## Files

* `data/xA3.txt.xz`: the base histogram, lines `nu0 delta0(nu0)` for nu0 <= 1535541 (1,479,287 lines). sha256 of
  the decompressed file `272bcc40d95a73215043ed38c552dac53ce82d68d2adc280ce96bd1e68d36bdb`; of the .xz
  `51aa143f287f202605e46142382aecfc88933880bd679e88ca56b279a8fdcf00`.
* `data/xA3.txt.S.json`: exact S. `data/x54.txt.S.json`: exact S of the 54-prime replay.
* `specs/specA3.json` (this certificate), `specs/spec54.json`, `specs/spec21.json` (replays).
* `logs/`: generator logs (`/usr/bin/time -v`); the 57-prime run took 2 h 26 min on 64 cores and 50 GB.
* `mkspec.py`: writes a spec for a group list and prints the smallest valid cut C.

## Validation

| run | primes | result | reference |
|---|---|---|---|
| 21-prime base | 39 | histogram sha256 `238174d0...` | **identical** to `../data/base_final.txt` (two independent checkers, Oct 2026) |
| 29-prime base | 54 | 1,434,945 bins, C_min = 1538836, S = 0.489481743486001047529993730534590109425... | Salavatov's #290: same bin count, same cut, S agrees in all 39 published digits |
| **32-prime base** | **57** | **S = 0.489473818147232759198...** | this certificate |

The 57-prime value was evaluated twice (on the generating machine and again locally from the downloaded file).

Reproduce:

    clang++ -O3 -std=c++17 -pthread c45x.cpp -o c45x      # or g++
    ./c45x 1535541 64 B:3,7,5,13 B:37,109,19,73 B:31,11,41,61,151 B:17,257,97,241 \
        A:29,113,43,127 A:8191,53,157,1613 A:23,89,397,683,2113 A:137,953 > xA3.txt
    python3 c45eval.py specs/specA3.json xA3.txt 1535541

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction (generator, evaluator and text);
the value is recomputed by the code above. It has not been verified line by line by a human expert.
