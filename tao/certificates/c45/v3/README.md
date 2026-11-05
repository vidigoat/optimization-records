# C45 (Romanoff's constant): upper bound 0.489440909905826 (34-prime base, 59 primes)

    upper density of {n odd : n = p + 2^k}  <  0.489440909905826

Exact certified rational sum: S = 0.489440909905825679514736174270792438529576472144417336768449...
(numerator and denominator in `data/xA5.txt.S.json`).

Comparisons:
- 3.29e-5 below our 57-prime certificate `../v2` (0.489473818147233).
- 4.08e-5 below the 54-prime bound 0.489481743486002 of R. Salavatov ([teorth/optimizationproblems#290](https://github.com/teorth/optimizationproblems/pull/290)).

## Inequality and method

The inequality, the group decomposition, the generator `c45x.cpp`, the evaluator `c45eval.py` and `mkspec.py` are the same as in `../v2`, unchanged byte for byte. See `../v2/README.md` for the full description.

## Prime set (59 primes)

* **Base, 34 primes**, with T0 = 12252240 = 2^4 3^2 5 7 11 13 17. It is the 32-prime base of `../v2`
  **+ 331 and 673**, of orders 30 and 48. Both orders divide T0, so T0 is unchanged.
  * 331 joins the B group with L = 60; lcm(60, 30) = 60.
  * 673 joins the B group with L = 48.
  * The group lcms, and so the pairwise-coprime quotient structure the generator needs, are unchanged.
* **Clusters:** the 25 cluster primes of Salavatov's certificate and of `../v2`, unchanged.

The groups:

    B:3,7,5,13   B:37,109,19,73   B:31,11,41,61,151,331   B:17,257,97,241,673
    A:29,113,43,127   A:8191,53,157,1613   A:23,89,397,683,2113   A:137,953

The smallest valid cut is C = 1528627, from `mkspec.py`.
- Only bins with nu0 <= C are written: 1,495,742 lines.
- The evaluator checks (C + 1) * H_min * kappa >= 1/2, so each class above C contributes exactly 1/(2M).

## Files

* `data/xA5.txt.xz`: the base histogram, lines `nu0 delta0(nu0)` for nu0 <= 1528627.
  * sha256 of the decompressed file: `612e6348a7af5344c5e7622e51ccaad58f0982d7d9e81b08da5375dc5a634b8d`.
  * sha256 of the .xz: `07bf15534bd824fb917c914adbfa6667f2965425b43c12d2dc1ce95455506527`.
* `data/xA5.txt.S.json`: the exact S.
* `specs/specA5.json`: this certificate's spec. `specs/spec54.json` is read by `mkspec.py`, for the clusters.
* `logs/xA5.err`: the generator log (`/usr/bin/time -v`). The run took 5 h 56 min on 64 cores and peaked at 169 GB.
* `logs/xA5.eval.remote`, `logs/xA5.eval.local`: two evaluations of S that agree, one on the generating machine and one run
  locally from the downloaded file. The local run took 13 s.

## Reproduce

    clang++ -O3 -std=c++17 -pthread c45x.cpp -o c45x      # or g++
    ./c45x 1528627 64 B:3,7,5,13 B:37,109,19,73 B:31,11,41,61,151,331 B:17,257,97,241,673 \
        A:29,113,43,127 A:8191,53,157,1613 A:23,89,397,683,2113 A:137,953 > xA5.txt
    python3 c45eval.py specs/specA5.json xA5.txt 1528627

Generation needs about 170 GB of RAM. Evaluating the published histogram needs only Python:

    xz -dk data/xA5.txt.xz && python3 c45eval.py specs/specA5.json data/xA5.txt 1528627

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction (generator, evaluator and text).
The value is recomputed by the code above. It has not been verified line by line by a human expert.
