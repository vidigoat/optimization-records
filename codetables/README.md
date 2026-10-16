# Record linear codes (codetables.de)

[codetables.de](http://www.codetables.de), maintained by Markus Grassl, is the standard table of the
best known linear error-correcting codes. For each field size q, length n and dimension k it lists
the largest known minimum distance d, and the best proven upper limit.

Seven new codes over GF(5) beat their entries:

| Code | Listed lower bound | New | Upper bound | Entry last changed |
|---|---|---|---|---|
| [27, 10] | 12 | **13** | 14 | 2001 |
| [64, 10] | 38 | **39** | 44 | 2003 |
| [81, 10] | 50 | **51** | 57 | 2008 |
| [82, 10] | 51 | **52** | 58 | 2003 |
| [67, 12] | 37 | **38** | 45 | 2003 |
| [68, 12] | 38 | **39** | 45 | 2003 |
| [69, 12] | 38 | **39** | 46 | 2003 |

They were found by Construction X on pairs of nested constacyclic codes (a smaller code inside a
larger one, both invariant under a twisted cyclic shift, joined to a short auxiliary code).

The files `gf5_<n>_<k>.txt` hold the generator matrices (one row per line, entries 0–4). To check:

    mkdir -p bin && cc -O3 -o bin/mindist codetables/mindist.c && python3 codetables/check_all.py

[`mindist.c`](mindist.c) computes the exact minimum distance by enumerating all codewords (5^10 ≈ 9.8 million for k = 10, 5^12 ≈ 244 million for k = 12)
 (a Gray-code walk, one row added per step). [`verify.py`](verify.py) also checks that the
rows are independent. Each distance was confirmed a second time with a separate NumPy enumeration.
