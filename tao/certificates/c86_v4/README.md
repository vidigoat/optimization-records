# C86 ≤ 1.8054838 (Schur–Siegel–Smyth trace constant, upper bound)

`certificate.json` (sha256 `7b8c14bd743593e61a2f75c0e020d100146da5e1162be71e212755537e70153f`) is in R. Pinho's format
([#298](https://github.com/teorth/optimizationproblems/pull/298)). Its bound:

    T/m = 260197628230434684107589673124053/144115188075815087401774532788224 = 1.80548373633979...

Comparisons:
- 1.03·10⁻⁵ below the 1.8054939963* of R. Pinho's #298 (10 Nov update). That update started from our `../c86_v3`
  density.
- 1.18·10⁻⁵ below our `../c86_v3` (1.8054955).

## What changed from `../c86_v3`

This continues the same search, from its saved state after one more pricing round (5 Nov). The run also uses the 31
auxiliary polynomials of Pinho's #298 round-5 certificate as extra columns; this is data only, and none of his code
was read or run. The cutting planes were run until no refined local minimum fell below half the margin. The converged
state was then exported with `search/export86.py`. The density has 7873 knots, and the certificate has 37 auxiliaries.
The search log is in `receipts/search_log.txt`.

## Check

`checkers/check86.py` (Arb, 128 bits) gives the verdict `VALID`; its output is in `receipts/check86.txt`. It checks:
- the structure and the degree condition (Σ b deg Q / m = 0.8992);
- that every ∫ log|Q| dμ is positive (smallest 8.4·10⁻⁸);
- that f > 0 on all 7871 support cells, over 112,583 intervals with 0 failures and smallest certified lower bound
  1.4·10⁻¹⁰;
- that the bound equals T/m exactly.

Run it with:

    python3 checkers/check86.py certificate.json --procs 8 --beat 1.8054939963

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checker above. It has not been verified line by line by a human expert.
