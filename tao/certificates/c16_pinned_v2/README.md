# C16 to 18 digits: 2.156225528154212956 ≤ C16 ≤ 2.15622552815421296

This continues `../c16_pinned` (14 digits) with the same two bounds, computed more precisely.
The two sides now meet to within 4.4·10⁻¹⁹.

| side | value | how |
|---|---|---|
| upper | **C16 ≤ 2.15622552815421296** | `verify2.py certificate.json`: largest interval-end value 2.15622552815421295831 ± 5·10⁻²¹ |
| lower | **C16 ≥ 2.156225528154212956** | `lower2.py 0x1.f42a1bd5af466p-4`: L ≥ 2.1562255281542129561 ± 3.2·10⁻²⁰ at δ = 3.92883616571834… |

## What changed from `../c16_pinned`

1. **G(μ) and H(μ) to ~10⁻²¹.** The row sums are the same exact (π/s)coth(πs) closed forms. The convex tail of rows is still
   bracketed between trapezoid and midpoint sums. With N = 200,000 explicit rows instead of 20,000, the bracket width
   drops from about 6·10⁻¹⁷ to about 6·10⁻²¹. Arb works at 192 bits.
2. **The lower bound is maximised over μ exactly.** The golden-section search compares Arb values, not floats. The
   earlier μ = 0.122110476195 was 2·10⁻¹⁷ short of the optimum, because float comparisons cannot resolve differences that small.
   The optimum is μ = 0x1.f42a1bd5af466p-4 = 0.12211047051364829.
3. **A finer upper certificate near the peak.** There are 383 δ-intervals covering [1, 2¹⁸]; the tail δ ≥ 2¹⁸ is handled exactly as before. On each
   interval, μ is chosen in double precision by balancing the two endpoint values:
   μ G(μ) = π (Λ(hi) − Λ(lo)) / (hi − lo), with Λ = log δ + log(1 + log δ), clipped to the stationary values μ*(lo) and μ*(hi).
   The difference Λ(hi) − Λ(lo) is taken in 50-digit arithmetic, which avoids cancellation on narrow intervals.
   Every endpoint is then checked in Arb. Only the checked bound is a claim; the double-precision choice of μ is not.

The inequality and the proof structure are unchanged; see `../c16_pinned/README.md` and the docstring of `verify2.py`.

## Reproduce (python-flint)

    python verify2.py certificate.json --crosscheck 21     # about 7 min on one core
    python lower2.py 0x1.f42a1bd5af466p-4                  # about 2 s

The outputs are in `verify2_output.txt` and `lower2_output.txt`. The brute-force lattice-sum cross-check is a sanity
check only: it agrees with the enclosed G(μ) to the accuracy of the brute-force sum.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The values are recomputed by the
scripts above. They have not been verified line by line by a human expert.
