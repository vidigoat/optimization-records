# C1c ≥ 0.4103608339 (Barnard–Steinerberger autocorrelation constant, lower bound)

`c1c_K98304_v5.json` is a step function on K = 98304 cells, with 98368 non-negative integer weights. Under the page's
exact rule, K min_{0≤k≤K} Σ_i w_i w_{i+k} / (Σ_i w_i)², its value is

    2556315663995378551725293568 / 6229433836731244234885554241 = 0.410360833904731...

The minimum is at k = 45081, with A_min = 26004187662713404863742 and S = 78926762487329. Both exact scorers agree:
`exact.c` (128-bit integers over every lag) and `exact2.py` (one big-integer product); their outputs are in this folder.

Run them with:

    python3 exact2.py c1c_K98304_v5.json
    cc -O3 -o exact exact.c && ./exact c1c_K98304_v5.json

Comparisons:
- 9.62e-7 above the candidate 0.4103598711* of A. Canli ([teorth/optimizationproblems#344](https://github.com/teorth/optimizationproblems/pull/344));
- 3.34e-7 above our `../c1c_v4` (0.4103604995*).

## How it was found

The run starts from Canli's #344 certificate (0.4103598711). `push.py` is the same script as in `../c1c_v4`: L-BFGS-B on
a soft minimum, re-softened to p = 100 and then doubled, with 5000 iterations per stage. Here the run continues to
p = 4·10⁹, and it kept improving through every stage from p ≈ 10⁶ onward. This run saved its weights as integers at
scale 2⁴⁰ (max weight 2⁴⁰). The copy of `push.py` here saves at 2⁵⁶, which loses less to rounding. The `int=` column
of `search_log.txt` is the exact score after rounding.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checkers above. It has not been verified line by line by a human expert.
