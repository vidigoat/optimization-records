# C1c ≥ 0.4103604995 (Barnard–Steinerberger autocorrelation constant, lower bound)

`c1c_K98304_v4.json` (sha256 `58dacf407fba7b293b044815a359b7b444d729c76b6e296d6c10e77f250416f6`) is a step function on
K = 98304 cells, with 98368 non-negative integer weights. Under the page's exact rule,
K min_{0≤k≤K} Σ_i w_i w_{i+k} / (Σ_i w_i)², its value is

    5282641762424933775753216 / 12873173145647069582596129 = 0.410360499517650353068...

The minimum is at k = 98288, with A_min = 26009100474178751093186 and S = 78934249869706. Both exact scorers agree:
`exact.c` (128-bit integers over every lag) and `exact2.py` (one big-integer product); their outputs are in this folder.

Comparisons:
- 6.28e-7 above the candidate 0.4103598711* of A. Canli ([teorth/optimizationproblems#344](https://github.com/teorth/optimizationproblems/pull/344));
- 8.48e-7 above our `../c1c_v3` (0.4103596512*).

## How it was found

The search starts from Canli's #344 certificate (sha256 `17c7640a…9d35`, 0.4103598711).
- `push.py` runs L-BFGS-B on a soft minimum, re-softened to p = 128 and then doubled, with 5000 iterations per stage.
- The gains came at p = 2·10⁶ and p = 4·10⁶; `search_log.txt` is the run's log.
- This is the p-continuation restart that #344 and #312 describe, applied with longer stages.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checkers above. It has not been verified line by line by a human expert.
