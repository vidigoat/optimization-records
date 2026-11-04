# C1c ≥ 0.4103540676 (Barnard–Steinerberger autocorrelation constant, lower bound)

The certificate `c1c_K98304_v2.json` is a step function f = Σ w_i 1_[i/K,(i+1)/K) on K = 98304 cells, with 98368 non-negative integer weights. The largest weight is 2⁴⁰.

It is scored exactly with the page's rule, K · min_{0≤k≤K} Σ_i w_i w_{i+k} / (Σ w_i)²:

    645458795582745374892957696 / 1572931393892132517694378321 = 0.41035406762757335811...   (minimum at k = 93282)

Previous bounds:
- 0.4103538918 (J. Julsgaard, github.com/Julsgaard/c1c-certificate, commit 91e47b1). That certificate was itself refined from our 26 Oct certificate, 0.4103432812.
- This bound is 1.76·10⁻⁷ above it.

Two independently written exact scorers agree:
- `exact.c` sums every lag in 128-bit integers: `cc -O3 -o exact exact.c && ./exact c1c_K98304_v2.json` (see `exact_c_output.txt`).
- `exact2.py` gets all lags from one big-integer product (Kronecker substitution), using only the standard library (see `exact2_output.txt`).

Both report A_min = 26263785627553115840371 and S = Σ w_i = 79320398231278.

## How it was found

The search starts from Julsgaard's K = 98304 certificate. `push.py` runs L-BFGS-B on the soft minimum in `refine.py`, with FFT gradients, in repeated cycles:
1. a softer sharpness (p = 10⁸),
2. then sharper ones (3·10⁸ and 10⁹),
3. each cycle restarting from the best point so far.

Re-softening lets the iterate leave the corner it reached at high p. Each cycle gained about 3·10⁻⁸, while one long run at fixed high sharpness stalls.

At this point about 78% of all lags lie within 10⁻⁷ of the minimum, and 29% lie within 10⁻⁸. So the step function is very close to a local optimum for this grid.

The search was done with Claude (Anthropic) agents under V. Patankar's direction.
