# C1c ≥ 0.4103596512 (Barnard–Steinerberger autocorrelation constant, lower bound)

The certificate `c1c_K98304_v3.json` (sha256 `e0786a59bd9d8598d9f4a86fcee40575ea108a97b184b332ef05bb6764caec66`) is a
step function f = Σ w_i 1_[i/K,(i+1)/K) on K = 98304 cells. It has 98368 non-negative integer weights, the largest
being 2⁴⁰. It is scored exactly with the page's rule, K · min_{0≤k≤K} Σ_i w_i w_{i+k} / (Σ w_i)²:

    640823574552569898220412928 / 1561614482809354047911104009 = 0.41035965124998351121...   (minimum at k = 98301)

with A_min = 26075178000999751717953 and S = Σ w_i = 79034536319494.

Comparison with earlier bounds:
- Previous bound: 0.4103540676 (`../c1c_v2`, ours; independently replayed in teorth/optimizationproblems#270).
- Before that: 0.4103538918 (J. Julsgaard, #318), refined from our 26 Oct certificate.
- This bound is 5.76·10⁻⁶ above #318 and 5.58·10⁻⁶ above `c1c_v2`.

Two independently written exact scorers agree:
- `exact.c`, which sums every lag in 128-bit integers;
- `exact2.py`, which gets all lags from one big-integer product, standard library only.

Their outputs are in `exact_c_output.txt` and `exact2_output.txt`.

    cc -O3 -o exact exact.c && ./exact c1c_K98304_v3.json
    python3 exact2.py c1c_K98304_v3.json

## How it was found

The search continued the cyclic re-softening of `../c1c_v2`, using `push.py` and `refine.py`: L-BFGS-B on a soft minimum, with
the sharpness cycled through p = 3·10⁷, 10⁸ and 3·10⁸. The change that mattered was **longer inner runs**: 20,000
L-BFGS iterations per stage instead of 5,000. That lane gained about 3·10⁻⁶ in a few hours, while the 5,000-iteration
lanes stalled near 0.4103565. The best lane (`--plist 1e7,3e7,1e8,3e8 --iters 20000`) was then restarted from the
leading point.

Found with Claude (Anthropic) agents under V. Patankar's direction; the value is recomputed exactly by the two scorers.
