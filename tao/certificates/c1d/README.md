# C1d ≤ 1.4546733458 (signed third autocorrelation constant, upper bound)

`c1d_n163840.json` is a signed step function with n = 163840 equal cells on [−1/4, 1/4], with integer cell values
`v_int` (largest |v_i| = 10¹²). By the knot formula (f∗f)(mh) = h (v∗v)_{m−1}, with h = 1/(2n), its value of
max_t |f∗f| / (∫f)² is 2n · max_m |(v∗v)_m| / (Σ v)², which is exactly

    646233848324886328286774558720 / 444246710235775351931155445761 = 1.454673345767513...

The peak is at index 165229, with |(v∗v)_165229| = 1972149195327411890523604 and S = Σ v = 666518349511681.

Two independent checks agree on the peak index, the peak value and the score:
- `exact1d.py` computes the whole autoconvolution exactly (Python standard library only). It writes v = P − Q with
  P, Q ≥ 0, so v∗v = P∗P − 2 P∗Q + Q∗Q, and each product is one big-integer multiplication;
- `check1d_fft.py` uses a float64 FFT to find every lag within max(100 × its error bound, 10⁻¹¹ × peak) of the
  peak (121,920 lags, because the plateau is very flat). It recomputes each of them exactly as an integer dot product,
  splitting every value into two halves so the partial sums fit in int64 (needs numpy).

Their outputs are in `exact1d_output.txt` and `check1d_fft_output.txt`. Run them with:

    python3 exact1d.py c1d_n163840.json
    python3 check1d_fft.py c1d_n163840.json

Comparison: 1.653e-9 below the candidate 1.454673347420953* (exactly
167771889099128792433071411029129035776/115333032942810266163095027979519535125) of A. Canli
([teorth/optimizationproblems#345](https://github.com/teorth/optimizationproblems/pull/345)).

## How it was found

The run starts from Canli's #345 certificate of 9 Nov (`c1d_upper_n163840.json`, sha256
`66d69522813cd2cf4104a52c97ebf98e7a6629ab59fe26d1c8bf495c944c444c`, 1.454673347420953). That certificate was refined
from Julsgaard's 81,920-cell certificate in
[#312](https://github.com/teorth/optimizationproblems/pull/312), which was itself refined from
[#269](https://github.com/teorth/optimizationproblems/pull/269).

`push1d.py` runs L-BFGS on the soft maximum of |2n (v∗v)_k / S²|. It restarts the p-continuation at p₀ = 10⁸ and
doubles p, with 3000 iterations per stage. This is the schedule that improved C1a in `../c1a_v2`. The command was:

    python push1d.py arh2.json out.json --split 1 --iters 3000 --rounds 4 \
      --plist 1e8,2e8,4e8,8e8,1.6e9,3.2e9,6.4e9,1.28e10,2.56e10,5.12e10,1.024e11,2.048e11,4.096e11,8.192e11

The first gain came at p = 1.6·10⁹, and every later stage of the first round improved the bound. `search_log.txt` is
the log of that round; its `int=` column is the exact score of the saved integer cell values.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checkers above. It has not been verified line by line by a human expert.
