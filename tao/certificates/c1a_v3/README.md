# C1a ≤ 1.5027427341 (Sidon set autocorrelation constant, upper bound)

`c1a_n262144_v3.json` is a non-negative step function with n = 262144 equal steps on [−1/4, 1/4], with integer heights
`v_int` (largest 2⁵⁶). Its score 2n · max_k (v∗v)_k / (Σ v)² is, exactly,

    116537598371962553887369519059399403896832 / 77549933019818026777245548252121251545225 = 1.502742734054730...

The peak is at index 262143, with (v∗v)_262143 = 222277828926015002989520109289931114 and S = Σ v = 278477886051689976515.

Two independent checks agree on the peak index, the peak value and the score:
- `exact1a.py` computes the whole autoconvolution exactly, as one big-integer product (Python standard library only);
- `check1a_fft.py` finds every lag within 10⁻¹⁰ of the peak by float64 FFT (79,895 lags, because the plateau is very
  flat), then recomputes each of them exactly as an integer dot product, splitting every height into three 19-bit limbs
  so the partial sums fit in int64 (needs numpy; about 4 minutes). This replaces the pure-Python dot products of
  `../c1a_v2/check1a_fft.py`, which would take about an hour on this many lags.

Their outputs are in `exact1a_output.txt` and `check1a_fft_output.txt`. Run them with:

    python3 exact1a.py c1a_n262144_v3.json
    python3 check1a_fft.py c1a_n262144_v3.json

Comparisons:
- 3.42e-9 below the candidate 1.5027427375* (exactly 1.50274273747657…) of A. Canli
  ([teorth/optimizationproblems#346](https://github.com/teorth/optimizationproblems/pull/346));
- 1.517e-9 below our `../c1a_v2` (1.502742735571633).

## How it was found

This is the run of `../c1a_v2`, continued: the same command, from Canli's #346 certificate, four stages further
(p up to 2.56·10¹⁰). Every stage from p = 8·10⁸ onward improved the bound (the log prints p to one significant digit). `search_log.txt` is the run's log up to
this certificate; its `int=` column is the exact score of the saved integer heights. The command was:

    python push1a.py arh_vint.json out.json --split 1 --iters 3000 --rounds 4 \
      --plist 1e8,2e8,4e8,8e8,1.6e9,3.2e9,6.4e9,1.28e10,2.56e10,5.12e10,1.024e11,2.048e11,4.096e11,8.192e11

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checkers above. It has not been verified line by line by a human expert.
