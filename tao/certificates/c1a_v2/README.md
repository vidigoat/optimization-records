# C1a ≤ 1.5027427356 (Sidon set autocorrelation constant, upper bound)

`c1a_n262144_v2.json` is a non-negative step function with n = 262144 equal steps on [−1/4, 1/4], with integer heights
`v_int` (largest 2⁵⁶). Its score 2n · max_k (v∗v)_k / (Σ v)² is, exactly,

    29134419116089766846043850668467284344832 / 19387496227028661445940960928888412813689 = 1.502742735571633181...

The peak is at index 96751, with (v∗v)_96751 = 222277977875440726059294514987695956 and S = Σ v = 278477979215798400234.

Two independent checks agree on the peak index, the peak value and the score:
- `exact1a.py` computes the whole autoconvolution exactly, as one big-integer product (Python standard library only);
- `check1a_fft.py` finds every lag within 10⁻¹⁰ of the peak by float64 FFT, which is about 100× its error bound of
  9·10⁻¹³. It then recomputes those 25 lags exactly as integer dot products (needs numpy).

Run them with:

    python3 exact1a.py c1a_n262144_v2.json
    python3 check1a_fft.py c1a_n262144_v2.json

Comparisons:
- 1.905e-9 below the candidate 1.5027427375* (exactly 1.50274273747657…) of A. Canli
  ([teorth/optimizationproblems#346](https://github.com/teorth/optimizationproblems/pull/346));
- 6.17e-7 below our `../c1a` (1.502743352178).

## How it was found

The run starts from Canli's #346 certificate (heights as published, 1.5027427374766). That certificate was itself
refined from our `../c1a` 262144-cell function.

`push1a.py` runs L-BFGS-B with v ≥ 0 on the soft maximum 2n‖v∗v‖_p / S². It restarts the p-continuation at
p₀ = 10⁸ and doubles p up to about 8·10¹¹, with 3000 iterations per stage. This is the schedule #346 describes, with
longer stages. The command was:

    python push1a.py arh_vint.json out.json --split 1 --iters 3000 --rounds 4 \
      --plist 1e8,2e8,4e8,8e8,1.6e9,3.2e9,6.4e9,1.28e10,2.56e10,5.12e10,1.024e11,2.048e11,4.096e11,8.192e11

The gain came at p = 2·10⁹. `search_log.txt` is the run's log; its `int=` column is the exact score of the saved
integer heights.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by
the checkers above. It has not been verified line by line by a human expert.
