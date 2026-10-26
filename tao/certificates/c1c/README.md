# C1c ≥ 0.410343281299 (Barnard–Steinerberger autocorrelation constant, lower bound)

Step function f = Σ w_i 1_[i/K,(i+1)/K) with K = 98304 and 98368 non-negative integer weights (`c1c_K98304.json`),
scored exactly as K · min_{0≤k≤K} Σ_i w_i w_{i+k} / (Σ w_i)²:

    2594500369901311586790309888 / 6322755819680014072464230569 = 0.41034328129923...   (minimum at k = 97663)

Previous best 0.4103333527 (Koshta, K = 24576). Found by splitting each cell of that certificate in two, twice,
and re-optimising with L-BFGS on a soft minimum (`refine.py`).

Two independent exact scorers agree on every digit:
- `exact.c`: every lag summed in 128-bit integers (`cc -O3 -o exact exact.c && ./exact c1c_K98304.json`);
- `exact2.py`: all lags at once by one big-integer product (Kronecker substitution), standard library only.
Both also reproduce the K = 24576 record value exactly.
