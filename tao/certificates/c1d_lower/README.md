# C1d ≥ 1.152 (signed autocorrelation constant, lower bound)

**Theorem.** Let f: ℝ → ℝ be integrable, supported on [−1/4, 1/4]. Then

    max_t |(f∗f)(t)|  ≥  1.152 · (∫ f)².

The previous lower bound was the trivial 1. The proof is half a page. Its only computer-assisted step is a
finite interval-arithmetic check of one explicit inequality on the boundary of a rectangle.

## Proof

Scale so that ∫ f = 1, and put g = f∗f. Then g is supported on [−1/2, 1/2] and ∫ g = 1. Let M = ‖g‖∞ ≥ 1, and
suppose M < ∞. Use the Fourier transform F(z) = ∫ f(x) e^{−2πizx} dx. F is entire because f has compact support,
and the transform of g is F².

**1. A non-negative remainder.** Put h = (M − g)·1_{[−1/2,1/2]}. Since |g| ≤ M,

    0 ≤ h ≤ 2M,      ∫ h = M − 1,      ĥ(z) = M·S(z) − F(z)²,     where S(z) = sin(πz)/(πz).

**2. A bound on the remainder.** For z = x + iy we have |ĥ(z)| ≤ ∫ h(t) e^{2π|y||t|} dt. Over all h with
0 ≤ h ≤ 2M and ∫ h = M − 1 on [−1/2, 1/2], the right side is largest when the mass sits at density 2M on an
end interval of length (M−1)/(2M). Writing out that integral,

    |M S(z) − F(z)²|  ≤  B_M(|y|) := (M/(π|y|)) (e^{π|y|} − e^{π|y|/M}),      B_M(0) = M − 1.

The ratio B_M(y)/M = (e^{πy} − e^{πy/M})/(πy) increases with M.

**3. Rouché.** Let R = [0.87, 1.23] × [−0.34, 0.34]. The only zero of S in R is z = 1, and it is simple. The
script `check_c1d_lower.py` covers ∂R by 80,000 small boxes. On each box it checks, in interval arithmetic,
that M₀|S(z)| > B_{M₀}(|Im z|) for M₀ = 1.152. Its output is in `check_output.txt`, where the smallest certified
ratio squared is 1.0107. By the monotonicity in step 2, the same strict inequality holds for every M in
[1, 1.152]. So if M ≤ 1.152, then on ∂R

    |F² − M S| < |M S|.

By Rouché's theorem, F² and M·S have the same number of zeros in R, counted with multiplicity. That number is 1.
But every zero of F² has even multiplicity, so the count for F² is even. This is a contradiction, so M > 1.152.
∎

**Remarks.**
- The argument uses only that the transform of f∗f is a square of an entire function. Realness of f is not used,
  so the bound also holds for complex f.
- On the real line, the same step 2 already gives |F(ξ)² − M·S(ξ)| ≤ M − 1. If M were close to 1, F² would be
  uniformly close to sinc, which has simple zeros. Rouché turns this into a quantitative bound.
- The method's limit is about 1.153. That is the saddle level of |S(z)|·M/B_M(|Im z|) around z = 1; the binding
  saddle is near 1.23 ± 0.34i. Using the phase of e^{−2πizt}, and not only its modulus, improves B_M by about
  1% there, which moves the bound by about 0.002.
- `sanity_random.py` evaluates |ĥ(z)| / B_M on random step functions. The bound held in every test: the largest
  ratio was 0.595, against a required ≤ 1. Its output is in `sanity_output.txt`.

## Running the checks

    python3 -m pip install mpmath
    python3 check_c1d_lower.py 1.152 20000      # about 8 s

The best upper bound is 1.4546733458 (our `../c1d`). The constant therefore lies in [1.152, 1.4547].

Found and written up with Claude (Anthropic) agents under V. Patankar's direction. The interval check is
mechanical. The proof is short, but it has not yet been checked by a human expert.
