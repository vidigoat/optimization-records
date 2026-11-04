# C16 = 2.15622552815421… (Brezis–Gallouët–Wainger remainder constant on T², determined to 14 digits)

**2.1562255281542128 ≤ C16 ≤ 2.15622552815422**

Previous bounds: [J2026] proved 2.1562255281542061 ≤ C16 (that bound is what the lower side was printed from), and [maximusrufus] proved C16 ≤ 2.1562256466.
- **Upper side:** improved by 1.2·10⁻⁷. It now meets the lower bound, leaving a window of 7·10⁻¹⁵.
- **Lower side:** our evaluation reproduces [J2026]. It is 7·10⁻¹⁵ higher only because μ is tuned more finely.

## Upper bound (`verify.py`, `certificate.json`)

The method is the same Cauchy–Schwarz dual bound as [maximusrufus]. For every μ > 0 and every zero-mean u ∈ H²(T²), T² = [0, 2π]²:

    4π ‖u‖∞² / ‖∇u‖² ≤ G(μ)(1 + μ δ(u)) / π,      G(μ) = Σ_{k≠0} 1 / (|k|²(1 + μ|k|²)).

Translation invariance lets us put the maximum of |u| at 0, so it is enough to bound |u(0)|.

U(δ, μ) = G(μ)(1+μδ)/π − ln δ − ln(1+ln δ) is convex in δ ≥ 1, so on each interval it is at most its value at one of the two endpoints.

The bound has two parts:
1. **δ ∈ [1, 2¹⁸].** `certificate.json` lists 353 intervals that cover [1, 2¹⁸] with no gaps, each with its own μ. All of them are exact binary floats.
2. **δ ≥ 2¹⁸.** Take μ = 1/(δ ln δ) and count lattice points: #{k ≠ 0 : |k| ≤ R} ≤ π(R + 1/√2)². Together with G(2⁻²⁰) this gives

       G(μ) ≤ G(μ₀) + π[ln(μ₀/μ) + π c/a₀ + c²/A₀]

   for μ ≤ μ₀ = 2⁻²⁰. The resulting bound on the tail is 2.0934.

**How G(μ) is enclosed** (`gsum.py`):
- G(μ) is written as Σ_{k≠0} [1/r − 1/(r + 1/μ)].
- Each row is summed exactly with Σ_y 1/(y² + s²) = (π/s) coth(πs), so nothing is truncated.
- Rows x ≥ N use the convex tail φ(x) = 1/x − 1/√(x² + 1/μ). It is bracketed between a trapezoid and a midpoint sum, which have closed-form integrals.

Everything is computed in Arb ball arithmetic (python-flint), at 160 bits.

**Verifier output** (`verify_output.txt`):

    [1, 2^18]: 353 intervals, 353 values of mu; max endpoint bound [2.1562255281542197 +/- 4.95e-17]
    delta >= 2^18: bound [2.09340515153 +/- 4.25e-12]
    CERTIFIED: C16 <= 2.15622552815422

**Sanity checks** (not part of the proof):
- A brute-force lattice sum over |k| ≤ 3000 agrees with the enclosed G(μ) to a relative error of 9·10⁻⁸, the accuracy of the brute-force sum.
- Rerunning the boxes near the peak with N = 40000 at 256 bits gives the same bound.

## Lower bound (`lower.py`)

The trial function is the untruncated conditional extremal of [BDZ2013], as used by [J2026]: c_k = 1/(|k|²(1 + μ|k|²)) with μ = 0.122110476195 (`0x1.f42a1d5c1aabdp-4`). Its sums are enclosed with the same row-sum method:
- S₀ = G
- S₁ = G − H/μ
- S₂ = H/μ², where H = Σ_{k≠0} (|k|² + 1/μ)⁻²

The result is δ = 3.928836059681948…, and

    L ≥ [2.156225528154213 +/- 1.68e-16]   (lower endpoint 2.15622552815421283…)

## Why the two sides meet

For a fixed δ, maximising |u(0)|² subject to ‖∇u‖² = 1 and ‖Δu‖² = δ is a quadratic problem with two constraints. Its maximiser has c_k ∝ 1/(|k|²(λ + ν|k|²)), the conditional-extremal family. The Cauchy–Schwarz weights are exactly the dual of that problem, so for each δ the upper bound equals the lower bound at the matching μ. That leaves only the numerical work: enough intervals near δ* ≈ 3.9288 and exact lattice sums.

## Reproduce

    python verify.py certificate.json --crosscheck 20     # about 15 s; prints CERTIFIED
    python lower.py 0.122110476195
    python build.py 2.15622552815422 certificate.json      # regenerates the intervals (about 15 s)

## Credit

- [BDZ2013] Bartuccelli, Deane, Zelik: the problem and the conditional extremal.
- [maximusrufus] (github.com/maximusrufus/c16-bgw-certificates): the Cauchy–Schwarz dual bound with these weights and the δ-convexity argument.
- [J2026] (github.com/Julsgaard/c16a-certificate): the untruncated family with tuned μ.

The search, certificate, checkers and this text were prepared with Claude (Anthropic) agents under V. Patankar's direction.
