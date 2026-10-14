# Proofs of the four bounds

Each bound comes from one explicit object. Each script in this folder recomputes its value; the
arguments below show why the value is a valid bound. Definitions follow the pages of Terence Tao's
[optimization constants repository](https://github.com/teorth/optimizationproblems).

## 16a. Brezis–Gallouët–Wainger remainder constant: L ≥ 2.156210

L is defined as max over δ ≥ 1 of 4πΘ(δ) − ln δ − ln(1 + ln δ), where Θ(δ) is the supremum of
|u(0)|² over zero-mean u on the torus with ‖∇u‖² = 1 and ‖Δu‖² = δ. Any single admissible u therefore
gives L ≥ 4π|u(0)|²/‖∇u‖² − ln δ(u) − ln(1 + ln δ(u)), provided δ(u) ≥ 1.

On T² = [0, 2π]² take u(x) = Σ c_k e^{ik·x} over the lattice points 0 < |k|² ≤ 1200², with
c_k = 1 / (|k|²(|k|² + μ)) and μ = 42/5. This is the conditional extremal of Bartuccelli, Deane and
Zelik (2013), truncated. u is real (c_k = c_{−k}) and has zero mean. By Parseval,
‖∇u‖² = (2π)² S₁ and ‖Δu‖² = (2π)² S₂, and u(0) = S₀, where S₀ = Σ c_k, S₁ = Σ |k|² c_k², S₂ = Σ |k|⁴ c_k².
So δ = S₂/S₁ ≈ 3.98730 ≥ 1, and

    L ≥ S₀² / (π S₁) − ln δ − ln(1 + ln δ) = 2.15621085…

`bgw.py` computes the three sums over n = |k|², weighted by the number of lattice points with
|k|² = n, in interval arithmetic (mpmath.iv, 40 digits), so the printed lower end is rigorous. The
interval has width about 10⁻³⁴. The torus normalisation agrees with the source: as the cutoff grows,
the value approaches the 2013 paper's numerical value 2.15627 (cutoff 50: 2.15174; 200: 2.15592;
800: 2.15621).

## 20a. Thin-shell constant (variance form): C ≥ 8

C is the least constant with Var(|X|²) ≤ C n for every isotropic log-concave X in Rⁿ. Take n = 1 and
X = E − 1, with E a standard exponential variable. Its density e^{−(x+1)} on (−1, ∞) is log-concave,
E X = 0 and E X² = 1, so X is isotropic. The k-th central moment of E is Σ_j C(k, j)(−1)^{k−j} j!,
which is 0, 1, 2, 9 for k = 1, 2, 3, 4. So Var(X²) = E X⁴ − (E X²)² = 9 − 1 = 8 and C ≥ 8. A product of
n independent copies is isotropic and log-concave with Var(|X|²) = 8n, so the ratio is 8 in every
dimension.

`elementary.py` computes the moments exactly.

## 20c. KLS constant: C ≥ √3

ψ_μ = 1/h_μ, where h_μ = inf over A of μ⁺(A) / min(μ(A), 1 − μ(A)). For any single set A,
ψ_μ ≥ min(μ(A), 1 − μ(A)) / μ⁺(A). Take n = 1 and μ uniform on [−√3, √3]: mean 0, variance 3/3 = 1,
log-concave. With A = (−∞, 0], μ(A) = 1/2 and μ⁺(A) is the density at 0, 1/(2√3). So ψ_μ ≥ √3 and
C ≥ √3 ≈ 1.732051, against the listed √(π/2) ≈ 1.253314.

## 70a. Reverse Brunn–Minkowski constant: C ≥ √(1/2 + 1/√π) ≈ 1.031596

C is the least constant such that every pair of symmetric convex bodies has some u with |det u| = 1
and |u(B₁) + B₂|^{1/n} ≤ C(|B₁|^{1/n} + |B₂|^{1/n}). Take n = 2, B₁ a square of side a = √π, and B₂
the unit disk. For any u with |det u| = 1, u(B₁) is a parallelogram P of area π. If its sides are p
and q with angle t, then π = pq sin t ≤ pq, so its perimeter 2(p + q) ≥ 4√(pq) ≥ 4√π. By Steiner's
formula, |P + B₂| = |P| + per(P)·1 + π ≥ 2π + 4√π for every u. Since |B₁|^{1/2} + |B₂|^{1/2} = 2√π,
every u gives a ratio of at least √(2π + 4√π) / (2√π) = √(1/2 + 1/√π).

`rbm.py` evaluates the constant in interval arithmetic and checks the inequality on 200,000 random
maps with determinant 1 (the smallest ratio found is 1.031603, at the square position).
