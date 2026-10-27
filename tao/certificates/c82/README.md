# C82 ≤ 0.2534350402 (essential minimum of the Zhang–Zagier height, upper bound)

`certificate.json` is in R. Pinho's format (github.com/rpinho/c82-height-certificate): Doche's limit-point
construction with 7 denominator factors (Doche, Grinsztajn, Épi; multiplicities re-optimised, D = 1,000,000,000,004)
and 35 numerator factors with rational exponents (deg A / D = 0.9417). 18 of Pinho's numerators keep positive weight;
**17 are new** (`cg*`, degrees 17–40).

    B = (1/D)[ Σ m_j log M(R_j∘T) + ∫₀¹ max(0, Σ q_i log|S_i(T(e^{2πit}))| − Σ m_j log|R_j(T(e^{2πit}))|) dt ],  T(z) = z(1−z)

**Certified value** B ≤ 0.25343504013 (previous record 0.2534446873, Pinho).

| check | certified upper bound |
|---|---|
| `vcheck.py` (own Arb checker), tol 1e-10 | 0.2534350401297 |
| `vcheck.py`, tol 1e-11, 192-bit | 0.2534350401282 |
| Pinho's `verify.py`, tol 1e-9 | 0.2534350401467 |
| Pinho's `independent_check.py`, tol 1e-8 | 0.2534350401334 |

`vcheck.py` reproduces Pinho's certified value 0.2534446873 on his own certificate, and rejects corrupted ones.
It composes each factor exactly with T in Z[z], encloses the full root multiset in Arb, computes the Mahler part by
Jensen, and bounds the integral over an exact dyadic tiling (whole-panel bounds, second-order Taylor bounds, and a
Lipschitz bound on panels where the integrand crosses zero). It also checks integrality, monic denominators,
positive multiplicities, q ≥ 0, deg A < D and coprimality.

**How the new factors were found** (`search/`). A plain lattice that makes |S|² small at the LP dual-measure points
only returns multiples of known factors. Instead, an iterated LLL step targets the dual-weighted average of log|S|:
start from a product S₀ of library factors, weight each dual-support point by its mass / |S₀|², LLL-reduce, keep the
best vector coprime to all denominators (a smaller weighted |S|² implies a smaller weighted log-average), factor it,
screen by exact reduced cost, and add it as an LP column (warm-started HiGHS). Each accepted factor gained
10⁻⁸–10⁻⁶; the search saturated after 17 factors.

    python3 vcheck.py certificate.json out.json 1e-10
