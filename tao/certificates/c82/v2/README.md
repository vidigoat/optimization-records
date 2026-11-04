# C82 ≤ 0.2524177045 (essential minimum of the Zhang–Zagier height, upper bound)

`certificate.json` is in R. Pinho's format (github.com/rpinho/c82-height-certificate). Same Doche limit-point
construction and same proof as Pinho's `degree44-denominators-v2` (open PR teorth/optimizationproblems#301,
0.253091834), with a new denominator library.

    B = (1/D)[ Σ m_j log M(R_j∘T) + ∫₀¹ max(0, Σ q_i log|S_i(T(e^{2πit}))| − Σ m_j log|R_j(T(e^{2πit}))|) dt ],  T(z) = z(1−z)

**Certified value: B ≤ 0.25241770445** (Pinho v2: 0.25309183355; improvement 6.74·10⁻⁴).
Mahler part 0.2524164715, integral part 1.233·10⁻⁶. SHA-256 of `certificate.json` in `SHA256SUMS`.

## What changed relative to Pinho v2

* **Denominators: all four are new.** Monic, irreducible, degree 104 (`D101_11_0`, `D111_0_0`, `D111_3_0`,
  `D111_3_1`), D = 999,999,999,960. None of Pinho's P9/MD* factors keeps positive weight. Normalised log-Mahler
  measures of the new factors are 0.252405–0.252442 per degree, against 0.25272 for P9 and 0.25308–0.25323 for
  the MD* factors. All four lie in one residue class mod X(X−1)P4: the class of a degree-104 product of a
  degree-44 MD-class factor and a degree-60 factor.
* **Numerators:** 24 with positive exponent (deg A / D = 0.9034). Doche's P1, P2, P3, P4 and P6; Pinho's MS192 and
  MS360; `cg2_4_0` from our 27 Oct certificate; and 16 new (`n*`, degrees 10–42, a few with leading coefficient
  2, 6 or 9).
* Multiplicities and exponents were re-optimised jointly by LP on a 2¹⁹-point grid, then rounded exactly
  (multiplicities ≈ 10¹²·μ; exponents in 2⁻²⁰ℤ, scaled down by 10⁻⁹ so that deg A < D strictly).

## How the denominators were found (`search/`)

This uses Pinho's reduced cost for a monic denominator, rc(R) = log M(R∘T) − Σ ν_k log|R(T(e^{2πit_k}))| − λ deg R,
and his residue-class restriction R = R₀ + X(X−1)P4·U, which keeps R(0), R(1) and R mod P4 fixed and avoids
collapsing onto the essential numerators X, X−1 and P4. Our implementation (`dsearch.py`, function `job2`)
differs in the lattice step:

1. The lattice of H·Xⁱ evaluations is weighted by the (uniform − ν) strata masses divided by |R_cur|², then
   LLL-reduced on its own with python-flint. The closest points to −R_cur come from **Schnorr–Euchner
   enumeration** (CVP, best 16 within a node budget), not from Kannan embedding. With embedding the IRLS step
   almost never moved, because the neighbours it offered cost +0.3 in rc. Enumeration finds lattice points
   within about 0.01 in rc, and the iteration then descends.
2. Starts are near-zero-rc library factors, or products of library factors of the target degree. Each start
   gets an initial weight profile from one of: randomly perturbed roots, noisy log-weights, or a power of the
   current weights. Some starts also take a unit-class shift mod P4 (Pinho's idea). IRLS iterations with exact
   1/|R_cur|² weights follow.
3. Candidates are scored by rc from FLINT root enclosures (in float), factored exactly, and added as LP
   columns.

The search alternated this with our earlier numerator search (`cg.py`, iterated dual-weighted LLL) on a
warm-started HiGHS LP (`alt.py`). It began at degree 44, as Pinho did. **Product starts at degree 60, 88 (44+44)
and 104 (44+60) turned out to be much richer:** right after the dual moves, up to several percent of starts give
a negative-cost factor. The bound went 0.2530918 → 0.25274 (degrees 44 and 60) → 0.25265 (88) → 0.25245 (104),
and reached 0.252418 at degree 104 with a 540-point quadrature. At degree 120 and above (composed degree ≥ 240)
the quadrature is too coarse, and nothing was found.

## Verification

| check | certified upper bound |
|---|---|
| `checkers/vcheck.py` (our Arb checker), tol 1e-11, 192-bit | 0.25241770444734 |
| Pinho's `verify_exact.py` (unmodified), 256-bit, tol 1e-11 | 0.25241770444732 |
| Pinho's `review/integrate_gauss.py` (unmodified), 256-bit, tol 1e-14, whole circle | 0.25241770444731 |
| Pinho's `review/audit_roots.py` (unmodified) | passed; normalised Mahler part 0.2524164715294 |
| Pinho's `review/audit_algebra.py` (only the pinned hash, D and deg A replaced) | passed; all 96 numerator/denominator pairs coprime by modular witness |

On Pinho's v2 certificate, `vcheck.py` gives 0.25309183355156; his `verify_exact.py` gives 0.25309183355151.

Negative controls: `vcheck.py` rejects four corrupted certificates (a non-monic denominator, a numerator
sharing a denominator factor, deg A ≥ D, a zero multiplicity; see `receipts/vcheck_negative_controls.txt`).
Pinho's `check_verifier.py` controls pass in this environment (python-flint 0.9.0).

To reproduce (python-flint 0.9.0):

    python checkers/vcheck.py certificate.json --tol 1e-11 --prec 192 --target 0.2524178
    python checkers/pinho_verify_exact.py certificate.json --target 0.2524178 --precision 256 --tol 1e-11 --out out.json

## Notes

* Search-grid values, reduced costs and LP values are heuristic. Only the checkers' outputs are claims.
* The search was limited by compute, not saturated. The 64-core lane was still finding degree-88 and degree-104
  factors worth 10⁻⁵ each when it was stopped. The last 15-core lane found two factors worth 2·10⁻⁵ in its
  final hour.
* Credit:
  * Doche: limit-point formula and the P-factors.
  * Grinsztajn: perturbative denominator factors.
  * Épi: Ed5 and unequal multiplicities.
  * R. Pinho: certificate format, joint LP, residue-class complementary-measure denominator search, and the
    MD/MS libraries.
  * vidigoat: the iterated numerator LLL search, this CVP denominator search, and this certificate.
