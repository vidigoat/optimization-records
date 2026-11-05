# C82 (Zhang–Zagier essential minimum): upper bound 0.2520591386

    C82  <=  B  <  0.2520591386,     B in [0.2520591385260044869143520 +/- 1.67e-26]   (checkers/vcheck.py, Arb, 192 bits)

That is 3.46e-4 below the candidate 0.252404806* of R. Pinho
([teorth/optimizationproblems#330](https://github.com/teorth/optimizationproblems/pull/330)), and 3.59e-4 below our
previous 0.2524177045* (`../v2`).

## Construction

Doche's limit-point construction (C. Doche, J. Théor. Nombres Bordeaux 13 (2001), 103–110), in R. Pinho's certificate
format:

    B = (1/D) [ sum_j m_j log M(R_j o T) + int_0^1 max(0, log|A(T(u))| - log|Q(T(u))|) dt ],   T(z) = z(1-z),

with Q = prod R_j^{m_j} (monic R_j), A = prod S_i^{q_i}, deg A < D = deg Q, and gcd(S_i, R_j) = 1 for q_i > 0.

* **Denominators:** 3 new monic degree-104 factors `D11_5_0`, `D11_5_1`, `D11_5_2`, with multiplicities
  5182713980, 3042236791, 1390433844 (D = 999,999,999,960). Their normalised Mahler measures log M(R∘T)/208 are
  0.251994, 0.252086 and 0.252106. All three are below every factor in the earlier certificates.
* **Numerators:** 21 factors, with deg A / D = 0.90528:
  * Doche's P1, P2, P3, P4, P6;
  * R. Pinho's MS360 (#301);
  * MS204_c1 and cg2_4_0_c2;
  * 13 found by our dual-weighted LLL numerator search.
* Split: Mahler part 0.2520394006582630, integral part ≤ 1.97379e-5.

## How it was found

We continued the CVP-enumeration denominator search of `../v2` (weighted-evaluation lattice, Schnorr–Euchner
enumeration, product starts). It alternated with the LLL numerator search on a warm-started joint LP. The starting
library was the union of the `../v2` search state and the six denominators of #330. None of the #330 factors appears in
the final certificate. After two rounds that added single factors, one round produced degree-104 factors with
log M/deg about 0.2521. Each later round then lowered the LP value by 5·10⁻⁵ to 1·10⁻⁴. The certificate is a fine-grid
(2^19) re-solve of the round-6 state, with multiplicities rounded to integers and numerator exponents rounded down
to dyadic rationals. The search was still improving when this snapshot was taken.

## Verification

    python checkers/vcheck.py certificate.json --tol 1e-11 --prec 192 --target 0.2524 --out out.json

This takes 8 minutes on 6 cores. Receipt: `receipts/vcheck_tol1e-11_prec192.json`.
* **Exact checks:** H1 (integer coefficients, monic denominators, positive integer multiplicities, q_i ≥ 0), H2 (deg A < D) and H3 (exact gcds).
* **Root enclosures:** Arb encloses every root of each composed factor, with multiplicities summing to the degree.
* **Mahler part:** computed by Jensen's formula from the enclosures.
* **Integral part:** an upper bound on an exact dyadic tiling of [0, 1/2], with rigorous remainders (405k panels).

The float value on a 2^20 grid is 0.2520591385373, consistent with the enclosure. `vcheck.py` is the checker of
`../v2`, unchanged. Arhan Canli's independent checker replayed `../v2`, the first version of this certificate, and #301.

## Credit

* Doche: the limit-point construction and the P-factors.
* R. Pinho: the certificate format, the joint LP, the residue-class denominator search, the degree-44 libraries
  (#229, #301), and the six denominators of #330 that were in the starting library.
* Grinsztajn and Épi/BakeLab: earlier factors.
* V. Patankar: the CVP denominator search, the LLL numerator search and this certificate.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction; the value is recomputed by the
checker above. It has not been verified line by line by a human expert.
