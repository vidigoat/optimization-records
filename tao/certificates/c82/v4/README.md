# C82 (Zhang–Zagier essential minimum): upper bound 0.2520079724

    C82  <=  B  <=  0.25200797232846  <  0.2520079724      (checkers/vcheck.py, Arb, 192 bits, quadrature tolerance 1e-11)

`vcheck.py` reports two parts, which add to the upper bound above:
- the Mahler part, as a two-sided Arb enclosure: [0.2520059412697726 ± 1.6e-26];
- the integral part, as an upper bound at quadrature tolerance 1e-11: 2.0311e-6.

The printed ball radius covers only Arb's arithmetic, not the quadrature tolerance (see `../v3/README.md`).

Comparisons:
- 9.95e-6 below the candidate 0.2520179248* of R. Pinho ([teorth/optimizationproblems#330](https://github.com/teorth/optimizationproblems/pull/330), 5 Nov update), which was built on our `../v3`;
- 5.12e-5 below our previous 0.2520591386* (`../v3`).

## Construction

Doche's limit-point construction, in R. Pinho's certificate format, exactly as in `../v3`:

    B = (1/D) [ sum_j m_j log M(R_j o T) + int_0^1 max(0, log|A(T(u))| - log|Q(T(u))|) dt ],   T(z) = z(1-z).

* **Denominators:** five monic degree-104 factors (D = 999,999,999,960).
  * `D11_5_0` comes from `../v3`.
  * `D11_9_0`, `D11_12_0`, `D11_12_1` and `D31_5_0` are new, found by our CVP-enumeration denominator search.

  Their log M(R∘T)/104 values are 0.251994, 0.252021, 0.252002, 0.252015 and 0.252002.

  | factor | multiplicity |
  |---|---|
  | `D11_5_0` | 1362219132 |
  | `D11_9_0` | 1354335534 |
  | `D11_12_0` | 2275789042 |
  | `D11_12_1` | 1761580524 |
  | `D31_5_0` | 2861460383 |
* **Numerators:** 27 factors, with deg A / D = 0.90185:
  * Doche's P1, P2, P3, P4 and P6;
  * R. Pinho's MS360 (#301);
  * 21 found by our dual-weighted LLL numerator search (`n*`).
* No factor of #330's 5 Nov update is used.

## How it was found

The CVP denominator search and the LLL numerator search of `../v3` continued to alternate on one warm-started joint LP.
- The starting library was the union of two `../v3` search lanes: 134 denominators and 1259 numerators.
- The certificate is a fine-grid (2^19) re-solve of the saved state, with multiplicities rounded to integers and
  exponents rounded down to dyadic rationals (`finalize2.py` of `../v3`).

## Verification

    python checkers/vcheck.py certificate.json --tol 1e-11 --prec 192 --target 0.2520079724 --out out.json

This takes 7 minutes on 6 cores. The receipt is `receipts/vcheck_tol1e-11_prec192.json`. `vcheck.py` is the checker of
`../v2` and `../v3`, unchanged. It checks:
- H1 to H3 exactly;
- the root enclosures of every composed factor;
- the Mahler part by Jensen's formula;
- the integral on an exact dyadic tiling with rigorous remainders.

The float value on a 2^20 grid is 0.2520079723536.

## Credit

* Doche: the limit-point construction and the P-factors.
* R. Pinho: the certificate format, the joint LP and MS360 (#301).
* Grinsztajn and Épi/BakeLab: earlier factors.
* V. Patankar: the CVP denominator search, the LLL numerator search and this certificate.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed by the
checker above. It has not been verified line by line by a human expert.
