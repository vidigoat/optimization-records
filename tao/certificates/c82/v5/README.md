# C82 (Zhang–Zagier essential minimum): upper bound 0.2519889614

    C82  <=  B  <=  0.25198896138704  <  0.2519889614      (checkers/vcheck.py, Arb, 192 bits, quadrature tolerance 1e-11)

`vcheck.py` reports two parts:
- the Mahler part, as a two-sided Arb enclosure: [0.2519878666177742 ± 4.7e-26];
- the integral part, as an upper bound at quadrature tolerance 1e-11: 1.0948e-6.

The printed radius covers Arb's arithmetic only; see `../v3/README.md`.

Comparisons:
- 3.38e-6 below the candidate 0.2519923410* of R. Pinho ([teorth/optimizationproblems#330](https://github.com/teorth/optimizationproblems/pull/330)), which was built on our `../v3` and `../v4`;
- 1.90e-5 below our `../v4` (0.2520079724*).

## Construction

Doche's limit-point construction in R. Pinho's certificate format, as in `../v3` and `../v4`.

* **Denominators:** eight monic degree-104 factors (D = 999,999,999,960).
  * Ours from earlier versions: `D11_5_0` (`../v3`), `D11_12_0` and `D31_5_0` (`../v4`).
  * New in this search: `D51_16_0`.
  * Four from R. Pinho's #330 certificate (prefixed `Pin330c_` here): `C82r6_8c76cd0c4e53`, `C82r6_935548ec48fe`, `C82r6_083ca9b3950a` and `C82r7_dd0a7c419abf`.
* **Numerators:** 25 factors, with deg A / D = 0.89586.
  * Doche's P1, P2, P3, P4 and P6.
  * 20 from our dual-weighted LLL numerator search: 6 earlier, and 14 new in this search (`n51_*`).

## How it was found

The pooled search state behind `../v4` (all our factors, plus the factors of #330's earlier rounds) was unioned with
the #330 certificate's factors. The CVP denominator search and the LLL numerator search of `../v3` then alternated on one
warm joint LP for seven rounds, on a 64-core cloud VM. The certificate is a fine-grid (2^19) re-solve of that state,
with integer multiplicities and exponents rounded down to dyadic rationals.

## Verification

    python checkers/vcheck.py certificate.json --tol 1e-11 --prec 192 --target 0.2519889615 --out out.json

The checker of `../v2`–`../v4` is unchanged. It checks:
- H1 to H3 exactly;
- root enclosures for all factors;
- the Mahler part by Jensen's formula;
- the integral on an exact dyadic tiling.

The receipt is in `receipts/`.

## Credit

* Doche: the construction and the P-factors.
* R. Pinho: the certificate format, the joint LP, and the four #330 denominators above.
* Grinsztajn and Épi/BakeLab: earlier factors.
* V. Patankar: the searches and this certificate.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed by the
checker above. It has not been verified line by line by a human expert.
