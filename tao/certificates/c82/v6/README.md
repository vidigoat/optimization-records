# C82 (Zhang–Zagier essential minimum): upper bound 0.2519694883

    C82  <=  B  <=  0.25196948827272  <  0.2519694883      (checkers/vcheck.py, Arb, 192 bits, quadrature tolerance 1e-11)

`vcheck.py` reports two parts:
- the Mahler part, as a two-sided Arb enclosure: [0.2519694078602313 ± 1.3e-26];
- the integral part, as an upper bound at quadrature tolerance 1e-11: 8.0412e-8.

The printed radius covers Arb's arithmetic only; see `../v3/README.md`. A second run at tolerance 1e-12 and 256 bits is
in `receipts/`.

Comparisons:
- 1.577e-7 below the candidate 0.2519696460\* of R. Pinho
  ([teorth/optimizationproblems#330](https://github.com/teorth/optimizationproblems/pull/330), 13 Nov), which was built
  on our `../v3`–`../v5`, and on which this certificate is built in turn;
- 1.947e-5 below our `../v5` (0.2519889614\*).

We replayed the #330 certificate (rpinho/math-certificates @ `977121e`, sha256 `0c4df95e…3077`) with our checker before
starting: B ≤ 0.2519696459872…, below its printed 0.2519696460 (`receipts/replay_pr330_977121e_vcheck.json`).

## Construction

Doche's limit-point construction in R. Pinho's certificate format, as in `../v3`–`../v5`.

* **Denominators:** the six monic degree-104 factors of the #330 certificate of 13 Nov, with re-optimised integer
  multiplicities (D = 1,000,000,000,064): `C82r8_b4ce16615ae6`, `C82r8_fcd3a4b2d94d`, `C82r8_1f8710f82b1d`,
  `C82r9_923bd2177a3d`, `C82r10_a6af275a730c`, `C82r10_76a523b0f880`. None is ours.
* **Numerators:** 37 factors, deg A / D = 0.89709.
  * Doche's P1, P2, P3, P4 and P6.
  * 15 from the #330 certificate (names `N82r*` kept).
  * 1 from our `../v4` (`n11_10_6`).
  * **16 new**, from our dual-weighted LLL numerator search today (`n71_*`, `n81_*`, `n91_*`, `n101_*`; degrees 24–42).
    One of them, `n101_2_0` (degree 42), gave 8.0e-8 on its own.

"New" means absent from the published #330 certificates and from all our earlier certificates; no novelty claim is
made beyond that. `search/provenance.py` reproduces this split by matching full coefficient arrays:

    python search/provenance.py certificate.json <#330 certificate.json> ../v5/certificate.json ../v4/certificate.json \
        ../v3/certificate.json ../v2/certificate.json ../certificate.json

## How it was found

The #330 certificate was unioned with our pooled search state behind `../v5` (158 denominators, 1440 numerators), so
the starting LP was #330's solution plus our library (0.2519696447 on the 2^17 grid). Our LLL numerator search of
`../v3` then ran on the warm joint LP, with 1792–1920 starts per round instead of 60 (`search/logs/U7`–`U11.log`). The
CVP denominator search of `../v3` (settings `DS_DEGS=104,88`, 100,000 enumeration nodes) ran for two rounds of 384
starts and found no denominator with negative reduced cost (best +0.004), so the gain is from numerators only. The search stopped when
four consecutive rounds found nothing. The certificate is a fine-grid (2^19) re-solve of state `U10`, with integer
multiplicities and exponents rounded down to dyadic rationals (`search/finalize2.py`).

The whole run took about three hours on one 64-core cloud VM.

## Verification

    python checkers/vcheck.py certificate.json --tol 1e-11 --prec 192 --target 0.2519694883 --out out.json

The checker of `../v2`–`../v5` is unchanged (same sha256). It checks:
- H1 to H3 exactly: integer coefficients, monic denominators, positive integer multiplicities, non-negative rational
  exponents (H1), deg A < D (H2), and exact gcds showing every numerator coprime to every denominator (H3);
- root enclosures for all factors, with the leading coefficients of the 7 non-monic numerators included;
- the Mahler part by Jensen's formula;
- the integral on an exact dyadic tiling.

The receipts are in `receipts/`.

## Credit

* Doche: the construction and the P-factors.
* R. Pinho: the certificate format, the joint LP, and the six denominators and fifteen numerators of his 13 Nov #330
  certificate used here.
* Grinsztajn and Épi/BakeLab: earlier factors.
* V. Patankar: the searches and this certificate.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed by the
checker above. It has not been verified line by line by a human expert.
