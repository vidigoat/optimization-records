# C86 ≤ 1.8055310732 (Schur–Siegel–Smyth trace constant, upper bound)

`certificate.json` (sha256 `09c05af07a4ef559983e2dc203d4c3e002c5b1f2248fde8b1c306fa2feb98d1f`) uses the certificate format
of R. Pinho ([teorth/optimizationproblems#298](https://github.com/teorth/optimizationproblems/pull/298)):
- exact rational edges and node densities (a piecewise-linear density, zero outside);
- auxiliary polynomials A with rational weights b;
- the bound T/m.

The bound:

    T/m = 195153337643980465525925899810927/108086391056860991292157728915456 = 1.8055310732071364...

This is 7.83e-5 below 1.80560940, Pinho's bound in #298.

## What it is

Smith's sufficiency criterion, in the Talebizadeh Sardari–Orloski construction, exactly as in #298:
- **Auxiliaries:** the same 28 auxiliary polynomials, of degree up to 21, from Pinho's certificate
  ([rpinho/math-certificates@0d82162](https://github.com/rpinho/math-certificates/tree/0d821621bab3d3d9408015f2679b89d5d5e9b379/c86/adaptive-cone)),
  plus 3 polynomials of degree 8, 12 and 15 found by dual-weighted lattice reduction (`search/extra_polys.json`). The 3 new ones end up with negligible weight (below 1e-11), so the auxiliary set is effectively Pinho's.
- **Density:** 5734 knots, with support [0.0215, 5.6185].

The gain over #298 comes from the density, in two steps:
1. **Re-optimizing Pinho's grid.** Each node of Pinho's grid becomes a triangular-density column, and the primal LP is solved to optimality. The LP requires f ≥ 2·10⁻⁶ at its sample points, with cutting planes as below. This alone gives 1.8055848345 (also certified: VALID, 0 failed cells).
2. **Refining the grid.** Half- and quarter-cell triangles are priced by their LP reduced cost, and the 300 most negative are added per round. After one round and re-converged cuts, the LP gives 1.8055310732.

The cutting planes run until f ≥ 10⁻⁶ at every local minimum of f. These minima are found on a 10⁶-point grid over [0.015, 5.65] and then located to ~10⁻¹² by repeated sub-gridding. The sub-gridding matters: a fixed grid misses the sharp V-shaped minima next to auxiliary roots.

## Checks

- `checkers/check86.py` gives a rigorous verdict, in Arb ball arithmetic at 128 bits. Its verdict is in `receipts/check86.txt`:
  - edges, densities and polynomials are exact; each Q is irreducible with positive leading coefficient, and b ≥ 0;
  - degree condition: Σ b deg Q < m, with ratio 0.898499;
  - ∫ log|Q| dμ > 0 for every Q (min 9.09·10⁻⁸);
  - f > 0 on every one of the 5732 support cells, via first-order Taylor models with rigorous remainders, bisected into 73,996 intervals. The minimum certified lower bound is 5.74·10⁻¹¹, and there are 0 failures;
  - the declared bound equals T/m exactly.

  Verdict: `VALID`, bound < 1.80560940.
- `checkers/fcheck86.py` is a float sanity check, not a proof: min f on 2·10⁶ points is 1.005·10⁻⁶.

Reproduce the check (python-flint ≥ 0.6):

    python3 checkers/check86.py certificate.json --procs 8 --beat 1.80560940

## Search code (`search/`)

- `denselp.py`: a dense primal-dual (Mehrotra) interior-point solver for this LP, using normal equations factored with the
  Accelerate/LAPACK Cholesky. It takes about 4 minutes per solve at 5,700 columns × 12,600 rows; HiGHS's IPM took about 25.
- `run_dense.py`: cutting planes with refined minima, triangle pricing, and auxiliary pricing by library reduced cost and dual-weighted LLL (`auxsearch.py`).
- `export86.py`: rounds to the exact certificate. Triangle weights become multiples of 10⁻¹⁵ and auxiliary weights are rounded down to multiples of 10⁻¹⁴.
- `receipts/search_log.txt`: the log of the run that produced this certificate. It was resumed from a converged state of Pinho's grid plus `extra_polys.json`.

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction (LP solver, search, certificate, checker and
text). The value is recomputed exactly by the checker above. It has not been verified line by line by a human expert.
