# C86 ≤ 1.8054956 (Schur–Siegel–Smyth trace constant, upper bound)

`certificate.json` (sha256 `182e76bfa2289c31ba262fbd4e515eee734ebd0c78a5c23f40f70db18072d397`) is in R. Pinho's format
([#298](https://github.com/teorth/optimizationproblems/pull/298)). Its bound:

    T/m = 780597977251274792019360777060007/432345564227458232572250425393152 = 1.80549551525085145...

This is 1.14e-4 below 1.80560940 (#298) and 3.56e-5 below our `../c86_v2` (1.8055311).

## What changed from `../c86_v2`

The same search was continued for two more rounds. Each round prices the 300 most negative half- and quarter-cell
triangles by LP reduced cost, adds them, re-solves, and re-converges the cutting planes. The density now has 6999 knots
(5734 in v2) and the same 28 auxiliaries of #298, plus 5 found by lattice reduction, all with negligible weight.

There are two bookkeeping fixes. Identical auxiliary columns found twice are merged in `export86.py` (their weights add),
and the search de-duplicates new auxiliaries. Rounds 3 and later ran on a cloud VM from the round-2 state
(`receipts/search_log_part1.txt` is the laptop part, `receipts/search_log_part2.txt` the rest).

## Checks

`checkers/check86.py` (Arb, 128 bits; receipt `receipts/check86.txt`) verifies:
- the exact structure;
- the degree condition (Σ b deg Q ≈ 0.8988 m < m);
- every ∫ log|Q| dμ > 0 (min 9.06·10⁻⁸);
- f > 0 on all 6997 support cells, using 96,103 intervals, with minimum certified lower bound 2.1·10⁻¹¹ and 0 failures;
- the bound equals T/m exactly.

Verdict: `VALID`.

    python3 checkers/check86.py certificate.json --procs 8 --beat 1.80560940

Found and prepared with Claude (Anthropic) agents under V. Patankar's direction. The value is recomputed exactly by the
checker above. It has not been verified line by line by a human expert.
