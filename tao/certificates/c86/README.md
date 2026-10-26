# C86 ≤ 1.808109995748 (Schur–Siegel–Smyth trace constant, upper bound)

Certificate `c86_deg8_1.808109995748.json` in the format of R. Pinho / D. Koshta (edges, densities, A, b, bound):
a step density on [11/200, 519/100] with cells of width 1/200 and 21 auxiliary polynomials, bound

    T/m = 144648799659767282439 / 79999999999999829200 = 1.80810999574...

Previous best 1.8084990 (Koshta, all 883 polynomials of degree ≤ 7). The gain comes from degree 8: `enum.py` lists
every irreducible monic integer polynomial with all roots in [0, 5.3] and trace ≤ 2·degree (it reproduces the
883 of degree ≤ 7 exactly: 3, 3, 7, 17, 51, 167, 635, and finds 2954 of degree 8). `lp.py` solves the discretised
Talebizadeh Sardari–Orloski primal LP over them with cutting planes (f is checked on a 200-points-per-cell grid
and violated points are added until none remain), then rounds to exact rationals.

Checked with Koshta's independent rigorous verifier `checkers/c86_trace.py`
(github.com/DennisKoshta/c86-sss-trace-certificate): `VALID 1.808109995748`, `BEATS 1.8084990`.

## Update: C86 ≤ 1.807768536863 (finer cells)

`c86_deg8_h400_1.807768536863.json`: the same LP with cells of width 1/400 (2054 cells), the same degree ≤ 8 library
(pruned to the 300 best by reduced cost on a coarse solve), cutting planes until f ≥ 10⁻⁶ on the fine grid. 19 polynomials
receive positive weight. Bound T/m = 578485931795944337659/319999999999999303200 = 1.80776853686...; Koshta's verifier
reports `VALID 1.807768536863`.

    python3 lp.py deg1to7_and_deg8.txt out.json --h 400 --sub 6 --margin 1e-6 --prune 300
