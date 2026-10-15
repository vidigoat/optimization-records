Four improved lower bounds and one improved upper bound, each from one explicit object or a short argument, with a script that recomputes it:
https://github.com/vidigoat/optimization-records/tree/main/tao (proofs in `tao/PROOFS.md`, scripts in `tao/certificates/`).

| const | side | was | now | object |
|---|---|---|---|---|
| 16a | lower | (β+π)/π ≈ 1.82283 | 2.156224 (rigorous) | truncated BDZ conditional extremal on [0,2π]², 0 < \|k\|² ≤ 4000², μ = 41/5, evaluated in interval arithmetic |
| 20a | lower | 2 | 8 | n = 1, X = E − 1 with E standard exponential: Var(X²) = 9 − 1 |
| 20c | lower | √(π/2) | √3 | n = 1, uniform on [−√3, √3], A = (−∞, 0] |
| 70a | lower | 1 | 1.201927 (rigorous) | cube and ball in n = 300,000; Szász's inequality shows the cube is already in its worst position, then Steiner's formula |
| 23b | upper | 1 | 1/2 | odd n: off-diagonal Gram entries are odd; Bhatia–Davis gives κ(n) ≥ 1 + (2 − √2)/√n |

**16a** turns the 2013 source's numerical value (≈ 2.15627) into a certified bound: the sums S₀, S₁,
S₂ over the lattice and the logarithms are evaluated with `mpmath.iv` at 40 digits, and the reported
value is the lower end of the interval. As the cutoff grows the value approaches 2.15627, which
confirms the torus normalisation matches the source.

**70a** uses a cube and a ball: Szász's inequality on the Gram matrix of u(cube) shows every intrinsic volume of
u(cube) is at least the cube's, so no linear map helps; the resulting bound increases with the dimension
(1.0316 for n = 2) and is evaluated in interval arithmetic at n = 300,000.

**20a, 20c** replace the Gaussian with a better explicit one-dimensional example; the arguments are in the
Comments column of each new row.

**23b** is an elementary argument for odd n (all off-diagonal Gram entries are odd, then the Bhatia–Davis
inequality); it shows κ(n) − 1 ≥ (2 − √2)/√n for every odd n ≥ 3, so the decay exponent is at most 1/2.

Each page gets a new row (earlier rows kept) and a reference entry; the README table and "Recent
progress" are updated. Inline math follows the escaping rules in CONTRIBUTING (`\_`, `\lvert`), and
each new row has three columns.

Found with search and proof assistance from Claude (Anthropic); every value is recomputed by the linked scripts
(`./check.sh` in the linked repository re-checks all of them in a few seconds).
