# Improved lower bounds for 16a, 20a, 20c and 70a

Four lower bounds, each from one explicit object, with a script that recomputes it:
https://github.com/vidigoat/optimization-records (see `tao/PROOFS.md`).

| const | side | was | now | object |
|---|---|---|---|---|
| 16a | lower | (β+π)/π ≈ 1.82283 | 2.156210 (rigorous) | truncated BDZ conditional extremal on [0,2π]², 0 < \|k\|² ≤ 1200², μ = 42/5, evaluated in interval arithmetic |
| 20a | lower | 2 | 8 | n = 1, X = E − 1 with E standard exponential: Var(X²) = 9 − 1 |
| 20c | lower | √(π/2) | √3 | n = 1, uniform on [−√3, √3], A = (−∞, 0] |
| 70a | lower | 1 | √(1/2 + 1/√π) ≈ 1.031596 | n = 2, square of side √π and unit disk; Steiner's formula and minimal parallelogram perimeter |

**16a** turns the 2013 source's numerical value (≈ 2.15627) into a certified bound: the sums S₀, S₁,
S₂ over the lattice and the logarithms are evaluated with `mpmath.iv` at 40 digits, and the reported
value is the lower end of the interval. As the cutoff grows the value approaches 2.15627, which
confirms the torus normalisation matches the source.

**20a, 20c, 70a** replace a classical example or trivial bound with a better explicit one; the full
arguments are short and in the Comments column of each new row.

Each page gets a new row (earlier rows kept) and a reference entry; the README table and "Recent
progress" are updated. Inline math follows the escaping rules in CONTRIBUTING (`\_`, `\lvert`), and
each new row has three columns.
