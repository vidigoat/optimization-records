# Tao's table of optimization constants

Five improved bounds for entries of Terence Tao's
[optimization constants repository](https://github.com/teorth/optimizationproblems).

| Entry | Constant | Side | Listed | New |
|---|---|---|---|---|
| [16a](https://teorth.github.io/optimizationproblems/constants/16a.html) | Brezis–Gallouët–Wainger remainder constant on the 2D torus | lower | 1.82283 | **2.156224** |
| [70a](https://teorth.github.io/optimizationproblems/constants/70a.html) | Reverse Brunn–Minkowski constant | lower | 1 | **1.201927** |
| [23b](https://teorth.github.io/optimizationproblems/constants/23b.html) | Decay exponent of the minimal condition number of ±1 matrices | upper | 1 | **1/2** |
| [20a](https://teorth.github.io/optimizationproblems/constants/20a.html) | Thin-shell constant (variance form) | lower | 2 | **8** |
| [20c](https://teorth.github.io/optimizationproblems/constants/20c.html) | KLS constant | lower | √(π/2) ≈ 1.2533 | **√3 ≈ 1.7321** |

The arguments are in [PROOFS.md](PROOFS.md). Each numerical value is recomputed by a script in
[`certificates/`](certificates), in interval arithmetic where rounding could matter:

| Script | Checks | Time |
|---|---|---|
| `certificates/bgw.py --K 4000 --mu 41/5` | 16a: sums over ~50 million lattice points, interval arithmetic | about a minute |
| `certificates/rbm_cube.py --n 300000 --t 0.013` | 70a: the cube-and-ball bound in dimension 300,000 | seconds |
| `certificates/rbm.py` | 70a in dimension 2 (square and disk), plus 200,000 random linear maps | seconds |
| `certificates/elementary.py` | 20a and 20c in exact arithmetic | instant |

23b is a short proof and needs no computation; PROOFS.md also lists a numerical sanity check.
