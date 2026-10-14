# Optimization records

New best-known bounds and world-record packings, each with a script that checks it.

**Vidit Patankar**, October 2026.

## Results

### Terence Tao's table of optimization constants

Four improved lower bounds for entries of the
[optimization constants repository](https://github.com/teorth/optimizationproblems).

| Entry | Constant | Listed lower bound | New lower bound | Object | Check |
|---|---|---|---|---|---|
| 16a | Brezis–Gallouët–Wainger remainder constant on the 2D torus | 1.82283 | **2.156224** (rigorous) | a trigonometric polynomial with about 50 million terms | `python3 tao/bgw.py --K 4000 --mu 41/5` |
| 20a | Thin-shell constant (variance form) | 2 | **8** | a centred exponential random variable | `python3 tao/elementary.py` |
| 20c | KLS constant | √(π/2) ≈ 1.2533 | **√3 ≈ 1.7321** | the uniform distribution on [−√3, √3], cut at 0 | `python3 tao/elementary.py` |
| 70a | Reverse Brunn–Minkowski constant | 1 | **1.20192** (rigorous) | a cube and a ball in 300,000 dimensions | `python3 tao/rbm_cube.py --n 300000 --t 0.013` |

16a turns the value the 2013 source computed numerically (≈ 2.15627) into a proved bound for the
first time: 2.156224, from an explicit function evaluated in interval arithmetic. 70a replaces the
trivial bound 1 with 1.20192, using Szász's inequality to show that a cube and a ball are already in
their worst relative position in every dimension. 20a and 20c replace a classical example in the
table with a better one. All four arguments are in
[tao/PROOFS.md](tao/PROOFS.md).

### Packomania world records

[Packomania](http://www.packomania.com) lists the best known packings of N equal circles or spheres
in a container. Many of its large entries are beaten by a larger packing in the same table: removing
objects from the packing of M > N gives a packing of N with a bigger radius. We first produced a
packing of that kind for every such entry, then re-optimised each one with a compiled optimiser
(penalty method with L-BFGS, several deletion strategies, compression moves) running on 64-core
machines. The optimised packings beat the published records by a median of about 0.1%, around 2–3
times the gain from deletion alone.

| Category | Entries improved |
|---|---|
| Circles in a circle (`cci`) | see [packomania/SUMMARY.tsv](packomania/SUMMARY.tsv) |
| Circles in a square (`csq`) | |
| Spheres in a cube (`scu`) | |
| Spheres in a sphere (`ssp`) | |

Every packing is checked in exact rational arithmetic: every object inside the container, no two
overlapping (touching allowed), and a radius strictly larger than the published one. To check them
all:

```
python3 packomania/verify_all.py packomania/packings.zip
```

## How the search worked

The search system works in turns on each target. Each turn, a program is written for the target
and run for up to 40 minutes on a cloud machine; the next turn reads what it printed and what the
checker said, and improves on it. Every candidate is judged only by an independent checker, never by
the program that produced it. The loop is in [search/chain.py](search/chain.py), and the
packing optimiser it produced is [packomania/optimiser.py](packomania/optimiser.py).

## Repository layout

- `tao/`: the four certificates and their proofs.
- `packomania/`: the exact checker, the optimiser, the sweep that ran it over every entry, and the
  summary of all improved packings.
- `search/`: the turn-by-turn search loop and the sandbox that runs candidate programs without network
  access.
- `submissions/`: the prepared edits for Tao's repository (`tao-pr/`) and the Packomania submission
  (`packomania/`).

## Requirements

Python 3.10+ and `mpmath` (`pip install mpmath`). The packing checker needs only the standard library.
