# Optimization records

**New best-known results in three public record tables of mathematics, each with a script that checks it.**

*Vidit Patankar · October 2026*

| | Table | Result |
|---|---|---|
| 📐 | [Terence Tao's optimization constants](https://github.com/teorth/optimizationproblems) | **5 constants improved**, including a first rigorous bound for 16a and a new cube-and-ball proof for 70a ([pull request #226](https://github.com/teorth/optimizationproblems/pull/226)) |
| 🔵 | [Packomania](http://www.packomania.com) (best packings of circles and spheres) | **2,072 packings** beat the published records, by up to 6.7% |
| 📡 | [codetables.de](http://www.codetables.de) (best error-correcting codes) | **7 record codes**, beating entries unchanged since 2001–2008 |

Everything here can be re-checked on a laptop:

```sh
pip install mpmath
./check.sh          # Tao certificates and the four codes, a few minutes
./check.sh --full   # also 16a at full strength and all 2,072 packings
```

---

## 📐 Tao's table of optimization constants

| Entry | Constant | Listed bound | New bound |
|---|---|---|---|
| 16a | Brezis–Gallouët–Wainger remainder constant on the 2D torus | ≥ 1.82283 | **≥ 2.156224** |
| 70a | Reverse Brunn–Minkowski constant | ≥ 1 | **≥ 1.201927** |
| 23b | Decay exponent of the minimal condition number of ±1 matrices | ≤ 1 | **≤ 1/2** |
| 20a | Thin-shell constant (variance form) | ≥ 2 | **≥ 8** |
| 20c | KLS constant | ≥ 1.2533 | **≥ 1.7321** |

- **16a.** Bartuccelli, Deane and Zelik (2013) computed L ≈ 2.15627 numerically, but the best proved
  bound was 1.82283. An explicit trigonometric polynomial, with its sums and logarithms evaluated in
  interval arithmetic, proves L ≥ 2.156224, within 0.00005 of the numerical value.
- **70a.** The listed lower bound was the trivial 1. For a cube and a ball, Szász's inequality on the
  Gram matrix shows no linear map can make their Minkowski sum smaller. The resulting bound grows with
  the dimension and reaches 1.2019.
- **23b.** For odd n, every off-diagonal entry of AᵀA is odd. The Bhatia–Davis inequality then forces
  κ(n) ≥ 1 + 0.58/√n, so the decay exponent is at most 1/2, not 1.
- **20a, 20c.** The table listed the Gaussian. A centred exponential variable gives 8, and the uniform
  distribution on an interval gives √3.

Proofs: [`tao/PROOFS.md`](tao/PROOFS.md). Certificates: [`tao/certificates/`](tao/certificates).

## 🔵 Packomania

More than 2,000 entries in four of Packomania's tables have a smaller radius than a packing with
*more* objects in the same table, so removing objects already beats the record. Starting from those
packings, a compiled optimiser on 64-core machines grew every radius further.

| Category | Entries improved | Median gain | Largest gain |
|---|---|---|---|
| Circles in a circle | 193 | 0.061% | 0.25% |
| Circles in a square | 1,813 | 1.703% | 6.71% |
| Spheres in a cube | 49 | 0.020% | 1.10% |
| Spheres in a sphere | 17 | 2.482% | 3.31% |
| **Total** | **2,072** | | |

Every packing is checked in exact rational arithmetic. Details and the full list are in
[`packomania/`](packomania); the packings are in the
[v1.0 release](https://github.com/vidigoat/optimization-records/releases/tag/v1.0).

## 📡 codetables.de

Seven linear codes over GF(5), found by Construction X on nested constacyclic codes:

| Code | Listed | New | Upper bound | Entry last changed |
|---|---|---|---|---|
| [27, 10] | 12 | **13** | 14 | 2001 |
| [64, 10] | 38 | **39** | 44 | 2003 |
| [81, 10] | 50 | **51** | 57 | 2008 |
| [82, 10] | 51 | **52** | 58 | 2003 |
| [67, 12] | 37 | **38** | 45 | 2003 |
| [68, 12] | 38 | **39** | 45 | 2003 |
| [69, 12] | 38 | **39** | 46 | 2003 |

Each minimum distance is computed exactly by enumerating every codeword. Details:
[`codetables/`](codetables).

---

## How the results were found

Each target got its own *chain* of turns. Each turn, a search program was written, ran for up to
40 minutes on a 64-core cloud machine without network access, and was judged only by an exact checker
written separately. The next turn read its output and improved on it. Programs were written with
Claude Opus through Anthropic's API. See [`search/`](search).

**What did not work.** The same system also tried to beat the record for the Shannon capacity of the
7-cycle (367 words in C₇⁵; about 1,500 programs, all stopping at exactly 367), a Turing machine for
Goldbach's conjecture with fewer than 25 states, finite counterexamples to several Erdős problems,
and matrix multiplication schemes with fewer multiplications. None of these improved a record.

## Layout

```
tao/            proofs and certificates for the five constants
packomania/     exact packing checker, optimiser, summary of all 2,072 entries
codetables/     the four generator matrices and an exact minimum-distance checker
search/         the turn-by-turn search loop, sandbox and budgeted API calls
submissions/    the prepared pull request for Tao's table and emails for the other two tables
check.sh        re-checks everything
```

MIT licence. To cite, see [`CITATION.cff`](CITATION.cff).
