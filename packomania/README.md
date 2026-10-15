# Packomania packings

[Packomania](http://www.packomania.com), maintained by Eckard Specht, lists the best known packings of
N equal circles or spheres in a container, for N up to 10,000.

**The observation.** In four of its tables, more than 2,000 entries have a smaller radius than a
packing with *more* objects in the same table. Deleting objects from that larger packing already
gives a better packing for the smaller N. We started from those packings and re-optimised each one.

**The result.** 2,072 packings with a larger radius than the published record:

| Category | Entries improved | Median gain | Largest gain |
|---|---|---|---|
| Circles in a circle (`cci`) | 193 | 0.061% | 0.25% |
| Circles in a square (`csq`) | 1813 | 1.703% | 6.71% |
| Spheres in a cube (`scu`) | 49 | 0.020% | 1.10% |
| Spheres in a sphere (`ssp`) | 17 | 2.482% | 3.31% |
| **Total** | **2072** | | |

Every entry, with the published and the new radius, is in [SUMMARY.tsv](SUMMARY.tsv). The packings
themselves (about 54 MB, in Packomania's `.pck` submission format) are attached to the
[v1.0 release](https://github.com/vidigoat/optimization-records/releases/tag/v1.0).

**Checking.** [`verify.py`](verify.py) checks a packing in exact rational arithmetic: every object
inside the container, no two overlapping (touching allowed). [`verify_all.py`](verify_all.py) runs it
on every packing in the zip and confirms each radius beats the published one:

    python3 packomania/verify_all.py packings.zip

**Code.** [`optimiser.py`](optimiser.py) is the optimiser (a compiled penalty method with L-BFGS,
several deletion strategies and compression moves, run with up to 60 processes). `tools/` holds the
scripts that found the dominated entries (`harvest.py`), ran the optimiser over all of them
(`sweep.py`) and built the submission (`package.py`); they are reference copies from the full
project and expect its layout.
