# C89 ≤ 0.0301388617424 (Ramsey multiplicity constant of K4, upper bound)

Step graphon with 3840 blocks of equal weight; red probability between blocks i, j is R_ij / 2³⁰
(`c89_stepgraphon_D30.json.gz`, schema of J. SK's certificate). Its monochromatic-K4 density is exactly

    31383471736748006978678023758588478740355521988220578341177988039 /
    1041295852679825698551191435837073366434387539971249757275291648000 = 0.03013886174230...

Previous best 0.0301388793133 (J. SK, a 120-step gradient refinement of Van Seters–Jothimani). Found by continuing
that refinement with L-BFGS on the same 15,739 symmetry classes of block pairs (`refine_lbfgs.py`).

**Exact evaluation** (`exact89.py`, written independently): the numerator Σ_i tr((diag(R_i)R)³) + the same for
B = 2³⁰ − R is computed over every block (no symmetry assumed) modulo 12 primes below 1.53·10⁶, where float64
products of residues are exact, combined by CRT, and confirmed by a 13th prime. Residues: `cand_part*.txt`
(about 6 minutes per prime on 64 cores). The same program reproduces J. SK's exact fraction
27220838459434583306277336063848530765605851167/903180180540490237232196590042258640992206848000 digit for
digit, and agrees with brute-force enumeration on small graphons.

    python3 exact89.py c89_stepgraphon_D30.json.gz part=0/4 64   # one quarter of the primes, 64 processes
    python3 exact89.py c89_stepgraphon_D30.json.gz residues=cand_part0.txt,cand_part1.txt,cand_part2.txt,cand_part3.txt

## Update: C89 ≤ 0.0301388584653 (second refinement round)

`c89_stepgraphon_D30_round2.json.gz`: the first-round graphon refined by a further L-BFGS round (same parametrisation),
rounded to multiples of 2⁻³⁰. Exact density, computed the same way (12 primes plus a check prime; residues in
`round2_part*.txt`):

    31383468324407329838135970970371081393990727093313271428891812741 /
    1041295852679825698551191435837073366434387539971249757275291648000 = 0.03013885846528...

    python3 exact89.py c89_stepgraphon_D30_round2.json.gz residues=round2_part0.txt,round2_part1.txt,round2_part2.txt,round2_part3.txt
