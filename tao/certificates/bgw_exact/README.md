# C16 >= 2.156225504429 (Brezis–Gallouët–Wainger constant, Tao's table entry 16a)

Same function as `../bgw.py`: on the torus [0, 2π]², u(x) = Σ c_k e^{ik·x} over 0 < |k|² ≤ K², with
c_k = 1/(n(n+μ)), n = |k|², μ = 41/5, and

    L ≥ S0²/(π S1) − ln δ − ln(1 + ln δ),   δ = S2/S1 ≥ 1,

with S0 = Σ c_k, S1 = Σ n c_k², S2 = Σ n² c_k². Only the cutoff changes: K = 40000 instead of 4000.

For a fixed set of frequencies this shape is optimal: setting the gradient of the bound to zero gives
c_n ∝ 1/(n(n+μ)) for some μ, so the only remaining freedom is the cutoff and μ (μ near 8.19–8.20 is best;
41/5 is within 2·10⁻¹¹ of the optimum). The gain from larger K decays like 12/K²; the limit of this
family is about 2.1562255120.

**Arithmetic.** At 5·10⁹ lattice points, floating-point interval sums are too loose, so `exact.c` floors
every term exactly with 128-bit integer division at scale 2^100 (terms are q/(n(qn+p)), q²/(n(qn+p)²) and
q²/(qn+p)² with μ = p/q), so each sum is known to within (number of terms)·2^−100. `combine.py` turns the
integer sums into exact rationals and evaluates the bound in interval arithmetic.

**Checks.** At K = 800 and K = 4000 this code reproduces the earlier interval-arithmetic certificate
(2.156206556532 and 2.156224753794, with 50,265,328 lattice points at K = 4000, the count found in the
independent replay on PR #226).

    ./run.sh          # K = 40000, about 10 minutes on 14 cores; prints C16 >= 2.15622550442903481
    ./run.sh 4000     # one second
    python3 combine.py sums_K40000.txt   # re-evaluate from the stored integer sums

Together with the upper bound C16 ≤ 2.1562256466* ([maximusrufus/c16-bgw-certificates]) this places C16 in an
interval of width 1.4·10⁻⁷.
