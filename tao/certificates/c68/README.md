# C68 ≥ 223/500 = 0.446 (Korenblum's constant, lower bound)

Unrefereed, computer-assisted. `PROOF.md` gives the analytic argument; the checkers verify the quantitative hypotheses
for the exact-rational certificate `certificate-c0446.json`.

The argument keeps Wikström's moment-duality criterion and Pinho's coefficient step, and sharpens the geometric input:

1. **Flux (period) lemma** (`PROOF.md`, Lemma 2). Single-valuedness of ω on the annulus forces the fluxes of the Green
   potentials of its zeros, plus that of the harmonic part, to sum to an integer; with a lower bound on the annulus
   Green function this gives `−log|ω(z₂)| ≥ G(z₂,z₁) + σ·sech²(A/2)` for zeros on the same circle.
2. **Half-diameter enclosing lemma** (Lemma 5): the image of each circle lies in a hyperbolic disk of radius D/2
   (hyperbolic Jung would give D/2 + log(2/√3)).

Together: `G(t) = (1 + sqrt(1 + Ψ/(2σ(t))))/2`, used in Pinho's coefficient step and Wikström's criterion. The
certificate has c = 223/500, R = 99/100, 8 atoms and Lebesgue measure on [R, 1]; minimum relative margins 6.99·10⁻⁶.
`certificate-c044605.json` proves 0.44605 the same way with smaller margins; the framework's limit is about 0.44607.

**Checkers** (independent of each other and of the optimiser): `check_c68.py` (Arb ball arithmetic, 320 bits) and
`check2_c68.py` (mpmath interval arithmetic). Both pass the new certificate, and both replay Pinho's 0.432 certificate
(`pinho-certificate-jung-0432.json`) with his margins. `mutation_tests.py` checks that both reject seven corrupted
certificates. `exploratory/` holds the search and numerical sanity scripts (not part of the proof).

    python3 check_c68.py certificate-c0446.json
    python3 check2_c68.py certificate-c0446.json
