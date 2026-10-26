# C68 (Korenblum's constant) >= 223/500 = 0.446

Unrefereed, computer-assisted. The analytic part is below; the checkers in this
directory verify only the quantitative hypotheses for the exact rational certificate.

The argument keeps Wikström's moment-duality criterion [Wik2026, Lemma 3.1] and
Pinho's coefficient step [Pin2026]. It changes the geometric input in two places:

1. **Flux (period) lemma.** Single-valuedness of `ω` on the annulus forces an extra
   term in the Schwarz-type bound, which roughly doubles `1-K`. Pinho's covering-strip
   bound uses only two lifts of the zero. The sharpened bound matches the explicit
   degree-2 extremal to about 1e-5 relative (see `sanity_green.py`, `sanity_enclosing.py`, `sanity_outer.py`; exploratory, not part of the proof).
2. **Half-diameter enclosing lemma.** For the pairwise profile coming from (1), the
   image of each circle lies in a hyperbolic disk of radius `D/2`. Hyperbolic Jung would
   give `D/2 + log(2/sqrt3)`. Jung's theorem is not used.

## Notation

- `0<c<1`, `w=-log c`, `A=pi^2/w`, `q=e^{-A}`, `Psi=(1+q)^2/(8q)`.
- `Ann = {c<|z|<1}`. `G(.,b)` is its Green function, with `G ~ -log|z-b|` and `G=0` on the boundary.
- `S={-w<Re ζ<0}` is the covering strip, with `z=e^ζ`. The map `E(ζ)=exp(iπ(ζ+w)/w)` sends `S` biholomorphically to the upper half plane.
- For `ζ=x+iy`, put `φ=(x+w)/w ∈ (0,1)`. Then `E(ζ)=e^{-πy/w}e^{iπφ}`, and `φ(r)=1+log(r)/w` for `|z|=r`.
- For `t=r^2`, put `σ(t)=sin^2(π log t/(2w)) = sin^2(π φ(r))`.
- `δ(a,b)=|a-b|/|1-\bar b a|` is the pseudohyperbolic distance, and `d=2 artanh δ` is the hyperbolic distance.

## Lemma 1 (Green function lower bounds)

Let `G_S(ζ,β) = -log|(E(ζ)-E(β))/(E(ζ)-\overline{E(β)})|`.

(a) For any two lifts `β0, β1` of `b`, `G(e^ζ,b) >= G_S(ζ,β0)+G_S(ζ,β1)`.

*Proof.* `G̃(ζ)=G(e^ζ,b)` is positive and superharmonic on `S`. The difference
`G̃-G_S(.,β0)-G_S(.,β1)` is superharmonic: it is harmonic at `β0, β1` and has `+∞`
singularities at the other lifts. It is bounded below, and it has `liminf >= 0` at every
boundary point of `S` in the Riemann sphere, because `G_S → 0` at `∂S` and at `±i∞`. The
minimum principle finishes the proof.

(b) `-log x >= (1-x^2)/2` and `1-|(Z-B)/(Z-\bar B)|^2 = 4 Im Z Im B/|Z-\bar B|^2`.

Since `|Z-\bar B| <= |Z|+|B|`,

    G_S(ζ,β) >= (1/2) sin(πφ_ζ) sin(πφ_β) sech^2(π(y_ζ-y_β)/(2w)).

Take the two lifts adjacent to `ζ` in the `y` direction. Their half-separations are `u` and `A-u` with `u ∈ [0,A]`. If `sech^2(A/2) <= 1/3` (checked), then

    sech^2 u + sech^2(A-u) >= 2 sech^2(A/2).

For `u` or `A-u` below `u0=artanh(1/sqrt3)`, one term alone is at least `2/3`. Otherwise both
arguments lie where `sech^2` is convex, and Jensen's inequality applies.

Hence, **for all `a,b ∈ Ann`**:

    G(a,b) >= sin(πφ_a) sin(πφ_b) sech^2(A/2).                         (1)

(c) Same-radius case (Pinho's two-zero bound). For `|a|=|b|=r` and angular separation
`θ ∈ [0,2π]`, set `v=Aθ/π` and `p(v)^2 = sinh^2(v/2)/(sinh^2(v/2)+σ)`. Then (a) gives

    G(a,b) >= -log p(v) - log p(2A-v).

## Lemma 2 (flux / period lemma)

Let `ω: Ann → D` be holomorphic, `ω ≢ 0`, with `ω(z1)=0` and `|z1|=|z2|=r` (zeros counted with multiplicity, extra copies of `z1` included among the `b`). Then

    -log|ω(z2)| >= G(z2,z1) + σ sech^2(A/2),       σ = sin^2(πφ(r)).

*Proof.* We work on `A_ε={c e^ε<|z|<e^{-ε}}`, choosing `ε` so that no zero lies on `∂A_ε`,
and let `ε→0` at the end. All constants below are continuous in `ε`, and `G_ε ↑ G`.

Put `u=-log|ω|`, a positive superharmonic function that is continuous up to the boundary.
On `A_ε` it splits as

    u = Σ_b G_ε(.,b) + h,

where the sum runs over the finitely many zeros (with multiplicity) and `h` is harmonic and
positive, continuous up to the boundary.

For a circle `|z|=ρ` carrying no zero, the flux of a function `f` is the derivative of its
circle mean with respect to `x = log ρ`:

    Fl(f) = d/dx (1/2π)∮ f.

The flux identities used are these.

- `Fl(log|ω|)` is the winding number of `ω` on that circle (Cauchy–Riemann), so `Fl(u) ∈ Z`.
- The circle mean of `G(.,b)` is piecewise linear in `x`, vanishes at both ends, and has a kink of `-1` at `log|b|`. So `Fl(G(.,b))` is `1-φ_b` or `-φ_b`, and in either case `|sin(π Fl)| = sin(πφ_b)`.
- The circle mean of `h` is linear in `x`, with slope `F_h = (m_o-m_i)/w`, where `m_o, m_i` are its means on the two boundary circles.

Take `ρ` slightly larger than `r`. The flux of `G(.,z1)` is `-φ(r)`, so
`Σ_{b≠z1} Fl(G_b) + F_h ≡ φ(r) (mod 1)`.

Next we bound each part at `z2`, writing `κ = sin(πφ(r)) sech^2(A/2)`.

- By (1), `G(z2,b) >= κ · sin(πφ_b) = κ · |sin(π Fl(G_b))|`.
- By the Poisson representation in strip coordinates, `h(z2) = (1/2π)∫_{∂} h ∂_n G(z2,.) dy`.
- Letting `b → ∂` in (1) gives `∂_n G(z2,.) >= κ π/w`.
- Hence `h(z2) >= (κ π/w)(m_o+m_i) >= κ π |F_h| >= κ |sin(π F_h)|`.

Finally, `|sin(π x)|` is subadditive and 1-periodic, so

    V(z2) := u(z2)-G(z2,z1) >= κ (Σ|sin π Fl(G_b)| + |sin π F_h|) >= κ |sin(π φ(r))| = σ sech^2(A/2).

**Corollary 3.** Take any holomorphic `ω: Ann→D`, apply a disk automorphism so that `ω(z1)=0`, and use Lemma 1(c).
For `z1, z2` on `|z|=r` at angular separation `θ`:

    δ(ω(z1),ω(z2)) <= p(v) p(2A-v) e^{-E},     E = σ sech^2(A/2),  v = Aθ/π.

Write `u=v/2`. Using `1/p^2 = 1+σ/sinh^2` and `e^{2E} >= 1+2E` gives:

    cosh d(ω(z1),ω(z2)) - 1 <= (2/σ) Φ(u),
    Φ(u) = 1/(csch^2 u + csch^2(A-u) + 2 sech^2(A/2)).                 (2)

## Lemma 4 (simplex bound)

If `8 cosh^4(A/6) <= Psi` (checked), then:

- `Φ(u1)+Φ(u2)+Φ(u3) <= Psi` whenever `u_i>=0` and `Σu_i = A`;
- `Φ(u) <= Psi/2` for every `u`.

*Proof.* Use `csch^2 a >= 4e^{-2a}` and `2 sech^2(A/2) = 8q/(1+q)^2`. Writing `u = A/2+x`, these give

    Φ(A/2+x) <= (Psi/2) sech^2 x.

Let `u3` be the smallest of the three, so `u3 <= A/3`. Then `x1+x2 = -u3`, so `max|x_i| >= u3/2`, and also `Φ(u3) <= sinh^2 u3`. Therefore

    ΣΦ <= Psi - (Psi/2) tanh^2(u3/2) + sinh^2 u3.

The bound `ΣΦ <= Psi` holds because `sinh^2 u3 <= (Psi/2) tanh^2(u3/2)` is equivalent to `8 cosh^4(u3/2) <= Psi`, and `u3 <= A/3`.

## Lemma 5 (half-diameter enclosing disk)

The compact set `Γ = ω(|z|=r)` lies in a closed hyperbolic disk of radius `ρ` with

    sinh^2 ρ <= Psi/(2σ).

*Proof.* Take the minimal enclosing disk, with centre `o`. The directions from `o` to the
contact points do not lie in any open half-plane, since otherwise moving `o` would shrink the
disk. So one of two cases holds.

- **Two antipodal contacts `p, q`.** Then `d(p,q)=2ρ`, so `2 sinh^2 ρ = cosh 2ρ - 1 <= (2/σ)Φ <= Psi/σ`.
- **Three contacts with consecutive central angles `γ_i < π`, `Σγ_i = 2π`.** The cosine rule gives `cosh d_ij - 1 = sinh^2ρ (1-cos γ_ij)`. Since `Σ cos γ_i <= -1` on that set, `Σ(1-cos γ) >= 4`. The three contact points are images of three points on the circle, whose arcs give `u`'s summing to `A`. So `4 sinh^2 ρ <= (2/σ) ΣΦ <= 2 Psi/σ` by (2) and Lemma 4.

## Lemma 6 (coefficient estimate, Pinho's step)

Suppose `f=Σa_k z^k, g=Σb_k z^k ∈ A^2(D)` with `|f| <= |g|` on `Ann` (if `g≡0` or `f/g` is a unimodular constant the conclusion is trivial; otherwise `ω=f/g` maps `Ann` into `D`), and write `Y(t)=Σ_{k>=1, |a_k|<|b_k|}(|b_k|^2-|a_k|^2) t^k`. Then

    Y(r^2) <= G(r^2)(M_g(r)^2 - M_f(r)^2),
    G(t) = (1+sqrt(1+Psi/(2σ(t))))/2.

*Proof.* Let `α` be the centre of the disk from Lemma 5, with pseudohyperbolic radius `tanh(ρ/2)`.
Then `1/(1-tanh^2(ρ/2)) = (1+cosh ρ)/2 <= G`. Set `F=f-αg` and `H=g-\bar α f`. These satisfy

    |H_k|^2-|F_k|^2 = (1-|α|^2)(|b_k|^2-|a_k|^2),

and the rest of the argument is Pinho's [Pin2026, Lemma 2] verbatim.

## Lemma 7 (Wikström's moment criterion, unchanged)

Let `β = 1-1/G` and let `μ` be a probability measure on `[c^2,1]` with no atom at 1. Suppose

- `∫t^k dμ >= 1/(k+1)` for `k>=0`;
- `∫β t^k dμ <= 1/(k+1)` for `k>=1`.

Then `c` is admissible.

## Certificate and what the checkers verify

The measure is `μ = Σ_j w_j δ_{t_j} + 1_{[R,1]}dt`, with exact rationals in
`certificate-c0446.json`. Its data:

- `c=223/500` and `R=99/100`;
- 8 atoms, with `t_1=c^2` and `t_8=R`;
- `Σw_j=R` exactly;
- a tail partition `R=τ_0<…<τ_N=1`, where `τ_i = 1-(1-R)(4/5)^i` for `i<N` and `τ_N=1`.

**Tail lower bound.** On `[R,1]`, `g=1/G` is decreasing, because `σ` decreases in `t` on `[c,1]`. On a piece `[a,b]`, using `-log t >= δ(1+δ_lo/2)` and `sin x >= x(1-U^2/6)`, one gets

    g(t) >= l_i (1-t),
    l_i = 2√2 b_i/(√2 b_i δ_hi + sqrt(2 b_i^2 δ_hi^2+Psi)),
    b_i = (π/2w)(1+δ_lo/2)(1-U_i^2/6),    U_i = -π log a/(2w) < √6.

The checkers verify, with outward-rounded interval arithmetic:

- the structural conditions (exact rationals);
- `sech^2(A/2) < 1/3`;
- `8cosh^4(A/6) < Psi`;
- `A_k > 0` for `1<=k<=5806`. For `k>=5807`, `w_8(k+1) >= R` makes the atom at `R` sufficient;
- `B_k > 0` for `1<=k<=1349`. For `k>=1350`, one tail piece `[τ_i,τ_{i+1}]` gives `l_i(1-τ_{i+1})(τ_{i+1}-τ_i)τ_i^k > Σ w_j t_j^k`. This propagates to all larger `k` because `t_j <= R <= τ_i`.

The minimum relative margins are `6.990e-6` for `A_k` and `6.982e-6` for `B_k`.

`certificate-c044605.json` proves `c=8921/20000=0.44605` the same way, with margins
`8.47e-7` and `8.16e-7`. Exploratory floating-point linear programs place the limit of this
whole framework at about `0.44607`.
