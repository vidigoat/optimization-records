# C41b ≤ 2√2 − 1 ≈ 1.8284271 (ambidextrous moving sofa, upper bound)

Previous upper bound: 2.2195, inherited from the ordinary sofa. Best lower bound: Romik's 1.64495521.

**Setting.** We use the standard formulation (Gerver; Kallus–Romik 2018). An ambidextrous sofa S, in its own frame,
lies in the strip H = ℝ × [0,1]. Both S and its mirror image are moving sofa shapes starting from H, with monotone
rotation paths whose total angles are β₁, β₂ ∈ (0, π/2].

**Proof.**

*Small turns.* If βᵢ < π/4, then S lies in H and in a unit strip at angle 90° − βᵢ to H. Their intersection has
area 1/cos βᵢ < √2.

*Both turns at least π/4.* Look at the two moments when the rotation equals π/4. At those moments, S lies in H, in
the left corridor turned by 45°, and in the right corridor turned by 45°. Switch to the coordinates s = x − y,
t = x + y. Areas double. Put r = √2. Then:

- H becomes the band D = {0 ≤ t − s ≤ 2}.
- The two corridors become axis-aligned L-shapes:
  - L₁ = {s ≥ a, t ≤ b} minus {s > a+r, t < b−r};
  - L₂ = {s ≤ d, t ≥ c} minus {s < d−r, t > c+r}.

We show that F = D ∩ L₁ ∩ L₂ has area at most 4√2 − 2. Write δ = d − a and η = b − c. Every point of F is of one
of four types:

- P: s ≤ a+r and t ≤ c+r. Then P lies in an r-square.
- Q: t ≥ b−r and s ≥ d−r. Then Q lies in an r-square.
- U: s ∈ [d−r, a+r] and t ∈ (c+r, b−r). This needs δ ≤ 2r and η > 2r.
- V: s ∈ (a+r, d−r) and t ∈ [b−r, c+r]. This needs δ > 2r and η ≤ 2r.

*Two lemmas.*

- (A) An r-square meets D in area at most 2√2 − 1.
- (B) A rectangle of width w ≤ r in s and height r in t meets D in area at most
  - g(w) = wr when w ≤ 2 − r, and
  - g(w) = wr − (w + r − 2)²/4 otherwise.

  In both ranges, g(w) ≤ √2 − 1 + w, with equality at w = r. Lemma A is the case w = r.

To prove (B): the band cuts off two corners, and the cut area is a convex function of where the cut falls. So the
total cut is smallest when both corners lose the same amount. The inequality g(w) ≤ √2 − 1 + w then follows because
h(w) = √2 − 1 + w − g(w) is convex on [2 − r, r], with h(r) = h′(r) = 0.

*Cases.*

1. δ > 2r and η > 2r: U and V are empty, so F lies in P ∪ Q, and area(F) ≤ 2(2√2 − 1) by (A).
2. δ ≤ 2r and η ≤ 2r: again U and V are empty, so area(F) ≤ 2(2√2 − 1).
3. δ ≤ 2r and η > 2r: V is empty.
   - If δ ≤ r, then F lies in a column of width δ, so its area is at most 2δ ≤ 2√2.
   - Otherwise put w = δ − r. Then F lies in the union of a column of width r − w and two w × r rectangles. Its
     area is at most 2(r − w) + 2g(w) ≤ 2r + 2√2 − 2 = 4√2 − 2.
4. δ > 2r and η ≤ 2r: the map (s,t) ↦ (−t,−s) fixes D, keeps the shapes of L₁ and L₂, and swaps δ and η. So this
   is case 3.

Dividing by 2 gives area(S) ≤ 2√2 − 1. ∎

**Remarks.**
- This is the ambidextrous analogue of Hammersley's 2√2. The bound is attained by the 45° relaxation: the extremal
  region is two hexagons of area √2 − 1/2 each.
- For motions that rotate the wrong way, a numerical search gives at most √2 in the three other angle cases.
- Using more angles, as Kallus and Romik do, should lower the bound further.
- `check_sofa.py` cross-checks the (s,t) reduction against direct polygon geometry; they agree to 1.4·10⁻¹⁵. It also
  searches the parameters at random, and the largest area it finds is 1.82815. Its output is in `check_output.txt`.

Found and written up with Claude (Anthropic) agents under V. Patankar's direction. The proof has not yet been checked
by a human expert.
