# C1a ≤ 1.502743356212 (Sidon set autocorrelation constant, upper bound)

Non-negative step function with n = 131072 equal steps on [−1/4, 1/4] (`c1a_n131072.json`, heights as decimal
strings read as exact rationals). Score 2n · max_k (a∗a)_k / (Σ a)², exactly:

    88396668012481333418540698869863464003678707923912272976102136217600000 /
    58823529411764719617703895028066775569135640608648734518121194077814337 = 1.5027433562121823...

Previous best 1.5027434868 (Atris, 65,536 steps). Found by splitting each step of that construction in two and
re-optimising with L-BFGS on a soft maximum (`refine_ab.py`). `exact_ab.py 1a c1a_n131072.json` computes the
full autoconvolution exactly with one big-integer product (standard library); it reproduces Atris's value
1.502743486709583... exactly on their construction.
