"""Two elementary improvements to lower bounds in Tao's table, checked in exact arithmetic.

20a (thin shell, variance form): C >= Var(|X|^2)/n for any isotropic log-concave X. The table's best
lower bound is 2 (the Gaussian). Take n = 1 and X = E - 1 with E ~ Exp(1): its density e^{-(x+1)} on
x > -1 is log-concave, E X = 0 and E X^2 = 1, so X is isotropic. Var(X^2) = E X^4 - 1, and the fourth
central moment of Exp(1) is 9, so C >= 8. (Products of n copies give the same ratio in every dimension.)

20c (KLS): C >= psi_mu >= min(mu(A), 1 - mu(A)) / mu^+(A) for any isotropic log-concave mu and set A.
The table's best lower bound is sqrt(pi/2) ~ 1.2533 (the Gaussian). Take n = 1, mu uniform on
[-sqrt 3, sqrt 3] (mean 0, variance 1, log-concave) and A = (-inf, 0]: mu(A) = 1/2 and mu^+(A) is the
density 1/(2 sqrt 3), so C >= sqrt 3 ~ 1.7321.
"""
from fractions import Fraction
from math import comb, factorial, sqrt, pi


def exp_central_moment(k):
    """E (E - 1)^k for E ~ Exp(1), exactly: sum_j C(k, j) (-1)^(k-j) E[E^j], with E[E^j] = j!."""
    return sum(comb(k, j) * (-1) ** (k - j) * factorial(j) for j in range(k + 1))


def thin_shell():
    m1, m2, m4 = exp_central_moment(1), exp_central_moment(2), exp_central_moment(4)
    assert (m1, m2) == (0, 1), "X = E - 1 must be isotropic"
    return Fraction(m4 - m2 ** 2)          # Var(X^2) / n with n = 1


def kls():
    # Uniform on [-a, a] with a = sqrt 3: variance a^2/3 = 1. Cut at 0.
    a2 = Fraction(3)
    assert a2 / 3 == 1
    mass = Fraction(1, 2)
    # boundary measure = density = 1/(2a); ratio = mass / density = a = sqrt 3
    return sqrt(3), "min(mu(A),1-mu(A)) / mu+(A) = (1/2) / (1/(2 sqrt 3)) = sqrt 3"


if __name__ == "__main__":
    v = thin_shell()
    print(f"20a: Var(|X|^2)/n = {v} for X = Exp(1) - 1, versus the table's 2  ->  improvement: {v > 2}")
    k, why = kls()
    print(f"20c: psi >= {k:.6f} ({why}), versus the table's sqrt(pi/2) = {sqrt(pi / 2):.6f}  ->  improvement: {k > sqrt(pi / 2)}")
