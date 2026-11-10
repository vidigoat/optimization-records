"""Piecewise-linear densities for the Schur-Siegel-Smyth trace problem (C86), floating-point search side.

A density is p = sum_k h_k phi_k, phi_k the hat with knots (e_{k-1}, e_k, e_{k+1}) and height 1 at e_k.
With H(w) = w^2 (log w / 2 - 3/4) (so H'' = log w), integrating by parts twice gives
    U_p(z) = int log|z - y| p(y) dy = Re sum_k s_k H(z - e_k),
where s_k are the slope jumps of p; this holds for complex z (principal log, z off the support's real line
or on it). mass(phi_k) = (e_{k+1} - e_{k-1}) / 2, mean(phi_k) = (e_{k-1} + e_k + e_{k+1}) / 3.
"""
import numpy as np


def H(w):
    w = np.asarray(w, dtype=complex)
    out = np.zeros_like(w)
    nz = w != 0
    out[nz] = w[nz] ** 2 * (np.log(w[nz]) / 2 - 0.75)
    return out


def hat_jumps(e):
    """Slope jumps of each interior hat: rows k = 1..n-2 of e, three knots each."""
    e = np.asarray(e, float)
    l, c, r = e[:-2], e[1:-1], e[2:]
    return np.stack([1 / (c - l), -1 / (c - l) - 1 / (r - c), 1 / (r - c)], axis=1), np.stack([l, c, r], axis=1)


def hat_potentials(e, z):
    """Matrix U[i, k] = U_{phi_k}(z_i) for interior hats k (knots e[k], e[k+1], e[k+2])."""
    s, kn = hat_jumps(e)
    z = np.asarray(z)
    out = np.zeros((len(z), len(s)))
    for j in range(3):
        out += (s[:, j][None, :] * H(z[:, None] - kn[:, j][None, :])).real
    return out


def hat_mass_mean(e):
    e = np.asarray(e, float)
    l, c, r = e[:-2], e[1:-1], e[2:]
    m = (r - l) / 2
    return m, m * (l + c + r) / 3


def poly_logabs(coeffs_low_first, x):
    """log|Q(x)| for real x (coefficients constant first)."""
    return np.log(np.abs(np.polynomial.polynomial.polyval(np.asarray(x, float), np.asarray(coeffs_low_first, float))))


def poly_integral_rows(coeffs_low_first, e):
    """Row I[k] = int log|Q| d phi_k  = mass_k log|lc| + sum_roots U_{phi_k}(alpha)."""
    c = np.asarray(coeffs_low_first, float)
    roots = np.polynomial.polynomial.polyroots(c)
    m, _ = hat_mass_mean(e)
    return m * np.log(abs(c[-1])) + hat_potentials(e, roots).sum(axis=0)
