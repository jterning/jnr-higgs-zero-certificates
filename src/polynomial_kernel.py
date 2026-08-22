"""Symbolic definition of the denominator-free CKLL nodal kernel.

This module mirrors the polynomial definitions used by interval_certificate.py.
It is provided separately so readers can inspect or reuse the exact kernel
without reading the Bernstein subdivision implementation.
"""
from __future__ import annotations

import sympy as sp

p, a, b, y, u, v = sp.symbols("p a b y u v", real=True)
VARS = (p, a, b, y, u, v)

M = 1 - p**2
Sr = (1 + p) * a
Si = (1 - p) * b
qr = u**2 - v**2 - Sr*u + Si*v + p
qi = 2*u*v - Sr*v - Si*u
G = sp.expand(
    qr**2 + qi**2
    + y * (
        y
        + sp.Rational(1, 2) * ((2*u - Sr)**2 + (2*v - Si)**2)
        + 1 + p**2
        - sp.Rational(1, 2) * (Sr**2 + Si**2)
    )
)
D0 = 1 + y + u**2 + v**2 - 2*(a*u + b*v)
A0 = 1 - y - u**2 - v**2
L2 = sp.expand(A0**2 + 4*y)
H = sp.expand(A0*G - 2*M*y*D0)

Gy, Gu, Gv = sp.diff(G, y), sp.diff(G, u), sp.diff(G, v)
Hy, Hu, Hv = sp.diff(H, y), sp.diff(H, u), sp.diff(H, v)
Cu = sp.expand(Gy*Hu - Gu*Hy)
Cv = sp.expand(Gy*Hv - Gv*Hy)

Ds = sp.expand(M * (1 - a**2 - b**2))
Qs = sp.expand(3 + p**2 - (1+p)**2*a**2 - (1-p)**2*b**2)
Ktilde = sp.expand(M*Ds*H - Qs*G**2)


def cap_polynomial(delta: sp.Rational) -> sp.Expr:
    """Return Cap = Qs - X_delta Ds, where Xi_J <= X_delta means Cap <= 0."""
    delta = sp.Rational(delta)
    X = (8 + 3*delta) / (2 + 3*delta)
    return sp.expand(Qs - X*Ds)


def tilt_polynomial(h0: sp.Rational) -> sp.Expr:
    """Return Nh; at a node, h >= h0 implies Nh >= 0."""
    h0 = sp.Rational(h0)
    return sp.expand(H*Qs - h0**2 * M*Ds*L2)
