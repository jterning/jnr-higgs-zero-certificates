#!/usr/bin/env python3
"""Exact check of the displayed equal-weight isosceles Sturm chain.

This is auxiliary to the sharp Xi_J >= 4 proof: the sharp proof itself uses
resultants and a Bernstein sign certificate, not this Sturm chain.  It is
included because the manuscript uses the Sturm computation to locate the two
isosceles reconnection points.
"""
from __future__ import annotations

import json
import sympy as sp


def sign(v):
    v = sp.sign(sp.simplify(v))
    if v == 1:
        return "+"
    if v == -1:
        return "-"
    return "0"


def variations(signs):
    nz = [s for s in signs if s != "0"]
    return sum(a != b for a, b in zip(nz, nz[1:]))


def compute():
    c = sp.symbols("c", real=True)
    W5 = 400*c**5 + 2544*c**4 + 5751*c**3 + 4658*c**2 - 351*c - 1338
    T = [
        W5,
        sp.diff(W5, c),
        60078*c**3 + 332159*c**2 + 552248*c + 260147,
        -(5403616793*c**2 + 12646417916*c + 5004643589),
        -(9647529398911*c + 16526296022516),
        -sp.Integer(1),
    ]

    positive_scales = []
    for i in range(len(T)-2):
        neg_rem = -sp.rem(T[i], T[i+1], domain=sp.QQ)
        ratio = sp.factor(neg_rem/T[i+2])
        assert not ratio.has(c)
        assert ratio > 0
        assert sp.expand(neg_rem - ratio*T[i+2]) == 0
        positive_scales.append(ratio)

    signs_minus = [sign(poly.subs(c, -1)) for poly in T]
    signs_plus = [sign(poly.subs(c, 1)) for poly in T]
    assert signs_minus == ["+", "-", "-", "+", "-", "-"]
    assert signs_plus == ["+", "+", "+", "-", "-", "-"]
    vminus = variations(signs_minus)
    vplus = variations(signs_plus)
    assert (vminus, vplus) == (3, 1)
    assert vminus - vplus == 2

    brackets = [
        (sp.Rational(-4551,5000), "+"),
        (sp.Rational(-9101,10000), "-"),
        (sp.Rational(2199,5000), "-"),
        (sp.Rational(4399,10000), "+"),
    ]
    bracket_signs = [[str(q), sign(W5.subs(c,q))] for q,_ in brackets]
    assert [s for _,s in bracket_signs] == [s for _,s in brackets]

    return {
        "status": "PASS",
        "sympy_version": sp.__version__,
        "sturm_positive_rescaling_factors": [str(v) for v in positive_scales],
        "endpoint_signs": {"c=-1": signs_minus, "c=+1": signs_plus},
        "variation_counts": {"c=-1": vminus, "c=+1": vplus, "roots_in_open_interval": vminus-vplus},
        "rational_bracket_signs": bracket_signs,
    }


def main():
    print(json.dumps(compute(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
