#!/usr/bin/env python3
"""Exact symbolic checks for the sharp aligned nodal bound Xi_J >= 4.

All arithmetic is exact SymPy arithmetic.  The script reconstructs the nodal
polynomial from the physical JNR coefficient identity, verifies the elimination
identity, the two resultants/factorizations, and the Bernstein sign certificate
used in the sharp nodal-bound proof.
"""
from __future__ import annotations

import json
from pathlib import Path
import sympy as sp


def bernstein_coefficients(poly: sp.Expr, var: sp.Symbol, degree: int):
    coeffs = sp.symbols(f"b0:{degree+1}")
    basis = sum(
        coeffs[k] * sp.binomial(degree, k) * var**k * (1-var)**(degree-k)
        for k in range(degree + 1)
    )
    eq_poly = sp.Poly(sp.expand(basis - poly), var)
    equations = [c for c in eq_poly.all_coeffs()]
    sol = sp.solve(equations, coeffs, dict=True)
    if len(sol) != 1:
        raise RuntimeError("Bernstein conversion was not unique")
    return [sp.factor(sol[0][c]) for c in coeffs]


def compute():
    p, x, C = sp.symbols("varpi x C", real=True)

    # u_* = sqrt(x) exp(i theta), C = cos^2(theta), so cos(2 theta)=2C-1.
    cos2 = 2*C - 1
    absS2 = sp.expand(
        x * ((2-x)**2 + p**2 + 2*p*(2-x)*cos2)
    )
    re_pbarS2 = sp.expand(
        p*x * (((2-x)**2 + p**2)*cos2 + 2*p*(2-x))
    )

    # Physical positive-weight JNR level mu_JNR(S,P)=N/D with P=varpi real.
    numerator = sp.expand(
        1 - absS2 - 2*p**2 + p**4 - p**2*absS2 + 2*re_pbarS2
    )
    denominator = sp.expand(3 + p**2 - absS2)

    delta_star = sp.expand((p+x)**2 - 4*p*x*C)
    mu_node = (1-x) * delta_star
    raw_residual = sp.factor(numerator - mu_node*denominator)
    F_nod = sp.cancel(raw_residual/(1-x))
    assert sp.denom(F_nod) == 1

    Delta4 = sp.expand(1 - p**2 - 4*(1-x)*delta_star)
    b = sp.expand(2*(p**2*(5*x-6) - 4*x**3 + 14*x**2 - 21*x + 12))
    E4 = sp.expand(
        (10-9*x)*p**4
        + (8*x**3 - 28*x**2 + 50*x - 32)*p**2
        + 16*x**5 - 80*x**4 + 152*x**3 - 132*x**2 + 39*x + 6
    )

    elimination_residual = sp.expand(
        16*(1-x)**2*F_nod - ((x-2)*Delta4**2 + b*Delta4 - E4)
    )
    assert elimination_residual == 0

    # The positivity estimate b/2 >= 2(1-x)^2(3-2x) follows from this exact factorization.
    b_gap = sp.factor(b/2 - 2*(1-x)**2*(3-2*x))
    expected_b_gap = (1-p**2)*(6-5*x)
    assert sp.expand(b_gap - expected_b_gap) == 0

    Delta41 = sp.expand(Delta4.subs(C, 1))
    H4 = 225*x**4 - 2860*x**3 + 8254*x**2 - 9352*x + 3781
    resultant_E4_Delta41 = sp.factor(sp.resultant(E4, Delta41, p))
    expected_res1 = 256*x**2*(1-x)**8*H4
    assert sp.expand(resultant_E4_Delta41 - expected_res1) == 0

    r = sp.symbols("r", nonnegative=True)
    H4_transform = sp.expand(sp.cancel((1+r)**4 * H4.subs(x, r/(1+r))))
    expected_transform = 48*r**4 + 716*r**3 + 2884*r**2 + 5772*r + 3781
    assert sp.expand(H4_transform - expected_transform) == 0

    # Initial derivative of Delta_{4,1} along the E4=0 branch.
    dpdx = -sp.diff(E4, x)/sp.diff(E4, p)
    total_derivative = sp.simplify(sp.diff(Delta41, x) + sp.diff(Delta41, p)*dpdx)
    initial_derivative = sp.simplify(total_derivative.subs({x: 0, p: 1/sp.sqrt(5)}))
    expected_initial = 4*(-69 + 14*sp.sqrt(5))/35
    assert sp.simplify(initial_derivative - expected_initial) == 0

    # Endpoint sector: delta=1-x, y=(varpi-x)/delta.
    d, y = sp.symbols("delta y", real=True)
    Pcert = sp.expand(4*d**2*y**2 + d*(1-y)**2 - 2*(1-y))
    Qcert = sp.expand(d**2*y**4 + d*y**3*(4-y) - 4*y**3 + 5*y**2 - 1)
    Theta4 = 2*y**4 - 73*y**3 + 103*y**2 - 35*y - 15
    resultant_PQ = sp.factor(sp.resultant(Pcert, Qcert, d))
    expected_res2 = -y**4*(y-1)*(5*y**2-1)*Theta4
    assert sp.expand(resultant_PQ - expected_res2) == 0

    bern = bernstein_coefficients(Theta4, y, 4)
    expected_bern = [sp.Integer(-15), -sp.Rational(95,4), -sp.Rational(46,3), sp.Integer(-8), sp.Integer(-18)]
    assert bern == expected_bern
    assert all(c < 0 for c in bern)

    # Monotonicity polynomials used after the boundary resultant.
    dP_dd = sp.factor(sp.diff(Pcert, d))
    dQ_dd = sp.factor(sp.diff(Qcert, d))
    assert sp.expand(dP_dd - (8*d*y**2 + (1-y)**2)) == 0
    assert sp.expand(dQ_dd - (2*d*y**4 + y**3*(4-y))) == 0

    # Exact positive test point on the connected boundary branch.
    delta0 = (sp.sqrt(65)-1)/8
    assert sp.simplify(Pcert.subs({d: delta0, y: sp.Rational(1,2)})) == 0
    q_test = sp.simplify(Qcert.subs({d: delta0, y: sp.Rational(1,2)}))
    expected_q_test = (27*sp.sqrt(65)-123)/512
    assert sp.simplify(q_test - expected_q_test) == 0
    assert q_test > 0

    return {
        "status": "PASS",
        "sympy_version": sp.__version__,
        "checks": {
            "reconstructed_F_nod_from_physical_JNR_identity": True,
            "main_elimination_identity": True,
            "b_positive_gap_factorization": str(b_gap),
            "resultant_E4_Delta41": str(resultant_E4_Delta41),
            "H4_positive_transform": str(H4_transform),
            "E4_branch_initial_derivative": str(initial_derivative),
            "resultant_P_Q": str(resultant_PQ),
            "Theta4_Bernstein_coefficients": [str(v) for v in bern],
            "boundary_positive_test": str(q_test),
        },
    }


def main():
    result = compute()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
