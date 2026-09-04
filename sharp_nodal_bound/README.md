# Exact symbolic checks for the sharp nodal bound

This directory regenerates the exact computer-algebra identities used in the
proof of the aligned charge-two nodal bound

\[
\Xi_{\rm J}\ge 4.
\]

The manuscript prints the proof and its algebraic ingredients.  These scripts
are an independent machine-readable reproduction of the exact expansions and
eliminations; they are not a numerical replacement for the proof.

## Contents

- `verify_sharp_nodal_bound.py`
  - reconstructs the physical nodal polynomial `F_nod` from the JNR
    coefficient identity;
  - verifies the main `Delta_4` elimination identity exactly;
  - regenerates `Res_varpi(E4, Delta_{4,1})` and its quartic factor;
  - verifies the positive rational transform of that quartic;
  - verifies the initial derivative of the `E4=0` branch;
  - regenerates the endpoint resultant `Res_delta(P,Q)`;
  - converts the endpoint quartic to the degree-four Bernstein basis and
    verifies that all five coefficients are strictly negative;
  - checks the exact boundary test point and monotonicity derivatives.
- `verify_isosceles_sturm.py`
  - verifies the displayed Sturm chain for the equal-weight isosceles
    reconnection quintic up to positive rescalings;
  - regenerates the endpoint sign sequences, variation counts and rational
    root brackets.
  - This Sturm computation is auxiliary geometry; it is not used in the proof
    of `Xi_J >= 4` itself.
- `outputs/` contains frozen JSON outputs from the pinned environment.
- `verify_all.py` reruns both computations and compares them with those frozen
  outputs.

## Environment

Reference environment:

- Python 3.13.5
- SymPy 1.14.0

No floating-point arithmetic is used in the algebraic checks.

Create the environment with

```bash
conda env create -f environment.yml
conda activate jnr-sharp-nodal-bound
```

or install the single Python dependency with

```bash
python -m pip install -r requirements.txt
```

## Replay

From this directory:

```bash
python verify_all.py
```

Expected final line:

```text
ALL SYMBOLIC CHECKS PASS
```

For verbose exact output, run the two component scripts directly.

## Manuscript mapping

The checks are keyed to equation labels rather than printed equation numbers,
so later typesetting changes do not alter the mathematical correspondence.
The central script reproduces the objects labelled in the manuscript as
`eq:Fnod-Delta4-E4-identity`, `eq:node-E4-Delta41-resultant`,
`eq:G4-positive-transform`, `eq:res-P-Q-short`, and `eq:Q4-Bernstein`.
The auxiliary Sturm script reproduces `eq:isosceles-wall-quintic`,
`eq:isosceles-sturm-chain`, `eq:isosceles-sturm-signs`, and
`eq:isosceles-W5-root-brackets`.
