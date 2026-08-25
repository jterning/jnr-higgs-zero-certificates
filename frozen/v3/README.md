# Frozen Bernstein certificate snapshot

This archive freezes the theorem-level near-aligned Bernstein certificate used by the accompanying charge-two dyonic-instanton manuscript.

## One-command verification

From the archive root:

```bash
python scripts/verify_exact_v3.py
```

The verifier reconstructs the five denominator-free polynomials and their initial tensor-product Bernstein coefficients from symbolic formulas, subdivides the deliberately enlarged rational box, and compares all invariant output metadata with the two frozen JSON results. Timing fields are ignored in the comparison.

## Certified cases

- `delta=1/50`, `h0=24/25`: `CERTIFIED_EMPTY`, 37,583 visited boxes.
- `delta=1/20`, `h0=9/10`: `CERTIFIED_EMPTY`, 30,159 visited boxes.

## Arithmetic assurance

Initial Bernstein coefficients are exact `fractions.Fraction` values. Binary64 values are used only as exact dyadic centers. The initial conversion error and every propagated de Casteljau error radius are exact rationals. At each pruning decision the selected binary64 minimum or maximum is converted back with `Fraction.from_float` and compared against the exact rational radius. The child radius is

`E_child = E_parent + degree_axis * 2^-51 * Mbound`.

Thus the sign decisions used to discard boxes are exact rational comparisons; floating-point arithmetic supplies only dyadic centers inside rigorously propagated envelopes.

## Files

- `src/jnr_near_aligned_rational_interval_certificate_v3.py` - symbolic kernel and exact-envelope subdivision engine.
- `scripts/verify_exact_v3.py` - one-command replay and frozen-output comparison.
- `certificates/delta_002/result.json` - exact frozen output for `delta=1/50`.
- `certificates/delta_005/result.json` - exact frozen output for `delta=1/20`.
- `requirements.txt`, `environment.yml` - pinned dependencies/reference environment.
- `RUNTIME.md` - commands and measured reference runtime.
- `PROVENANCE.md` - snapshot provenance and scope.
- `SHA256SUMS` - archive-file hashes.
