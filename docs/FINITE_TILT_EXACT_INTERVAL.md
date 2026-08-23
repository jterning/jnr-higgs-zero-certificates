# Finite-tilt exact-rational interval certificate

This verifier supplies the rigorous finite-tilt exclusion used in the near-aligned charge theorem.

## Certified feasibility problem

It excludes simultaneous solutions on the enlarged outer box

- `p in [0,1/2]`,
- `a,b,u,v in [-1,1]`,
- `y in [0,1]`,

of

```text
Ktilde = 0,
Cu     = 0,
Cv     = 0,
Cap   <= 0,
Nh    >= 0.
```

`Cap` encodes `Xi_J <= X_delta`; `Nh` encodes `h >= h0`.  The search box is larger than the physical Schur/JNR domain, so certified emptiness of this polynomial problem implies the physical exclusion used in the manuscript.

## Arithmetic model

1. Monomial coefficients and root-box tensor-product Bernstein coefficients are constructed exactly over `Q` using `fractions.Fraction`.
2. Binary64 numbers are used only as de Casteljau center values.  Every stored center is interpreted as its exact dyadic rational with `Fraction.from_float`.
3. Every polynomial state carries an exact rational uniform radius `E`.  After a half-box split in coordinate `i`,

   ```text
   E_child = E_parent + d_i * 2^(-51) * B,
   ```

   where `d_i` is the Bernstein degree in that coordinate and `B` is an exact power-of-two bound on all stored centers.
4. The center-bound invariant and `E < 1` are checked after every split.
5. All pruning signs are decided by exact rational comparisons of interval endpoints with zero.  No floating-point inequality certifies a sign.

## Reproduce

From the repository root:

```bash
python jnr_near_aligned_rational_interval_certificate_v3.py \
  --case 0.02 --prefix audit002

python jnr_near_aligned_rational_interval_certificate_v3.py \
  --case 0.05 --prefix audit005
```

or run both:

```bash
python jnr_near_aligned_rational_interval_certificate_v3.py \
  --case both --prefix audit
```

The script requires Python 3, NumPy, and SymPy.

## Frozen certified outputs

### `delta = 1/50`, `h0 = 24/25`

- `X_delta = 403/103`
- status: `CERTIFIED_EMPTY`
- visited boxes: `37583`
- maximum depths `(p,a,b,y,u,v) = (7,7,7,6,6,6)`
- theorem width: `epsilon_delta >= 103/7350`
- smallest exact certified pruning gap:
  `3875097268924783/720575940379279360000`
- largest exact interval radius:
  `28807/65970697666560`

### `delta = 1/20`, `h0 = 9/10`

- `X_delta = 163/43`
- status: `CERTIFIED_EMPTY`
- visited boxes: `30159`
- maximum depths `(p,a,b,y,u,v) = (6,6,6,6,6,5)`
- theorem width: `epsilon_delta >= 43/1140`
- smallest exact certified pruning gap:
  `361761209966267/14411518807585587200`
- largest exact interval radius:
  `26407/65970697666560`

The JSON files under `certificates/` retain the exact rational quantities, polynomial hashes, pruning counts, build metadata, and subdivision-depth data.

## Repository policy

This verifier supersedes the older binary64-only finite-tilt certificate as the theorem-level implementation.  Retain the previous implementation only as historical/legacy code if desired; do not cite it as the proof certificate for the current manuscript.
