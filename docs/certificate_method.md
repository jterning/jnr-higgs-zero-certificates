# Certificate method

The certificates rule out ordinary rank-loss nodes in an enlarged compact box
using the denominator-free reduced CKLL nodal kernel in variables
`(p,a,b,y,u,v)`.

The simultaneous feasibility conditions are

- `Ktilde = 0`,
- `Cu = 0`,
- `Cv = 0`,
- `Cap <= 0`, which contains the spin-gap cap `Xi_J <= X_delta`, and
- `Nh >= 0`, which is necessary for a node with `h >= h0`.

The outer domain is deliberately larger than the physical Schur/JNR domain:

```text
p in [0, 1/2]
a,b,u,v in [-1,1]
y in [0,1]
```

Each polynomial is converted to a tensor-product Bernstein representation on
this box. Initial Bernstein coefficients are constructed exactly as rational
numbers. They are rounded once to binary64 and the exact conversion error is
recorded. During half-box de Casteljau subdivision, each coefficient tensor
carries a conservative absolute roundoff envelope enlarged by
`degree_axis * 2^-51 * Mbound` per split. A box is pruned only when every
exact Bernstein coefficient is certified to have the sign needed to exclude
one of the feasibility conditions.

The frozen results are:

- `delta = 0.02`, `h0 = 24/25`: `CERTIFIED_EMPTY`, yielding
  `rho_delta >= 7/25` and `epsilon_delta >= 103/7350`.
- `delta = 0.05`, `h0 = 9/10`: `CERTIFIED_EMPTY`, yielding
  `rho_delta >= sqrt(19)/10` and `epsilon_delta >= 43/1140`.

The archived JSON files contain polynomial hashes, box counts, pruning counts,
depth statistics, strict Bernstein sign margins, and propagated coefficient
error bounds.
