# JNR Higgs-zero interval certificates

This repository contains the rigorous computer-assisted interval certificates
used to make the near-aligned charge bounds effective in the charge-two JNR
dyonic-instanton Higgs-zero analysis.

The theorem-level finite-tilt verifier is
`jnr_near_aligned_rational_interval_certificate_v3.py`. It constructs the
initial Bernstein data exactly over the rationals and permits every pruning
step only after a rigorous exact-rational interval sign test. The earlier
certificate engine in `src/` is retained as an independent cross-check.

The certificates use the denominator-free reduced CKLL node kernel in the
variables `(p,a,b,y,u,v)` and multivariate Bernstein bounds on a deliberately
enlarged compact box. The two frozen results are:

| delta | h0 | status | rho lower bound | epsilon lower bound |
|---:|---:|---|---:|---:|
| 0.02 | 24/25 = 0.96 | `CERTIFIED_EMPTY` | 7/25 = 0.28 | 103/7350 ≈ 0.0140136 |
| 0.05 | 9/10 = 0.90 | `CERTIFIED_EMPTY` | sqrt(19)/10 ≈ 0.435890 | 43/1140 ≈ 0.0377193 |

These imply, on the regular JNR stratum, the effective sufficient conditions

```text
j_max/Q_E >= 1.686666...,  j_min/Q_E > 0.9859863946
```

and

```text
j_max/Q_E >= 1.716666...,  j_min/Q_E > 0.9622807018
```

respectively.

## Repository layout

```text
jnr_near_aligned_rational_interval_certificate_v3.py
    authoritative theorem-level finite-tilt verifier
certificates/jnr_rational_interval_v3_delta002.json
certificates/jnr_rational_interval_v3_delta005.json
    frozen exact-rational v3 certificates
scripts/verify_certificates.py
    stable public replay entry point for all paper-level certificates
scripts/verify_exact_v3.py
    internal finite-tilt replay; reruns both v3 certificates and compares
    invariant metadata
docs/FINITE_TILT_EXACT_INTERVAL.md
    exact-rational arithmetic and certificate specification
frozen/v3/environment.yml
    immutable environment specification for the paper snapshot
sharp_nodal_bound/
    exact symbolic checks supporting the sharp nodal bound
src/
  interval_certificate.py   certificate engine and exact polynomial definitions
  polynomial_kernel.py      standalone symbolic kernel for inspection/reuse
certificates/
  delta_002/result.json     frozen delta=0.02 certificate
  delta_005/result.json     frozen delta=0.05 certificate
scripts/
  run_delta_002.py          recompute the delta=0.02 certificate
  run_delta_005.py          recompute the delta=0.05 certificate
  verify_all.py             recompute both and compare with frozen results
docs/
  certificate_method.md     method and arithmetic-assurance notes
environment.yml             archived development Conda environment
requirements.txt            pip dependency pins
CITATION.cff                 citation metadata
.zenodo.json                 Zenodo release metadata
LICENSE                      BSD-3-Clause
SHA256SUMS                   hashes of archived source and certificate files
```

## Environment

The paper snapshot uses:

```text
Python 3.13.5
NumPy  2.3.5
SymPy  1.14.0
```

The immutable Conda/Mamba specification for that snapshot is
`frozen/v3/environment.yml`:

```bash
mamba env create -f frozen/v3/environment.yml
mamba activate jnr-higgs-zero-certificates-v3
```

The top-level `environment.yml` records the archived development environment.
Alternatively, with an existing Python 3.13 environment:

```bash
python -m pip install -r requirements.txt
```

## Verification

The stable public replay command is

```bash
python scripts/verify_certificates.py
```

It replays both the exact symbolic sharp-nodal checks and the exact-rational
finite-tilt certificates. The Make target is an alias for the same public
interface:

```bash
make verify
```

Individual paper-level certificate families can be replayed through the same
stable entry point:

```bash
python scripts/verify_certificates.py nodal
python scripts/verify_certificates.py finite-tilt
```

The finite-tilt replay reruns both exact-rational v3 interval trees and exits
with status 0 only if both return `CERTIFIED_EMPTY` and their
machine-independent certificate metadata agree with the frozen JSON files.
Wall-clock timing is ignored.

The earlier independent finite-tilt implementation is retained and can be
checked with

```bash
make legacy-verify
```

Individual legacy runs remain available with

```bash
make delta002
make delta005
```

## Arithmetic assurance

Initial Bernstein coefficients in the v3 verifier are constructed exactly as
Python `Fraction` objects. Numerical coefficient centers used during
subdivision are interpreted as their exact binary64 dyadic rationals, and each
carries an exact rational enclosure radius. De Casteljau subdivision
propagates a conservative rational roundoff bound. A box is discarded only
when an exact rational comparison proves that the relevant interval is
strictly separated from zero.

For the `delta=1/50` certificate, the smallest certified rational sign margin
is about `5.37778e-6`, while the largest exact rational error radius encountered
is about `4.36664e-10`. For `delta=1/20`, the corresponding values are about
`2.51022e-5` and `4.00284e-10`.

These decimal values are diagnostics only; the pruning decisions themselves
use the exact rational interval endpoints.

See `docs/FINITE_TILT_EXACT_INTERVAL.md` and the frozen v3 JSON files for the
full theorem-level arithmetic specification and certificate metadata.
`docs/certificate_method.md` documents the earlier independent implementation.

## Sharp nodal-bound symbolic verification

Exact symbolic checks supporting the sharp nodal bound
\(\Xi_J \ge 4\) are provided in `sharp_nodal_bound/`.

The stable replay command is

```bash
python scripts/verify_certificates.py nodal
```

It reproduces the nodal elimination identity, resultant factorizations,
Sturm checks, and exact Bernstein sign checks used in the proof.

## Citation and archival release

For a public release, create a tagged GitHub release (for example
`v1.0.0-paper`) and archive that release with Zenodo. `CITATION.cff` and
`.zenodo.json` are included so GitHub and Zenodo can ingest the metadata. Add
the final repository URL and Zenodo DOI to those files after the first public
release if desired.

## License

BSD-3-Clause. See `LICENSE`.
