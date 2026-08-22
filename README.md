# JNR Higgs-zero interval certificates

This repository contains the rigorous computer-assisted interval certificates
used to make the near-aligned charge bounds effective in the charge-two JNR
dyonic-instanton Higgs-zero analysis.

The certificates use the denominator-free reduced CKLL node kernel in the
variables `(p,a,b,y,u,v)` and multivariate Bernstein bounds on a deliberately
enlarged compact box. The two frozen results are:

| delta | h0 | status | rho lower bound | epsilon lower bound |
|---:|---:|---|---:|---:|
| 0.02 | 24/25 = 0.96 | `CERTIFIED_EMPTY` | 7/25 = 0.28 | 103/7350 ≈ 0.0140136 |
| 0.05 | 9/10 = 0.90 | `CERTIFIED_EMPTY` | sqrt(19)/10 ≈ 0.435890 | 43/1140 ≈ 0.0377193 |

These imply, on the regular JNR chart, the effective sufficient conditions

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
environment.yml             exact Conda environment used for the archived run
requirements.txt            pip dependency pins
CITATION.cff                 citation metadata
.zenodo.json                 Zenodo release metadata
LICENSE                      BSD-3-Clause
SHA256SUMS                   hashes of archived source and certificate files
```

## Environment

The archived run used:

```text
Python 3.13.5
NumPy  2.3.5
SymPy  1.14.0
```

With Conda/Mamba:

```bash
mamba env create -f environment.yml
mamba activate jnr-higgs-zero-certificates
```

or with an existing Python 3.13 environment:

```bash
python -m pip install -r requirements.txt
```

## One-command verification

From the repository root:

```bash
python scripts/verify_all.py
```

This recomputes both interval trees and exits with status 0 only if both runs
return `CERTIFIED_EMPTY` and their invariant certificate metadata match the
frozen JSON files. Timing fields are intentionally ignored. On the archived
machine the two subdivisions take tens of seconds in total.

Individual runs are also available:

```bash
python scripts/run_delta_002.py
python scripts/run_delta_005.py
```

## Arithmetic assurance

Initial Bernstein coefficients are constructed exactly as Python `Fraction`
objects. Each is rounded once to binary64 and its exact conversion error is
computed. Every de Casteljau subdivision propagates a conservative absolute
roundoff envelope using `2^-51`, four times the standard binary64 unit
roundoff. A box is discarded only when every exact Bernstein coefficient of a
relevant polynomial is separated from zero by more than that envelope.

For the archived runs, the smallest strict pruning margins are about
`5.38e-6` and `2.51e-5`, while the largest propagated coefficient-error bounds
are about `3.30e-10` and `3.02e-10`, respectively.

See `docs/certificate_method.md` and the frozen JSON files for the full
certificate metadata.

## Citation and archival release

For a public release, create a tagged GitHub release (for example
`v1.0.0-paper`) and archive that release with Zenodo. `CITATION.cff` and
`.zenodo.json` are included so GitHub and Zenodo can ingest the metadata. Add
the final repository URL and Zenodo DOI to those files after the first public
release if desired.

## License

BSD-3-Clause. See `LICENSE`.
