.PHONY: verify finite-tilt nodal exact-v3 legacy-verify delta002 delta005 hashes

# Stable public verification interface: replays all paper-level certificates.
verify:
	python scripts/verify_certificates.py

finite-tilt:
	python scripts/verify_certificates.py finite-tilt

nodal:
	python scripts/verify_certificates.py nodal

# Internal theorem-level finite-tilt verifier, retained for direct use.
exact-v3:
	python scripts/verify_exact_v3.py

# Retained independent/legacy finite-tilt implementation.
legacy-verify:
	python scripts/verify_all.py

delta002:
	python scripts/run_delta_002.py

delta005:
	python scripts/run_delta_005.py

hashes:
	sha256sum \
	  jnr_near_aligned_rational_interval_certificate_v3.py \
	  certificates/jnr_rational_interval_v3_delta002.json \
	  certificates/jnr_rational_interval_v3_delta005.json \
	  scripts/verify_certificates.py \
	  scripts/verify_exact_v3.py \
	  sharp_nodal_bound/verify_all.py \
	  docs/FINITE_TILT_EXACT_INTERVAL.md \
	  src/interval_certificate.py \
	  src/polynomial_kernel.py \
	  certificates/delta_002/result.json \
	  certificates/delta_005/result.json \
	  scripts/verify_all.py \
	  scripts/run_delta_002.py \
	  scripts/run_delta_005.py \
	  README.md \
	  Makefile \
	  docs/certificate_method.md \
	  frozen/v3/environment.yml \
	  environment.yml \
	  requirements.txt \
	  CITATION.cff \
	  .zenodo.json \
	  LICENSE > SHA256SUMS
