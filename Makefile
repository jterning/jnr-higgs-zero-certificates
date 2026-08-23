.PHONY: verify exact-v3 legacy-verify delta002 delta005 hashes

# Canonical theorem-level verification.
verify: exact-v3

exact-v3:
	python scripts/verify_exact_v3.py

# Retained independent/legacy implementation.
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
	  scripts/verify_exact_v3.py \
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
	  environment.yml \
	  requirements.txt \
	  CITATION.cff \
	  .zenodo.json \
	  LICENSE > SHA256SUMS
