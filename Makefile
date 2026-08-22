.PHONY: verify delta002 delta005 hashes

verify:
	python scripts/verify_all.py

delta002:
	python scripts/run_delta_002.py

delta005:
	python scripts/run_delta_005.py

hashes:
	sha256sum src/interval_certificate.py src/polynomial_kernel.py \
	  certificates/delta_002/result.json certificates/delta_005/result.json \
	  scripts/verify_all.py scripts/run_delta_002.py scripts/run_delta_005.py \
	  README.md docs/certificate_method.md environment.yml requirements.txt \
	  CITATION.cff .zenodo.json LICENSE > SHA256SUMS
