#!/usr/bin/env python3

import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERIFIER = ROOT / "jnr_near_aligned_rational_interval_certificate_v3.py"

CASES = [
    (
        "delta002",
        ROOT / "certificates" / "jnr_rational_interval_v3_delta002.json",
    ),
    (
        "delta005",
        ROOT / "certificates" / "jnr_rational_interval_v3_delta005.json",
    ),
]


def load(path):
    with path.open() as f:
        return json.load(f)


def invariant_part(data):
    """Remove machine-dependent timing fields recursively."""
    if isinstance(data, dict):
        return {
            key: invariant_part(value)
            for key, value in data.items()
            if key not in {"seconds", "subdivision_seconds"}
        }
    if isinstance(data, list):
        return [invariant_part(value) for value in data]
    return data


with tempfile.TemporaryDirectory(prefix="jnr_exact_v3_") as td:
    prefix = Path(td) / "audit"

    subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--case",
            "both",
            "--prefix",
            str(prefix),
        ],
        cwd=ROOT,
        check=True,
    )

    ok = True

    for suffix, frozen_path in CASES:
        generated_path = Path(f"{prefix}_{suffix}.json")

        generated = invariant_part(load(generated_path))
        frozen = invariant_part(load(frozen_path))

        if generated != frozen:
            ok = False
            print(f"FAIL: {suffix}", file=sys.stderr)

            keys = sorted(set(generated) | set(frozen))
            for key in keys:
                if generated.get(key) != frozen.get(key):
                    print(
                        f"  {key}:\n"
                        f"    generated = {generated.get(key)!r}\n"
                        f"    frozen    = {frozen.get(key)!r}",
                        file=sys.stderr,
                    )
        else:
            print(f"PASS: {suffix} matches frozen invariant certificate")

    if not ok:
        raise SystemExit(1)

print("PASS: all exact-rational finite-tilt certificates verified")
