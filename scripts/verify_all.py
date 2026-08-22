#!/usr/bin/env python3
"""Recompute both interval certificates and compare with frozen results.

Each case is run in a fresh Python process to keep the peak memory use bounded.
Exit status 0 means both runs returned CERTIFIED_EMPTY and all invariant
certificate metadata matched the archived JSON files. Timing fields are ignored.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CASES = [
    ("delta_002", ROOT / "scripts" / "run_delta_002.py"),
    ("delta_005", ROOT / "scripts" / "run_delta_005.py"),
]

EXACT_KEYS = [
    "status", "delta", "X_delta", "h0", "outer_domain", "statement",
    "visited_boxes", "prune_counts", "max_axis_depths", "deepest_depth_vector",
    "hashes", "arithmetic_note", "bernstein_note",
]
FLOAT_KEYS = [
    "delta_float", "X_delta_float", "h0_float", "rho_delta_lower_bound",
    "epsilon_delta_lower_bound",
]
DICT_FLOAT_KEYS = [
    "min_certified_Bernstein_sign_margin", "max_coefficient_error_bound_seen",
]
BUILD_EXACT = ["degree", "shape"]
BUILD_FLOAT = ["max_initial_exact_to_float_error", "roundoff_bound_initial", "Mbound"]


def close(a: float, b: float) -> bool:
    return math.isclose(float(a), float(b), rel_tol=2e-12, abs_tol=2e-14)


def compare(case: str, got: dict, ref: dict) -> list[str]:
    errors: list[str] = []
    for k in EXACT_KEYS:
        if got.get(k) != ref.get(k):
            errors.append(f"{case}: exact mismatch in {k}")
    for k in FLOAT_KEYS:
        if not close(got.get(k), ref.get(k)):
            errors.append(f"{case}: float mismatch in {k}: {got.get(k)} != {ref.get(k)}")
    for k in DICT_FLOAT_KEYS:
        gd, rd = got.get(k, {}), ref.get(k, {})
        if set(gd) != set(rd):
            errors.append(f"{case}: key mismatch in {k}")
            continue
        for sub in gd:
            if gd[sub] is None or rd[sub] is None:
                if gd[sub] != rd[sub]:
                    errors.append(f"{case}: mismatch in {k}.{sub}")
            elif not close(gd[sub], rd[sub]):
                errors.append(f"{case}: float mismatch in {k}.{sub}")
    gb, rb = got.get("build", {}), ref.get("build", {})
    if set(gb) != set(rb):
        errors.append(f"{case}: build polynomial set mismatch")
    else:
        for poly in gb:
            for k in BUILD_EXACT:
                if gb[poly].get(k) != rb[poly].get(k):
                    errors.append(f"{case}: mismatch in build.{poly}.{k}")
            for k in BUILD_FLOAT:
                if not close(gb[poly].get(k), rb[poly].get(k)):
                    errors.append(f"{case}: float mismatch in build.{poly}.{k}")
    return errors


def main() -> int:
    all_errors: list[str] = []
    for case, runner in CASES:
        print(f"[verify] recomputing {case} in a fresh process", flush=True)
        subprocess.run([sys.executable, str(runner)], cwd=ROOT, check=True)
        case_dir = ROOT / "certificates" / case
        ref = json.loads((case_dir / "result.json").read_text())
        generated_path = case_dir / "result.generated.json"
        got = json.loads(generated_path.read_text())
        errors = compare(case, got, ref)
        if got.get("status") != "CERTIFIED_EMPTY":
            errors.append(f"{case}: status is not CERTIFIED_EMPTY")
        if errors:
            print(f"[FAIL] {case}")
            for err in errors:
                print(f"  - {err}")
        else:
            print(
                f"[PASS] {case}: CERTIFIED_EMPTY, "
                f"visited_boxes={got['visited_boxes']}, "
                f"epsilon>={got['epsilon_delta_lower_bound']:.15g}"
            )
        generated_path.unlink(missing_ok=True)
        all_errors.extend(errors)

    if all_errors:
        print(f"\nVerification failed with {len(all_errors)} mismatch(es).")
        return 1
    print("\nAll interval certificates reproduced and matched the frozen results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
