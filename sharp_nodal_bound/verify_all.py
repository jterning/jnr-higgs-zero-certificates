#!/usr/bin/env python3
"""Run all exact symbolic checks and compare with frozen JSON outputs."""
from __future__ import annotations
import json
from pathlib import Path
from verify_sharp_nodal_bound import compute as compute_sharp
from verify_isosceles_sturm import compute as compute_sturm

ROOT = Path(__file__).resolve().parent
CASES = [
    ("outputs/sharp_nodal_bound.json", compute_sharp),
    ("outputs/isosceles_sturm.json", compute_sturm),
]

for rel, fn in CASES:
    got = fn()
    expected = json.loads((ROOT/rel).read_text())
    if got != expected:
        raise SystemExit(f"MISMATCH: {rel}")
    print(f"PASS {rel}")
print("ALL SYMBOLIC CHECKS PASS")
