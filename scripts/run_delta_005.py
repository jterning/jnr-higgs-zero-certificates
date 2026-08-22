#!/usr/bin/env python3
from __future__ import annotations
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from interval_certificate import build_common, run_certificate

out = ROOT / "certificates" / "delta_005" / "result.generated.json"
common = build_common()
result = run_certificate(Fraction(1, 20), Fraction(9, 10), out, common)
print(f"{result['status']}: {out}")
