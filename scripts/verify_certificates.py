#!/usr/bin/env python3
"""Stable public replay entry point for jnr-higgs-zero-certificates.

This wrapper intentionally decouples the archival/public command from
versioned internal verification scripts.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TARGETS = {
    "nodal": ROOT / "sharp_nodal_bound" / "verify_all.py",
    "finite-tilt": ROOT / "scripts" / "verify_exact_v3.py",
}


def run_target(name: str) -> None:
    target = TARGETS[name]
    if not target.is_file():
        raise SystemExit(f"missing verifier: {target}")
    subprocess.run([sys.executable, str(target)], cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "target",
        nargs="?",
        choices=("all", "nodal", "finite-tilt"),
        default="all",
        help="certificate family to replay (default: all)",
    )
    args = parser.parse_args()

    names = ("nodal", "finite-tilt") if args.target == "all" else (args.target,)
    for name in names:
        print(f"[verify-certificates] {name}")
        run_target(name)


if __name__ == "__main__":
    main()
