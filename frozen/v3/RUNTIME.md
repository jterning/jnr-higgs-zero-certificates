# Reference environment and runtime

Reference environment:

- Python 3.13.5
- NumPy 2.3.5
- SymPy 1.14.0

Canonical command from the archive root:

```bash
python scripts/verify_exact_v3.py
```

The two cases can also be replayed separately:

```bash
python scripts/run_delta_002.py
python scripts/run_delta_005.py
```

Timings are machine-dependent.  The frozen JSON files record the theorem-level
polynomial-build and subdivision times.  Fresh standalone replays in the
reference environment on 2026-08-24 gave:

- `delta=1/50`: `CERTIFIED_EMPTY`, 37,583 visited boxes; 20.116 s subdivision,
  26.316 s verifier wall time, 28.12 s process wall time.
- `delta=1/20`: `CERTIFIED_EMPTY`, 30,159 visited boxes; 16.228 s subdivision,
  22.525 s verifier wall time, 24.17 s process wall time.

The exact certificate metadata compared by the replay deliberately excludes
timing fields.
