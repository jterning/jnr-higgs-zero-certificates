#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys, time
ROOT=Path(__file__).resolve().parents[1]
t0=time.time()
for name in ('run_delta_002.py','run_delta_005.py'):
    p=subprocess.run([sys.executable,str(ROOT/'scripts'/name)],cwd=ROOT)
    if p.returncode:
        raise SystemExit(p.returncode)
print(f'ALL CERTIFICATES PASS; wall_seconds={time.time()-t0:.3f}',flush=True)
