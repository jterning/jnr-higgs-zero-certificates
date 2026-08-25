from __future__ import annotations
import json, sys, tempfile, time
from pathlib import Path
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from jnr_near_aligned_rational_interval_certificate_v3 import build_common, run_certificate
DROP={'seconds','subdivision_seconds'}
def strip_timing(x):
    if isinstance(x,dict): return {k:strip_timing(v) for k,v in x.items() if k not in DROP}
    if isinstance(x,list): return [strip_timing(v) for v in x]
    return x

def verify(delta,h0,relative):
    t0=time.time(); common=build_common()
    frozen=ROOT/relative
    with tempfile.TemporaryDirectory() as td:
        fresh=run_certificate(delta,h0,Path(td)/'fresh.json',common)
    ref=json.loads(frozen.read_text())
    if strip_timing(fresh)!=strip_timing(ref):
        print('MISMATCH',frozen,file=sys.stderr); return 1
    print(f'PASS {frozen.parent.name}: status={fresh["status"]}, boxes={fresh["visited_boxes"]}, subdivision={fresh["subdivision_seconds"]:.3f}s, wall={time.time()-t0:.3f}s',flush=True)
    return 0
