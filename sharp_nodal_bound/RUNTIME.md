# Reference runtime

Reference environment:

- Python 3.13.5
- SymPy 1.14.0

Command:

```bash
python verify_all.py
```

Reference replay on the container used to prepare this archive:

```text
real 1.78 s
user 2.07 s
sys  0.16 s
```

Runtime is not part of the certificate: correctness is determined by exact
symbolic equality checks and comparison with the frozen JSON outputs.
