"""Regenerate double-precision LAPACK references for the exact float32 inputs.

Requires NumPy. Run from the repository root; the test suite only reads the JSON.
"""
import json
from pathlib import Path
import numpy as np

rng = np.random.default_rng(20260919)
cases = []


def rounded(a, real):
    return a.real.astype(np.float32).astype(np.complex128) if real else a.astype(np.complex64).astype(np.complex128)


def packed(a):
    return [[float(z.real), float(z.imag)] for z in a.ravel()]


for real in (True, False):
    for scale in (1e-30, 1.0, 1e30):
        for m, n in ((1, 7), (7, 1), (5, 8), (8, 5), (12, 12)):
            a = rounded(scale * (rng.normal(size=(m, n)) + 1j * rng.normal(size=(m, n))), real)
            values = np.linalg.svd(a, compute_uv=False)
            cases.append(dict(name=f'svd-{real}-{m}x{n}-{scale}', operation='SVD', real=real,
                              rows=m, columns=n, a=packed(a), spectrum=packed(values)))
        for operation in ('EVD', 'Hermitian', 'GEVD'):
            n = 9
            a = rounded(scale * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))), real)
            if operation == 'Hermitian':
                a = rounded((a + a.conj().T) / 2, real)
                values = np.linalg.eigh(a)[0]
            else:
                values = np.linalg.eigvals(a)
            case = dict(name=f'{operation}-{real}-{scale}', operation=operation, real=real,
                        rows=n, columns=n, a=packed(a))
            if operation == 'GEVD':
                # Independent A/B units; B is well conditioned so solving is a suitable reference here.
                b = rounded((np.eye(n) * n + rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / scale, real)
                values = np.linalg.eigvals(np.linalg.solve(b, a))
                case['b'] = packed(b)
            case['spectrum'] = packed(values)
            cases.append(case)

target = Path('tests/UMapx.Tests/Data/decomposition-unification.json')
target.write_text(json.dumps(dict(generator=f'NumPy {np.__version__}; numpy.linalg (LAPACK), complex128 on exact float32 inputs', cases=cases), indent=2) + '\n')
print(f'Wrote {len(cases)} cases to {target}')
