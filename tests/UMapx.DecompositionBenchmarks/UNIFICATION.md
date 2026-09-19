# Real and complex decomposition algorithms

The real implementations are the algorithmic baseline. Real and complex kernels
retain their own scalar arithmetic and storage; shared numerical iterations
operate through column-transformation callbacks, outside elementwise loops.

## Implementation

| Area | Implementation after unification |
| --- | --- |
| SVD | The original real Golub–Kahan QR iteration is extracted into `DiagonalizeBidiagonal`, used by both domains. Complex input is reduced by Householder transformations, with diagonal phases making the bidiagonal problem real. Reflectors reuse the private input buffer, and only economy factors are constructed. |
| Symmetric/Hermitian EVD | The original real `tql2` algorithm supplies the shared `DiagonalizeTridiagonal` QL iteration. The complex Hermitian path uses rank-two Householder reduction and phase normalization to a real symmetric tridiagonal problem. Both domains share the disconnected-component partition and ascending spectrum ordering. |
| Householder | The complex tridiagonal reduction uses the Hermitian counterpart of the real symmetric rank-two update. Its matrix-vector work buffer is reused. |
| Schur / general complex EVD | Complex Schur uses implicit single-shift QR with direct accumulation of two-row rotations, following the real Hessenberg/implicit-QR strategy. The real double-shift kernel and its 2-by-2 blocks remain specialized. General complex EVD obtains eigenvectors by triangular back-substitution. |
| QZ / GEVD | Complex Hessenberg-triangular reduction follows the real `qzhes` stages, applying each reflector directly to both inputs. Both QZ domains and real/complex GEVD normalize A and B independently and restore their separate units. Existing real single/double-shift and complex single-shift iteration/back-substitution specializations are retained. |
| Basis completion | Complex completion uses the real helper's first sufficiently large candidate and positive-norm failure criterion. |
| Polar / GSVD | Both use the shared SVD iteration, including GSVD's complementary-subspace refinement. |

Other decomposition algorithms already follow the same main scheme in both
domains. This change does not route real matrices through complex storage and
does not attempt to force real quasi-triangular Schur/QZ factors into complex
triangular output types.

## Observable behavior

- Public method signatures and default argument values remain unchanged.
- Complex SVD and GSVD `iterations` now count maximum QR sweeps per singular
  value, matching the real algorithm. Polar forwards that same limit to SVD.
  Previously these complex paths counted cyclic Jacobi sweeps. An explicitly
  small limit can therefore behave differently.
- Exactly Hermitian EVD inputs use QL and return ascending real eigenvalues.
  Near-Hermitian inputs continue to use the general path, preserving small
  imaginary eigenvalues.
- Signs, phases, bases in repeated/zero subspaces, eigenvalue order in general
  problems, and homogeneous alpha/beta normalizations may differ. Validation
  uses mathematical identities and spectra rather than elementwise equality of
  nonunique factors.
- LDL/UDL remain unpivoted, and NMF remains real-only.

## Validation

- Full Release test suite: **18,174 passed, zero failed or skipped**.
- This includes 48 independent NumPy/LAPACK reference cases generated from the
  exact float32 input entries, evaluated in complex128. Cases cover real and
  complex SVD, general/Hermitian EVD, and GEVD with independently scaled inputs.
- Additional tests cover tall/wide SVD, zero singular subspaces, Moore–Penrose
  reconstruction, repeated and signed Hermitian spectra, interleaved disconnected
  blocks with very different scales, singular pencils/infinite eigenvalues, and
  real/complex basis completion at numerical breakdown.
- Existing decomposition regressions cover input preservation, finite outputs,
  structural contracts, unitary factors, tiny imaginary eigenvalues, scale
  separation, convergence limits, and independent concurrent calls.

References are stored in `tests/UMapx.Tests/Data/decomposition-unification.json`.
Regenerate from the repository root with NumPy installed:

```text
python tests/reference/generate_decomposition_unification.py
```

NumPy is only needed to regenerate fixtures, not to build or run the .NET tests.
GEVD references use `eigvals(solve(B,A))` only for well-conditioned nonsingular B;
singular pencils are checked separately against homogeneous equations.

## Performance

Measured on Windows x64, Intel Family 6 Model 183, .NET 8.0.31, Release,
with tiered compilation and ReadyToRun disabled. Each version runs in its own
process; input generation and validation are outside timing. Values are medians
of seven batches following warmup. Allocations are total bytes per public call,
including the small common benchmark wrapper; they are not peak live memory.
Baseline is the working tree before this change (saved Release assembly).

These measurements describe these seeded matrices on this machine. They do not
establish universal speedups or worst-case convergence bounds. Small timing
differences should be treated as measurement variability.

### Complex square matrices

| Method | Size | Before, ms | After, ms | Before/after | Bytes before | Bytes after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SVD | 32 | 3.881 | 0.729 | 5.32x | 66,216 | 78,296 |
| Householder | 32 | 0.127 | 0.127 | 1.01x | 57,992 | 58,528 |
| EVD-SPD | 32 | 6.266 | 0.480 | 13.07x | 1,946,776 | 51,848 |
| EVD | 32 | 9.170 | 1.785 | 5.14x | 2,998,640 | 100,856 |
| Schur | 32 | 9.249 | 1.662 | 5.56x | 2,955,776 | 57,992 |
| QZ | 32 | 4.164 | 4.195 | 0.99x | 166,664 | 107,856 |
| GEVD | 32 | 4.308 | 4.216 | 1.02x | 176,784 | 117,976 |
| GSVD | 32 | 4.634 | 1.548 | 2.99x | 472,744 | 493,064 |
| QR | 32 | 0.097 | 0.097 | 1.00x | 117,336 | 117,336 |
| SVD | 128 | 320.186 | 42.481 | 7.54x | 1,050,408 | 1,196,120 |
| Householder | 128 | 7.162 | 6.796 | 1.05x | 919,688 | 921,760 |
| EVD-SPD | 128 | 2435.638 | 28.207 | 86.35x | 104,657,904 | 795,368 |
| EVD | 128 | 3710.467 | 106.049 | 34.99x | 159,201,616 | 1,582,328 |
| Schur | 128 | 3682.377 | 99.439 | 37.03x | 158,542,128 | 919,688 |
| QZ | 128 | 249.313 | 252.356 | 0.99x | 2,631,176 | 1,708,368 |
| GEVD | 128 | 257.765 | 253.622 | 1.02x | 2,769,936 | 1,847,128 |
| GSVD | 128 | 361.023 | 83.812 | 4.31x | 7,490,440 | 7,741,928 |
| QR | 128 | 5.152 | 5.135 | 1.00x | 1,844,568 | 1,844,568 |

### Real square matrices

| Method | Size | Before, ms | After, ms | Before/after | Bytes before | Bytes after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SVD | 32 | 0.159 | 0.135 | 1.18x | 27,568 | 27,864 |
| Householder | 32 | 0.026 | 0.026 | 1.00x | 36,152 | 36,152 |
| EVD-SPD | 32 | 0.069 | 0.062 | 1.13x | 28,072 | 28,400 |
| EVD | 32 | 0.325 | 0.328 | 0.99x | 27,984 | 27,984 |
| Schur | 32 | 0.255 | 0.257 | 0.99x | 37,976 | 37,976 |
| QZ | 32 | 0.430 | 0.428 | 1.00x | 67,200 | 67,200 |
| GEVD | 32 | 0.408 | 0.400 | 1.02x | 38,680 | 38,680 |
| GSVD | 32 | 0.364 | 0.314 | 1.16x | 251,952 | 252,544 |
| QR | 32 | 0.026 | 0.027 | 0.99x | 46,576 | 46,576 |
| SVD | 128 | 8.687 | 6.900 | 1.26x | 404,272 | 404,568 |
| Householder | 128 | 1.105 | 1.111 | 0.99x | 537,656 | 537,656 |
| EVD-SPD | 128 | 3.706 | 3.056 | 1.21x | 406,024 | 406,352 |
| EVD | 128 | 17.852 | 17.682 | 1.01x | 405,840 | 405,840 |
| Schur | 128 | 11.946 | 11.846 | 1.01x | 544,088 | 544,088 |
| QZ | 128 | 25.236 | 25.223 | 1.00x | 1,005,312 | 1,005,312 |
| GEVD | 128 | 23.861 | 23.690 | 1.01x | 546,712 | 546,712 |
| GSVD | 128 | 18.807 | 15.880 | 1.18x | 3,652,688 | 3,653,280 |
| QR | 128 | 1.364 | 1.354 | 1.01x | 677,104 | 677,104 |

### Rectangular complex SVD

| Shape | Before, ms | After, ms | Before/after | Bytes before | Bytes after |
| --- | ---: | ---: | ---: | ---: | ---: |
| 257 x 17 | 5.667 | 1.247 | 4.54x | 182,200 | 252,000 |
| 17 x 257 | 5.732 | 1.306 | 4.39x | 252,144 | 321,984 |

All 38 before/after scenarios passed the benchmark's validation. The maximum
current relative residual was 6.74e-08; the maximum reported
orthogonality error was 3.88e-08.

Complex SVD has higher total allocations than Jacobi on these cases (about 14%
at 128 x 128 and 28–38% for the rectangular cases), despite reusing the matrix
buffer; temporary reflector vectors account for the remaining tradeoff.
Complex Schur and Hermitian EVD substantially reduce both time and allocations.
QZ and GEVD retain their specialized iteration kernels; their main measured
complex-domain improvement is reduced allocation during the initial reduction.

Raw results for this local run are in
`artifacts/decomposition-unification/results-final.json`. Use `Compare.ps1` to
reproduce comparisons with saved before/after assemblies.

## Subsequent workspace removal

The structural refactoring uses commit `3b735b0` (the completed algorithm
unification) as its baseline. The measurements above describe the earlier
algorithm change; the measurements below describe only workspace removal.

The four nested `RealWorkspace` classes in SVD, EVD, Schur and GEVD have been
removed. Real and complex entry points and factorization drivers sit alongside
each other, followed by their numerical helpers. Temporary arrays belong to
each call and are passed explicitly to static helpers. The shared SVD QR and
symmetric/Hermitian EVD QL iterations are retained. Schur and GEVD write directly
to the result arrays, and symmetric EVD no longer allocates an unused Hessenberg
matrix. No public signatures, defaults, numerical recurrences or ordering rules
changed in this refactoring.

Validation on .NET 8.0.31, Release:

- All **18,174 tests passed**, with no failures or skipped tests.
- A direct comparison against the baseline covered **798 calls** to real and
  complex SVD, EVD, Schur, GEVD, QZ, Polar and GSVD: 658 successful calls produced
  bitwise-identical factors/spectra; 140 exception outcomes matched in type,
  message and parameter name. Inputs were unchanged in both versions.
- Cases included dense, symmetric/Hermitian, diagonal, zero, rank-deficient,
  independently scaled block and rectangular matrices, invalid arguments and
  the 48 existing NumPy/LAPACK fixtures. Public decomposition API signatures,
  parameter names and defaults also matched.
- Source comparison confirmed that 18 numerical kernels were unchanged after
  accounting for movement, comments and explicit local-buffer declarations.

All 18 benchmark scenarios passed reconstruction and orthogonality checks:
six real methods at sizes 32 and 128, and six complex controls at size 128.
The real 128-by-128 results were:

| Method | Before, ms | After, ms | Bytes before | Bytes after |
| --- | ---: | ---: | ---: | ---: |
| SVD | 6.897 | 6.791 | 404,568 | 404,512 |
| EVD-SPD | 3.060 | 3.010 | 406,352 | 271,096 |
| EVD | 17.739 | 13.190 | 405,840 | 405,768 |
| Schur | 11.898 | 11.892 | 544,088 | 404,720 |
| GEVD | 23.812 | 23.719 | 546,712 | 475,928 |
| QZ | 23.892 | 25.163 | 1,005,312 | 1,005,312 |

The general real EVD improvement repeated at 17.742 to 13.215 ms. QZ repeated
at 23.771 to 25.225 ms, but a control comparing the baseline assembly with
itself also varied from 25.200 to 23.821 ms. These short QZ timings therefore
do not isolate a refactoring effect. Complex timings differed by at most 1.4%
and allocations were identical in all six controls. Allocation reductions in
real EVD-SPD, Schur and GEVD were approximately 33%, 26% and 13%, respectively.

Local comparison inputs, audit runner, logs and raw timings are preserved in
`artifacts/decomposition-layout/`. The existing `Compare.ps1` reproduces timing
comparisons using `before.dll` from that directory and the current Release DLL.
