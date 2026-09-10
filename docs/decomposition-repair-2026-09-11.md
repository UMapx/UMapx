# Matrix decomposition repairs — September 11, 2026

**B06 is closed: all 15 previously failing decomposition cases now pass.**
The complete Windows Release run contains **15,203 cases: 15,035 passed,
168 failed, none skipped**. All 14,707 previous test IDs remain. The comparison
records 15 resolved failures, no regressions, and **496 added passing tests**.
The remaining failures belong to B07–B12.

| Original B06 family | Resolved cases |
| --- | ---: |
| Zero-matrix Schur termination | 3 |
| SVD reconstruction and pseudoinverse | 9 |
| Generalized eigenvectors | 3 |
| **Total** | **15** |

## Changes

- [Schur.cs](../sources/Decomposition/Schur.cs): exact-zero subdiagonals deflate
  even when the norm is zero. QR iteration has a finite limit and reports
  nonconvergence. The initial QR step no longer tests a stale zero shift as if
  it were the current reflector norm. Hessenberg reduction and QR accumulation
  use scaled double work buffers. Converged real pairs and entries below the
  first subdiagonal are stored with their structural zeros.
- [SVD.cs](../sources/Decomposition/SVD.cs): the cancellation rotation now updates
  the auxiliary diagonal and every row of U, including row zero. Splitting
  stops when the canceled element is negligible. The routine checks convergence
  after the last permitted sweep and throws if unfinished. Scaled double work
  buffers avoid overflow/underflow in intermediate shifts and retain isolated
  small singular values alongside large ones.
- `SVD.P` uses the numerical-rank threshold
  `max(rows, columns) * 2^-23 * max(S)`. Singular values below or equal to this
  threshold contribute zero; retained components are accumulated in double
  precision. Zero and rank-deficient matrices therefore have finite
  pseudoinverses when the mathematical result fits the public float type.
- [GEVD.cs](../sources/Decomposition/GEVD.cs): right transformations start from
  the identity. QZ stages use scaled double work buffers and respect the
  convergence status before back-substitution. The tolerance has a binary64
  roundoff floor. Alpha and beta are restored to the input scale, while complex
  eigenvector columns retain their joint normalization.
- The shared QZ entry point also uses the corrected work buffers and removes
  the saved epsilon from the bottom-left scratch entry before returning T.
  Both QZ reconstruction equations and strict triangular structure are tested.

All helpers stay in their primary class files and have English XML comments.
The public signatures, defaults, and single-precision output types are unchanged.
Empty/nonfinite matrices, null inputs, invalid SVD iteration limits, and NaN
deflation tolerances now produce explicit argument errors. Failed iterations
produce `InvalidOperationException` instead of hanging or returning partial factors.

The cancellation and QZ bookkeeping were checked against the original
[EISPACK SVD](https://www.netlib.org/eispack/svd.f),
[QZHES](https://www.netlib.org/eispack/qzhes.f),
[QZIT](https://www.netlib.org/eispack/qzit.f), and
[QZVEC](https://www.netlib.org/eispack/qzvec.f) routines.

## Validation and limits

[DecompositionRepairTests.cs](../tests/UMapx.Tests/DecompositionRepairTests.cs)
adds 496 cases covering:

- Rectangular SVDs, zero matrices, rank-one/rank-two inputs, repeated singular
  values, and deterministic dense inputs over scales from `1e-30` to `1e30`.
- Reconstruction, orthogonality, sorted nonnegative singular values, all four
  Penrose equations, and the independent closed-form rank-one pseudoinverse.
- The relative pseudoinverse threshold and explicit iteration-limit behavior.
- Schur zero/diagonal/dense/rank-one matrices, nilpotent Jordan blocks, real
  representations of complex pairs, and zero/default/explicit tolerances.
- Generalized eigen-equations, known diagonal eigenvalues, repeated eigenvalues,
  complex-pair normalization, and common pencil scales from `1e-20` to `1e20`.
- A singular B through homogeneous eigen-equations, and finite alpha/beta for
  input matrices expressed in very different units.
- Preservation of isolated `1e-30` and `1e30` entries in the same decomposition,
  input immutability, invalid-input handling, and QZ scratch-storage cleanup.

Checks use independent double-precision products and relative Frobenius
residuals, without an absolute tolerance that would hide errors in tiny matrices.
The original 218 decomposition tests also pass; existing acceptance tolerances
and subprocess deadlines were retained. The full run checks dependent Polar,
GSVD, QZ, and other algorithms, with no previously passing case regressed.

This closes the tested B06 defects, not every possible input domain. Numerical
rank truncation deliberately changes `SVD.P` near singularity. Outputs remain
float, so results outside its representable range cannot be returned as finite
numbers. Infinite generalized eigenvalues require the homogeneous alpha/beta
representation; a singular pencil need not have a unique eigensystem. Double
work buffers also increase temporary memory use; no performance claim is made.

## Evidence

- [Test-ID comparison](audit-decomposition/comparison.json).
- [Current summary and source/input hashes](audit-decomposition/summary.json).
- [Source inventory](audit-decomposition/source-inventory.md).
- [Remaining failure assignments](audit-decomposition/repair-blocks.json).
- [Previous snapshot](audit-consolidation/summary.json).
- [Updated repair plan](remaining-repair-blocks-2026-09-10.md).

The 427 source files have reported execution coverage of **83.55% of lines**
and **77.07% of branches**. Coverage is not a mathematical correctness guarantee.

```powershell
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/b06-repair/final
```

The full command exits with status 1 because the 168 remaining expectations stay enabled.
