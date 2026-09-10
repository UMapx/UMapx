# Approximation, interpolation, and local-filter repairs — September 11, 2026

**B05 is closed: all 15 assigned failures now pass.** The complete Windows
Release audit contains **14,707 cases: 14,524 passed, 183 failed, none skipped**.
All 14,112 previous test IDs remain. There are **no regressions, no removed tests,
and 595 added passing cases**. Existing assertions, reference data, and tolerances
were retained.

| Scope | Previous failures | Current failures | Resolved |
| --- | ---: | ---: | ---: |
| Padé coefficient systems | 2 | 0 | 2 |
| Bilinear grid interpolation | 2 | 0 | 2 |
| Bicubic array/bitmap resizing | 3 | 0 | 3 |
| Local averaging | 2 | 0 | 2 |
| Morphology ranks | 4 | 0 | 4 |
| Matrix half-turn rotations | 2 | 0 | 2 |
| Other blocks | 183 | 183 | 0 |
| **Total** | **198** | **183** | **15** |

Counts refer to test cases, not independent defects. The resolved ID set is exactly
the B05 assignment set in the previous snapshot. B06–B12 counts did not change.
Previously repaired arithmetic, special functions, matrix operations and
distributions retain their passing results.

## Changes and mathematical contracts

- [Padé](../sources/Analysis/Pade.cs): coefficients at negative Taylor indices
  are zero when constructing the denominator system. Both real and complex
  approximants support denominator degree greater than numerator degree.
  The existing positive-degree restriction and linear solver are retained.
- [Bilinear interpolation](../sources/Analysis/Interpolation.cs): clamp each
  coordinate to the grid and interpolate within its boundary cell. A query on
  an x edge still interpolates along y; y queries outside the grid no longer
  extrapolate. Differences and intermediate interpolation use double precision.
- [Array resizing](../sources/Core/Matrice.cs) and
  [bitmap resizing](../sources/Imaging/Resize.cs): bicubic source coordinates are
  `(destinationIndex + 0.5) * sourceLength / destinationLength - 0.5`.
  The integer anchor uses floor, including negative border coordinates. The
  existing cubic kernel and replicated border samples are retained. This fixes
  the half-sample translation at equal sizes and the truncated cubic support
  near an enlarged image's boundary. Vector and complex overloads use the same
  coordinate convention. Bilinear and nearest-neighbor resize conventions are unchanged.
- [Local mean filters](../sources/Core/LinealgOptions.cs): one rolling-window
  implementation replaces inconsistent index updates in all twelve real/complex,
  vector/matrix, weighted/unweighted axis filters. Samples enter and leave once,
  so work is linear in the line length. Accumulation uses double-precision complex
  arithmetic. Weighted normalization uses the actual weight sum, without the
  former fixed `1e-8` bias that distorted small weights.
- [Morphology sorting](../sources/Core/LinealgOptions.cs): consume the existing
  zero-based rank directly. The median and dilation no longer select the element
  one position below their correct rank. Replicated boundary windows and the
  existing window-size interpretation are retained.
- [Matrix rotation coefficients](../sources/Core/Matrice.Interpolation.cs):
  reduce angles modulo 360 and return exact sine/cosine values at multiples of
  90 degrees. All three real/complex interpolation modes share the coefficients.
  Exact lattice coordinates no longer fall below an integer because of trig
  rounding, which previously broke nearest-neighbor half turns.

Public signatures and defaults are unchanged. Numerical outputs at the repaired
boundaries intentionally change. The averaging parameters denote **window
lengths**, despite the former "radius" wording. A length `r` covers indices
`i-floor(r/2)` through `i+floor((r-1)/2)`, clipped to available samples and
renormalized. Even windows have one extra sample on the left. Lengths below two
and axes shorter than two retain the existing identity behavior.

The public weighted matrix mean retains its **separable horizontal-then-vertical
operation**, using the supplied weights in each pass. It is not a single weighted
sum over a rectangular neighborhood. A zero total weight returns zero in an
active filtering pass, preserving the existing all-zero-weight result.

## Validation

[ApproximationRepairTests](../tests/UMapx.Tests/ApproximationRepairTests.cs) adds
**595 passing cases**:

| Independent check | Cases |
| --- | ---: |
| Padé Taylor-product equations and exponential coefficient ratios, real/complex degrees up to `[4/5]` | 40 |
| Mixed-linear surfaces on nonuniform grids, all edges/corners, clamping and large finite coordinate differences | 73 |
| Direct local window sums, odd/even/oversized windows, empty/singleton vectors, small weights | 140 |
| Independent separable matrix means, real/complex, weighted/unweighted, rectangular/empty shapes | 96 |
| Sorted replicated morphology windows, duplicate values, median/erosion/dilation | 45 |
| Cubic Hermite reference interpolation for arrays/vectors, singleton dimensions, enlargement and reduction | 30 |
| Exact orthogonal rotation permutations, negative angles, odd/even sizes, all interpolation modes | 168 |
| Independent ARGB channel interpolation for bitmap resizing | 3 |

The cubic reference uses a Hermite polynomial with centered endpoint slopes;
it does not call the library's cubic kernel or resize routine. Local-mean
references sum each window directly rather than using a rolling update.
Matrix and complex expectations accumulate in `System.Numerics.Complex` double
precision. Bitmap comparisons allow one channel level for byte quantization.
No existing tolerance was relaxed and no failing expectation was disabled.

## Evidence and reproduction

```powershell
dotnet test tests/UMapx.Tests -c Release --no-restore -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~ApproximationRepairTests"
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/b05-repair/verification
```

The recorded complete run is `artifacts/math-audit/b05-repair/final-confirmed`.
It exits with **1** because all 183 remaining failures stay enabled.

- [Run summary and source/input hashes](audit-approximation/summary.json).
- [Every-test-ID comparison with the previous 14,112-case run](audit-approximation/comparison.json).
- [Remaining failures](audit-approximation/failures.json), [failing families](audit-approximation/failures.md), and [block ownership](audit-approximation/repair-blocks.json).
- [Source inventory and coverage](audit-approximation/source-inventory.md) and [uncovered methods](audit-approximation/uncovered-methods.json).
- [Previous block assignments](audit-matrix-distribution/repair-blocks.json).

Execution coverage is **83.26% of lines** (27,806/33,393) and **76.49% of branches**
(10,529/13,764). The line total decreased because duplicated filtering loops were
replaced by a shared implementation. Coverage is not a correctness percentage.
All 437 production source hashes match the recorded snapshot. Source, tests,
tools and report text contain no Cyrillic characters.

## Remaining work and limits

The remaining blocks are B06 15, B07 57, B08 18, B09 23, B10 21, B11 44 and B12 5.
**B06 is next: matrix decompositions**, starting with the zero-matrix Schur
termination defect; keep its subprocess timeout enabled.

These tests cover selected finite inputs and the repaired algorithms, not all
floating-point values. Padé conditioning and singular systems retain the existing
solver behavior. General-angle rotation sampling and nonfinite filtering inputs
are not exhaustively audited by this repair. Bicubic reduction retains the
existing interpolation kernel without adding an antialiasing filter.
