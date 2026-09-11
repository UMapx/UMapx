# UMapx mathematical audit tests

The audit covers every source area. The current source inventory contains 427 C#
files after [helper consolidation](../../docs/helper-consolidation-2026-09-11.md).
See the [B07–B10 repair report](../../docs/b07-b10-repair-2026-09-11.md)
and [remaining repair blocks](../../docs/remaining-repair-blocks-2026-09-10.md).
The [decomposition repair report](../../docs/decomposition-repair-2026-09-11.md),
[approximation repair report](../../docs/approximation-repair-2026-09-11.md),
[matrix and distribution repair report](../../docs/matrix-distribution-repair-2026-09-10.md),
[arithmetic repair report](../../docs/arithmetic-repair-2026-09-10.md),
[special-function repair report](../../docs/special-functions-repair-2026-09-10.md),
and [expanded baseline report](../../docs/math-audit-expanded-2026-09-10.md) remain historical records.

The current complete run contains **15,622 cases: 15,573 passed, 49 failed, none
skipped**. All 15,203 previous cases remain: 119 B07–B10 failures now pass and no
passing cases regressed. All 419 added cases pass. B01–B10 are closed; B11 (44)
and B12 (5) remain open. Failing tests remain enabled and expect the mathematical
answer; the full command exits with status 1.

Execution coverage is **83.69% of lines** and **77.40% of branches**. These figures
include failing and contract tests. They are not a correctness percentage, and
this suite does not establish absence of errors. The report explicitly lists
unexecuted lines/methods, unsupported APIs, and incomplete parameter domains.

## Run the complete audit

Requirements: Windows, .NET SDK 8 or later, the .NET 8 runtime, and restored NuGet
packages. Bitmap and rendering tests use System.Drawing.Common and require
Windows. SupportedOSPlatform attributes document that restriction but do not
skip tests automatically on other platforms. No camera, live network stream, or
screen capture is used.

From the repository root:

```powershell
dotnet test sources/UMapx.sln -c Release -p:GeneratePackageOnBuild=false
```

For coverage and a compact evidence export, also install Python 3 (standard
library only) and run:

```powershell
./tools/Run-MathAudit.ps1
# After restoring dependencies:
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/run-two
```

The script preserves the failing test exit code. It creates TRX, Cobertura and
coverlet JSON results, then generates source inventory, failures, and a summary
under the selected results directory. It does not overwrite the report snapshot
in `docs/audit`. Both portable PDB settings are essential because the library's
normal Release configuration disables debug symbols.

Equivalent coverage command without Python:

```powershell
dotnet test sources/UMapx.sln -c Release -p:GeneratePackageOnBuild=false -p:DebugType=portable -p:DebugSymbols=true --collect "XPlat Code Coverage" --settings tests/UMapx.Tests/coverage.runsettings --logger "trx;LogFileName=full-audit.trx" --results-directory artifacts/math-audit/run
```

For a particular mathematical area or counterexample:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Decomposition"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~MeshEvaluationUsesBothIndependentCoordinates"
```

Available categories: `Identity`, `Regression`, `Reference`, `Core`, `Matrix`,
`Analysis`, `ColorSpace`, `Decomposition`, `Distance`, `Distribution`, `Window`,
`WindowTransform`, `Transform`, `Wavelet`, `Response`, `Imaging`, `Geometry`,
`Video`, and `Contract`. Category totals and test-family counts are available in
[the run summary](../../docs/audit-b07-b10/summary.json).

A numeric-only filter for environments without Windows bitmap support is:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category!=Imaging&Category!=Geometry&Category!=Video&Category!=Contract"
```

This filter has not been validated as a full cross-platform compatibility test.
The complete reported result is from Windows with .NET SDK 10.0.401.

## Independent evidence

- Direct DFT and transform matrices, double-precision matrix products,
  reconstruction residuals, orthogonality, and all four Penrose equations.
- BigInteger arithmetic, independent scalar and complex equations, exact index
  mappings, and component/stride tests with guarded memory.
- 1,029 arithmetic references at 100 decimal digits, boundary cases spanning the
  float range, cubic/quadratic residuals and Vieta identities, exact integer
  arithmetic at signed limits, pseudoprimes, and large-semiprime factorization.
- 3,727 special-function reference cases, 1,991 distribution reference cases,
  and 16 high-precision Hankel matrix fixtures. Fixtures use mpmath 1.3.0 at
  40, 60, or 80 decimal digits, with library inputs rounded to binary32 first.
  The dedicated special-function suites also contain 22 boundary/identity cases.
- Independent PDF integration for distribution moments and entropy. Median/mode
  consistency checks supplement these references; they are not independent
  proofs when they call the library's own CDF/PDF.
- 1,611 additional distribution references at 70 decimal digits, including small
  tails, moment-existence boundaries, direct entropy/moment integrals, and exact
  discrete probability sums. A further 55 checks cover modes, scaling and support.
- 316 matrix/distance repair cases cover clipping, empty/rectangular arrays,
  extreme shifts, mixed diagonal products, Hermitian statistics and contingency counts.
- 595 approximation/filter repair cases cover Padé Taylor equations, clamped
  grid interpolation, direct window sums, sorted morphology windows, independent
  cubic Hermite resampling, exact orthogonal rotations, and bitmap channels.
- 496 decomposition repair cases cover matrix rank and scaling, all four Penrose
  equations, Schur deflation and complex blocks, generalized eigen-equations,
  homogeneous eigenvalues, QZ structure, and explicit convergence/input contracts.
- 419 B07–B10 repair cases cover biorthogonal impulses and vanishing moments,
  independent Meyer quadrature, narrow windows, high-order Bessel zeros, known
  IIR poles, complex grid sample influences, direct local-Laplacian remapping,
  threshold components, cone kernel mass, and neutral/chromatic color round trips.
- Window formulas, wavelet analysis/synthesis coefficients, impulses, filter
  transfer polynomials, image pixel equations and neutral/constant invariants.
- Real/complex vectors, rectangular matrices, even/odd sizes, singular inputs,
  large finite values, and selected parameter boundaries.

Random test inputs use fixed seeds. Several production algorithms have internal
random initialization without a public seed, so their tests use stated residual
budgets. A dedicated nonparallel collection protects the global SIMD flag.
Prime termination regressions, large factorization cases, and the repaired Schur zero-matrix termination cases run in `UMapx.AuditProbe` with a five-second
process deadline and process-tree termination; they do not hang the test host.
Subprocess execution is not added to the parent coverlet coverage totals.

## Tolerances and conventions

`NumericAssert` uses mixed absolute and relative bounds and rejects nonfinite
results when a finite result is expected. Scalar defaults are `2e-6 + 2e-5*abs(x)`.
Original special-function references use `2e-4 + 2e-4*abs(x)`. New special-function
references use `1e-7 + 2e-5*abs(x)`, replacing the absolute term with `1e-44` for
nonzero reference magnitudes below `1e-4` to test small tails. Distribution density/CDF
references normally use `3e-4 + 3e-4*abs(x)`; integrated moments use `2e-3` for
both terms. Individual tests tighten or relax these explicitly according to the
operation, conditioning, approximation, or pixel quantization. Hankel matrix
entries use `5e-4 + 5e-4*abs(x)`.

The distribution repair references use `2*float.Epsilon + 2e-5*abs(x)`, with
exact discrete medians and explicit NaN/infinity checks. This keeps small tails
under a relative accuracy budget. Existing reference tolerances were retained.

These are audit acceptance budgets, not an existing library accuracy guarantee.
Absolute tolerances near zero do not prove relative accuracy. Loosening a budget
requires a documented mathematical reason, not merely a failing test.

The test sources record nonstandard but intentional conventions: Fresnel
integrals use cos(t*t)/sin(t*t), Struve H/L are not Hankel functions, matrix
convolution uses correlation orientation, complex 2D transforms conjugate the
right-side matrix, and some wavelet prototypes do not promise perfect
reconstruction. Geometric entropy uses bits in the existing API. Wrapped-Cauchy
variance is circular variance. These differences are not silently called bugs.

Complex LogGamma uses analytic continuation with a negative-real-axis cut; the
previously excluded record now runs. Gerf uses the entire continuation of
`n!/sqrt(pi) * integral(exp(-t^n), 0, x)`. Unsupported distribution getters and undefined mode
sets are identified separately; a contract check returning successfully is not
proof that an unimplemented numerical operation works.

Complex sample statistics use squared magnitudes and Hermitian covariance.
PowerNormal and PowerLognormal follow the NIST survival-power laws. Discrete
median getters select the lower median; ChiSquare median uses CDF inversion.
See the current repair report for compatibility details and numerical limits.

Local mean parameters are window lengths. Windows are clipped and renormalized;
even lengths have one extra sample on the left. Weighted matrix means retain
separate horizontal and vertical passes. Bicubic resizing aligns sample centers
using `(index+0.5)*sourceLength/destinationLength-0.5` and floors the anchor.

SVD pseudoinversion discards singular values at or below
`max(rows, columns) * 2^-23 * max(S)`. Schur and GEVD apply a binary64 roundoff
floor to the requested relative tolerance. Their work buffers use scaled double
arithmetic; public outputs remain float. Iteration failure is explicit.

B07–B10 preserve existing public signatures. Biorthogonal convolution accumulates
in double; bands remain float. Complex bilateral grids use magnitude guidance
with shared weights for both components. Complex Under/Over thresholding is
componentwise; Abs compares magnitude. Local-Laplacian processing preserves its
base level, interpolates lookup tables, and remains unsupported for complex data.
XYZ retains nonnegative relative tristimulus values above one.

## Reference generation

Python and mpmath are not required to run the C# tests: JSON fixtures are embedded.
To regenerate them in a Python environment with mpmath 1.3.0:

```powershell
python -m pip install mpmath==1.3.0
python -X utf8 tests/UMapx.Tests/Data/generate_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_extended_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_special_repair_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_arithmetic_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_distributions.py
python -X utf8 tests/UMapx.Tests/Data/generate_distribution_repair_reference.py
```

The generators retain their explicit domains and exclusion rules. Original special-function
fixtures omit poles, non-real answers for real APIs, nonfinite values and finite
magnitudes above `1e35`. The repair generator permits magnitudes through `3e38`
for single-precision results and `1e300` for its selected double-returning APIs.
Distribution fixtures include selected divergent moments
as explicit Infinity/NaN expectations. Arithmetic fixtures omit the three singular
reciprocal-atan inputs (0, +i, -i) and the two reciprocal-hyperbolic zero poles;
zero conventions and real poles have separate tests. See the [mpmath documentation](https://mpmath.org/doc/1.3.0/).

To refresh the checked-in report evidence after an intentional new audit, pass
actual paths from that run:

```powershell
python -X utf8 tools/summarize_audit.py --trx artifacts/math-audit/run/full-audit.trx --coverage artifacts/math-audit/run/RESULT-ID/coverage.cobertura.xml --output docs/audit
```

Replace `RESULT-ID` with the collector's directory name. Update the explanatory
report and these counts together; generated evidence includes source and input
hashes to expose stale snapshots. Do not treat a reduced set of failing tests as
proof that every corresponding root cause has been fixed.

The [B07–B10 reference generator](../../tools/generate_b07_b10_references.py)
requires mpmath 1.3.0 and evaluates at 70 decimal digits. Run it from the repository
root with mpmath on Python's import path. It generates the committed
`Data/b07-b10-repair.json`: 38 Meyer values from independent spectral quadrature
and 24 Hankel matrices from ordered positive Bessel zeros. Tests need no Python
or network access. New Meyer bounds are `2e-7 + 2e-7*abs(x)`; new Hankel bounds
are `3e-5 + 3e-5*abs(x)`. Direct local-Laplacian references use an absolute 3e-4
lookup-interpolation budget for the tested widths. Earlier tolerances are unchanged.
