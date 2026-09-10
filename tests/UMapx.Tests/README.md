# UMapx mathematical audit tests

The expanded audit covers every source area and inventories all 426 C# source
files. See the [expanded report](../../docs/math-audit-expanded-2026-09-10.md) for
counterexamples, causes, repair priorities, conventions, and remaining gaps.

At source commit `0a1e916` (version 7.5.1.5), the final run contains **8,408 cases:
7,535 passed, 873 failed, none skipped**. Two complete runs produced the same
counts. The original 840-case suite remains included. Production algorithms have
not been changed. Failing tests are enabled and expect the mathematical answer;
the test command intentionally exits with status 1 while defects remain.

Execution coverage is **82.56% of lines** and **74.42% of branches**. These figures
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
[the run summary](../../docs/audit/summary.json).

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
- 1,653 special-function reference cases, 1,991 distribution reference cases,
  and 16 high-precision Hankel matrix fixtures. Fixtures use mpmath 1.3.0 at
  40 or 60 decimal digits, with library inputs rounded to binary32 first.
- Independent PDF integration for distribution moments and entropy. Median/mode
  consistency checks supplement these references; they are not independent
  proofs when they call the library's own CDF/PDF.
- Window formulas, wavelet analysis/synthesis coefficients, impulses, filter
  transfer polynomials, image pixel equations and neutral/constant invariants.
- Real/complex vectors, rectangular matrices, even/odd sizes, singular inputs,
  large finite values, and selected parameter boundaries.

Random test inputs use fixed seeds. Several production algorithms have internal
random initialization without a public seed, so their tests use stated residual
budgets. A dedicated nonparallel collection protects the global SIMD flag.
Known nonterminating prime/Schur cases run in `UMapx.AuditProbe` with a five-second
process deadline and process-tree termination; they do not hang the test host.
Subprocess execution is not added to the parent coverlet coverage totals.

## Tolerances and conventions

`NumericAssert` uses mixed absolute and relative bounds and rejects nonfinite
results when a finite result is expected. Scalar defaults are `2e-6 + 2e-5*abs(x)`.
Special-function references use `2e-4 + 2e-4*abs(x)`. Distribution density/CDF
references normally use `3e-4 + 3e-4*abs(x)`; integrated moments use `2e-3` for
both terms. Individual tests tighten or relax these explicitly according to the
operation, conditioning, approximation, or pixel quantization. Hankel matrix
entries use `5e-4 + 5e-4*abs(x)`.

These are audit acceptance budgets, not an existing library accuracy guarantee.
Absolute tolerances near zero do not prove relative accuracy. Loosening a budget
requires a documented mathematical reason, not merely a failing test.

The test sources record nonstandard but intentional conventions: Fresnel
integrals use cos(t*t)/sin(t*t), Struve H/L are not Hankel functions, matrix
convolution uses correlation orientation, complex 2D transforms conjugate the
right-side matrix, and some wavelet prototypes do not promise perfect
reconstruction. Geometric entropy uses bits in the existing API. Wrapped-Cauchy
variance is circular variance. These differences are not silently called bugs.

One left-half-plane complex LogGamma record is excluded because the API does not
specify its unwrapped branch. Complex generalized-erf continuation remains an
explicit contract question. Unsupported distribution getters and undefined mode
sets are identified separately; a contract check returning successfully is not
proof that an unimplemented numerical operation works.

## Reference generation

Python and mpmath are not required to run the C# tests: JSON fixtures are embedded.
To regenerate them in a Python environment with mpmath 1.3.0:

```powershell
python -m pip install mpmath==1.3.0
python -X utf8 tests/UMapx.Tests/Data/generate_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_extended_reference.py
python -X utf8 tests/UMapx.Tests/Data/generate_distributions.py
```

The generators retain their explicit domains and exclusion rules. Special-function
fixtures omit poles, non-real answers for real APIs, nonfinite values and finite
magnitudes above `1e35`. Distribution fixtures include selected divergent moments
as explicit Infinity/NaN expectations. See the [mpmath documentation](https://mpmath.org/doc/1.3.0/).

To refresh the checked-in report evidence after an intentional new audit, pass
actual paths from that run:

```powershell
python -X utf8 tools/summarize_audit.py --trx artifacts/math-audit/run/full-audit.trx --coverage artifacts/math-audit/run/RESULT-ID/coverage.cobertura.xml --output docs/audit
```

Replace `RESULT-ID` with the collector's directory name. Update the explanatory
report and these counts together; generated evidence includes source and input
hashes to expose stale snapshots. Do not treat a reduced set of failing tests as
proof that every corresponding root cause has been fixed.
