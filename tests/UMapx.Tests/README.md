# UMapx tests

The xUnit suite exercises numerical algorithms, transforms, imaging, rendering,
video parsing, and public API contracts. It includes independent reference
values, mathematical identities, and regression cases for previously observed
failures.

## Run the tests

Requirements for the complete suite:

- Windows, because bitmap and rendering tests use System.Drawing.Common.
- .NET SDK 8 or later and the .NET 8 runtime.
- NuGet access for the initial restore, or already restored dependencies.

Run commands from the repository root:

```powershell
dotnet test UMapx.sln -c Release -p:GeneratePackageOnBuild=false
```

The command builds the test project, the probe executable, and the library. Package generation is
disabled because testing does not require a NuGet package. Add `--no-restore`
only when dependencies have already been restored.

To collect coverage and a TRX results file:

```powershell
dotnet test UMapx.sln -c Release -p:GeneratePackageOnBuild=false -p:DebugType=portable -p:DebugSymbols=true --collect "XPlat Code Coverage" --settings tests/UMapx.Tests/coverage.runsettings --logger "trx;LogFileName=tests.trx" --results-directory artifacts/tests/coverage
```

Both PDB settings are required: the library's normal Release configuration
disables debug symbols. [coverage.runsettings](coverage.runsettings) collects
Cobertura and JSON coverage for the UMapx assembly, excluding the test assembly
and generated `obj` files. Results are written to `artifacts/tests/coverage`;
generated artifacts are not tracked in Git.

Coverage measures executed code, not numerical correctness. A passing run
validates the included cases; it does not establish that every supported input
is correct. Contract tests for unsupported operations do not establish that
those operations are implemented.

To select an area, a test class, or a specific regression:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Decomposition"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~TimeoutStreamRepairTests"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~IsolatedLargeEigenvalueDoesNotCorruptSmallBlockEigenvectors"
```

For reproducible decomposition timings across library versions, use the
[standalone comparison runner](../UMapx.DecompositionBenchmarks/README.md).
`DecompositionPerformanceRepairTests` covers numerical accuracy of the optimized
real kernels, QZ accumulation, Lanczos reorthogonalization, and NMF workspace
reuse. Wall-clock performance thresholds are kept outside the unit suite.

Categories describe subject areas, not operating-system compatibility.
For example, [ApproximationRepairTests.cs](ApproximationRepairTests.cs) contains
bitmap tests under `Category=Analysis`. Excluding `Imaging`, `Geometry`,
`Video`, and `Contract` therefore does not produce a guaranteed portable suite.
`SupportedOSPlatform` attributes document restrictions; they do not
automatically skip tests.

## Projects and coverage areas

[UMapx.Tests.csproj](UMapx.Tests.csproj) targets .NET 8 and references the
`netstandard2.0` library. [UMapx.AuditProbe](../UMapx.AuditProbe/Program.cs)
is a separate executable for operations whose termination must be bounded.

| Categories | Main checks |
| --- | --- |
| `Core`, `Matrix`, `Analysis`, `Distance` | Scalar and complex arithmetic, number theory, containers, matrix operations, approximation, interpolation, calculus, and distances. |
| `Reference`, `Identity`, `Regression` | Special-function reference values, mathematical identities, branch conventions, and boundary regressions. |
| `Decomposition` | Real and complex factorizations, reconstruction, orthogonality, eigenvector equations, pseudoinverses, rank, scaling, and convergence contracts. |
| `Distribution` | Densities, CDFs, support, moments, entropy, medians, modes, and parameter boundaries. |
| `Transform`, `Wavelet`, `Window`, `WindowTransform`, `Response` | Direct transform references, reconstruction, window formulas, wavelet coefficients, framed transforms, and filter responses. |
| `ColorSpace`, `Imaging`, `Geometry` | Color conversions, pixel equations, bitmap composition, stride and padding, depth maps, tensors, and rendering. |
| `Video` | MIME boundaries, partial reads, JPEG framing, synthetic video sources, stream deadlines, and exception propagation. |
| `Contract` | Public API behavior, invalid inputs, and explicitly unsupported operations. |

Video tests use synthetic streams and generated images. They do not exercise
live cameras, external MJPEG servers, or screen capture.

## Numerical checks and reproducibility

Independent checks include direct DFT formulas, double-precision matrix
products, all four Penrose equations, BigInteger arithmetic, exact index
mappings, high-precision fixtures, scalar pixel formulas, guarded image buffers,
and numerical integration of distribution densities. Checks that reuse the
library's own CDF or PDF establish consistency rather than independent accuracy.

[NumericAssert.cs](NumericAssert.cs) uses a mixed absolute and relative bound.
The scalar default is `2e-6 + 2e-5 * abs(expected)`; finite expectations reject
nonfinite results. Tests override these tolerances according to the operation,
conditioning, and quantization. Some matrix tests use relative residuals rather
than scalar tolerances.

[SpecialFunctionRepairTests.cs](SpecialFunctionRepairTests.cs) and
[DistributionRepairTests.cs](DistributionRepairTests.cs) tighten absolute
bounds for small nonzero results so that returning zero cannot pass merely
because of a large absolute tolerance. Undefined or infinite expected values
have explicit checks. Test tolerances are acceptance criteria for those cases,
not a general accuracy guarantee. Do not relax a tolerance just to hide a
failure.

Branch choices, normalization, matrix orientation, and parameter conventions
are recorded beside the relevant assertions and in the fixture generators.
Use those sources when extending a test, especially for complex functions,
wavelets, statistical distributions, and image interpolation.

Random test inputs use fixed seeds where generated by the suite. Some
production algorithms initialize randomly without exposing a seed; their tests
use residual bounds. [TestCulture.cs](TestCulture.cs) sets invariant culture,
and the `SIMD audit` collection in
[MatrixFilterAuditTests.cs](MatrixFilterAuditTests.cs) disables parallel
execution while changing the global SIMD setting.

[AuditProcess.cs](AuditProcess.cs) runs selected primality, factorization,
Schur, and eigenvalue termination checks in the probe process with a
five-second deadline and process-tree termination. Child-process execution is
not included in the parent test host's coverage totals.

## Reference data

The JSON files in [Data](Data) are embedded resources. Python and mpmath are
not needed to run the C# tests.

To regenerate fixtures, use Python with mpmath 1.3.0 and run the required script
from the repository root:

```powershell
python -m pip install mpmath==1.3.0
python -X utf8 tests/UMapx.Tests/Data/generate_reference.py
```

| Generator in `tests/UMapx.Tests/Data` | Output files | Decimal precision |
| --- | --- | --- |
| [generate_reference.py](Data/generate_reference.py) | `special-functions.json` | 60 |
| [generate_extended_reference.py](Data/generate_extended_reference.py) | `special-functions-extended.json`, `hankel.json` | 60 |
| [generate_special_repair_reference.py](Data/generate_special_repair_reference.py) | `special-functions-repair.json` | 80 |
| [generate_arithmetic_reference.py](Data/generate_arithmetic_reference.py) | `arithmetic-repair.json` | 100 |
| [generate_distributions.py](Data/generate_distributions.py) | `distributions.json`, `distribution-catalog.json` | 40 |
| [generate_distribution_repair_reference.py](Data/generate_distribution_repair_reference.py) | `distribution-repair.json` | 70 |
| [generate_meyer_hankel_references.py](Data/generate_meyer_hankel_references.py) | `meyer-hankel.json` | 70 |

Generators specify their input rounding, parameter domains, exclusions, and
handling of poles or divergent results. Preserve those conventions when
regenerating data. Inspect fixture changes and run the affected tests before
accepting regenerated values.
