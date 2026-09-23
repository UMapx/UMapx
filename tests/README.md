# UMapx tests

The xUnit suite covers numerical algorithms, matrix decompositions, transforms,
imaging, rendering and public API contracts. It uses independent
reference values, mathematical identities and regression cases.

## Run the tests

The complete suite requires Windows, .NET SDK 8 or later, and the .NET 8 runtime.
Bitmap and rendering tests use System.Drawing.Common. The first restore requires
NuGet access unless the dependencies are already cached.

Run from the repository root:

```powershell
dotnet test UMapx.sln -c Release -p:GeneratePackageOnBuild=false
```

This builds the library and test project without producing a
NuGet package. Add `--no-restore` when dependencies have already been restored.
The test output reports the current case count and results.

To select an area or a test class:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Decomposition"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~NumberTheoryRepairTests"
```

Categories group subjects, not operating systems. Bitmap tests also occur under
`Category=Analysis`; excluding imaging categories does not make the remaining
suite portable. `SupportedOSPlatform` attributes do not automatically skip tests.

## Coverage

```powershell
dotnet test UMapx.sln -c Release -p:GeneratePackageOnBuild=false -p:DebugType=portable -p:DebugSymbols=true --collect "XPlat Code Coverage" --settings tests/UMapx.Tests/coverage.runsettings --logger "trx;LogFileName=tests.trx" --results-directory artifacts/tests/coverage
```

Both PDB settings are needed because the library's Release configuration disables
debug symbols. [coverage.runsettings](UMapx.Tests/coverage.runsettings) collects
Cobertura and JSON reports for UMapx, excluding test assemblies and generated
`obj` files. Coverage records executed code; numerical accuracy is checked by
the assertions and reference data.

## Test conventions

- [NumericAssert.cs](UMapx.Tests/NumericAssert.cs) uses an absolute and relative
  tolerance: `2e-6 + 2e-5 * abs(expected)` by default for scalars. Finite
  expectations reject nonfinite results. Individual tests select bounds for the
  operation, conditioning and quantization; matrix tests also use relative
  residuals. Preserve the documented tolerances when extending tests.
- Decomposition checks cover reconstruction, orthogonality, eigenvector
  equations, pseudoinverses, singular pencils, scale separation and convergence
  limits. Compare identities and spectra when factors have nonunique signs,
  phases or bases for repeated eigenvalues.
- Random inputs use fixed seeds. [TestCulture.cs](UMapx.Tests/TestCulture.cs)
  selects invariant culture. Tests changing the global SIMD setting use the
  nonparallel `SIMD audit` collection.
- Imaging tests use generated bitmaps and check lock ownership, disposal and
  recovery from failures.

## Reference data

JSON files in [Data](UMapx.Tests/Data) are embedded resources. Python is needed
only to regenerate them, not to build or run the C# tests.

For special-function, distribution and arithmetic fixtures, install
`mpmath==1.3.0` and run the required generator from the repository root:

```powershell
python -m pip install mpmath==1.3.0
python -X utf8 tests/UMapx.Tests/Data/generate_reference.py
```

| Generator in `tests/UMapx.Tests/Data` | Output files | Decimal precision |
| --- | --- | --- |
| [generate_reference.py](UMapx.Tests/Data/generate_reference.py) | `special-functions.json` | 60 |
| [generate_extended_reference.py](UMapx.Tests/Data/generate_extended_reference.py) | `special-functions-extended.json`, `hankel.json` | 60 |
| [generate_special_repair_reference.py](UMapx.Tests/Data/generate_special_repair_reference.py) | `special-functions-repair.json` | 80 |
| [generate_arithmetic_reference.py](UMapx.Tests/Data/generate_arithmetic_reference.py) | `arithmetic-repair.json` | 100 |
| [generate_distributions.py](UMapx.Tests/Data/generate_distributions.py) | `distributions.json`, `distribution-catalog.json` | 40 |
| [generate_distribution_repair_reference.py](UMapx.Tests/Data/generate_distribution_repair_reference.py) | `distribution-repair.json` | 70 |
| [generate_meyer_hankel_references.py](UMapx.Tests/Data/generate_meyer_hankel_references.py) | `meyer-hankel.json` | 70 |

For decomposition references, use NumPy (the stored fixtures were generated with
version 2.3.5):

```powershell
python -m pip install numpy==2.3.5
python tests/reference/generate_decomposition_unification.py
```

[The decomposition generator](reference/generate_decomposition_unification.py)
evaluates exact float32 inputs in complex128 using NumPy/LAPACK and writes
[decomposition-unification.json](UMapx.Tests/Data/decomposition-unification.json).
GEVD references use `eigvals(solve(B,A))`
only for well-conditioned, nonsingular B. Singular pencils are tested separately
using homogeneous eigenvalue equations.

Generators define rounding, branch and parameter conventions. Preserve those
conventions, inspect fixture changes and run the affected tests before accepting
regenerated values.
