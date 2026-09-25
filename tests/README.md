# UMapx tests

Run commands from the repository root:

```shell
dotnet test tests/UMapx.Tests.csproj -c Debug
dotnet test tests/UMapx.Tests.csproj -c Release
dotnet test tests/UMapx.Tests.csproj --filter "Category=Geometry"
dotnet test tests/UMapx.Tests.csproj --collect:"XPlat Code Coverage" --settings tests/coverage.runsettings
```

## Organization

Tests live in `Tests/`, grouped by library namespace. Each directory name matches
the corresponding `UMapx` namespace.

| Directory | Scope | Category |
| --- | --- | --- |
| `Analysis` | Integration, differentiation, roots, interpolation and approximation | `Analysis` |
| `Colorspace` | Color conversions and component preservation | `Colorspace` |
| `Core` | Scalar arithmetic, number theory, special functions, kernels, heaps, console diagnostics, geometry and matrices | `Core`, `Geometry`, `Matrix` |
| `Decomposition` | Matrix factors, eigenvalues, numerical edge cases and independent references | `Decomposition` |
| `Distance` | Real, complex and boolean distances | `Distance` |
| `Distribution` | Probability references, moments, modes, medians and time-frequency kernels | `Distribution` |
| `Response` | FIR and IIR responses and stability | `Response` |
| `Transform` | Direct, fast, windowed and multichannel transforms and their filters | `Transform` |
| `Wavelet` | Filter banks, reconstruction and analytic wavelets | `Wavelet` |
| `Window` | Window formulas, dimensions and boundaries | `Window` |

Name test classes after their subject and test methods after the behavior they verify.
Keep regression cases with the relevant subject. Shared helpers and their tests live
directly in `Tests/`; test classes should not provide utilities for other test classes.
Tests of shared helpers use the `Infrastructure` category.

`GeometryRoundingTests` holds shared rounding rules and exceptional conversions.
Type-specific tests retain randomized component checks. Hash tests check equality
and dictionary behavior; unequal values are allowed to have the same hash.

Tests that change console output or SIMD settings use nonparallel collections and
restore the previous state. `TestCulture` initializes a consistent default culture;
formatting tests explicitly select and restore their own culture.

## Reference data

`Data/*.json` files are embedded in the test assembly. Normal test runs read the
checked-in values and do not require Python. The adjacent generators use `mpmath`
for high-precision references or NumPy/LAPACK for decomposition references.

| Generator | Output |
| --- | --- |
| `generate_arithmetic_reference.py` | `arithmetic.json` |
| `generate_decompositions.py` | `decomposition.json` |
| `generate_distributions.py` | `distributions.json`, `distribution-catalog.json` |
| `generate_distribution_edge_cases.py` | `distribution-edge-cases.json` |
| `generate_special_functions.py` | `special-functions.json` |
| `generate_extended_references.py` | `special-functions-extended.json`, `hankel.json` |
| `generate_special_function_edge_cases.py` | `special-functions-edge-cases.json` |
| `generate_meyer_hankel_references.py` | `meyer-hankel.json` |

General reference sets and edge-case sets retain their own accuracy thresholds.
Small probabilities and representable tails require relative accuracy; merging
fixtures must not replace these checks with a larger absolute tolerance.
