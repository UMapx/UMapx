# UMapx Mathematical Tests

This project was added during the September 10, 2026 audit. It checks mathematical
identities, specific counterexamples, and special function values against an
independent implementation. See the [audit report](../../docs/math-audit-2026-09-10.md)
for details.

Results for source commit `9d623dcdede20411f220dee976e62cf4c8cf5728`:
**840 tests, 615 passed, 225 failed**. The build succeeds; the test exit code is 1.
Failing tests assert the correct mathematical results and reproduce the issues
found during the audit. They are enabled and do not encode the library's incorrect
behavior as the expected result. The number of failing tests is not the number
of independent defects.

## Running the Tests

Run from the repository root with .NET SDK 8 or later and the .NET 8 runtime installed:

```powershell
dotnet test sources/UMapx.sln -c Release -p:GeneratePackageOnBuild=false
```

To run individual categories:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Identity"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Regression"
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "Category=Reference"
```

| Category | Purpose | Passed | Failed |
| --- | --- | ---: | ---: |
| Identity | FFT, matrix decompositions, transform inversion, integration, interpolation, FIR, and distances | 48 | 0 |
| Regression | Minimal counterexamples for the 31 issue groups in the report | 0 | 43 |
| Reference | Special function comparisons against mpmath | 567 | 182 |

The FFT test uses an independent direct DFT implemented with `System.Numerics.Complex`.
Matrix products are computed independently with `double` accumulation.
Random inputs use fixed seeds. Invertibility alone cannot establish that a formula
is correct: matching errors in the forward and inverse methods can cancel out.
The FFT tests therefore also compare the output against the mathematical definition.

## Reference Values

`Data/special-functions.json` contains 750 records computed with mpmath 1.3.0
at 60 decimal digits of precision. All non-integer arguments are first rounded to
IEEE 754 binary32, so the reference implementation receives the same input as UMapx.
The resulting values are stored in JSON as double-precision numbers.

These records produce 749 tests. One complex `LogGamma` case in the left half-plane
is excluded during data enumeration because the API does not specify the continuous
branch of the logarithm of the gamma function. A difference of `2*pi*i` alone is
not treated as a proven defect here. This exclusion is explicit in the test source.

The Fresnel references follow the UMapx convention: integrals of `cos(t*t)` and
`sin(t*t)`. The normalized mpmath functions are adjusted by a change of variable.
`H` and `L` are compared against Struve functions, not Hankel functions.

The reference acceptance criterion is
`abs(actual - expected) <= 2e-4 + 2e-4*abs(expected)`.
For complex values, the error is measured by its magnitude. This is the criterion
chosen for the audit, not a previously published accuracy guarantee from the library.
Near zero, it permits an absolute error of `2e-4` and does not establish a small
relative error. Regression tests use tighter tolerances; the absolute tolerance is
reduced further for the very small value of `Beta(20,30)`. `NaN` and infinite results
always fail when the reference value is finite.

Python is not required to build or run the tests: the JSON file is embedded in the
test assembly as a resource. To optionally regenerate the data in a Python
environment with mpmath 1.3.0:

```powershell
python -m pip install mpmath==1.3.0
python tests/UMapx.Tests/Data/generate_reference.py
```

The generator excludes undefined, infinite, and excessively large reference values
(magnitude greater than `1e35`), as well as non-real results for real-valued APIs.
The data therefore does not cover poles, the entire `float` range, or every side
of complex branch cuts. For details about the independent reference implementation,
see the [mpmath documentation](https://mpmath.org/doc/1.3.0/).
