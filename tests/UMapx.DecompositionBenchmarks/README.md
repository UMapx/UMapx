# Decomposition benchmarks

This standalone .NET 8 runner compares Release builds of UMapx. It loads one
assembly per process and has no compiled dependency on the current library.
Use PowerShell 7 and run commands from the repository root.

## Run a comparison

Build both library versions in Release, then build the runner:

```powershell
dotnet build tests/UMapx.DecompositionBenchmarks -c Release
./tests/UMapx.DecompositionBenchmarks/Compare.ps1 `
  -PreviousAssembly path/to/previous/UMapx.dll `
  -CurrentAssembly sources/bin/Release/netstandard2.0/UMapx.dll `
  -OutputFile artifacts/decomposition-comparison.jsonl
```

Use a new output filename for every run. Check JSON `error` records and reported
residuals before interpreting timings. Each operation has a 60-second process
timeout; older versions may fail on inputs accepted by the current library.

| Option | Meaning |
| --- | --- |
| `-Sizes 32,128` | Column counts; matrices are square unless `-Rows` is set. |
| `-Rows 257 -Sizes 17 -Names SVD` | A tall 257-by-17 SVD case. Reverse the dimensions for a wide case. |
| `-Names QR,Cholesky,SVD` | Select algorithms instead of the default set. |
| `-Names Lanczos-Full` | Enable the explicit second reorthogonalization pass. |
| `-Domain real` or `-Domain complex` | Compare static APIs in the selected scalar domain. Also supply a supported `-Names` list. |

Without `-Domain`, the runner compares real inputs and supports both the legacy
constructor/property API and the current static tuple API. With `-Domain`, both
assemblies must expose the static API. Supported names in that mode are `SVD`,
`Polar`, `GSVD`, `Householder`, `QR`, `Schur`, `EVD`, `EVD-SPD`, `QZ` and `GEVD`:

```powershell
./tests/UMapx.DecompositionBenchmarks/Compare.ps1 `
  -PreviousAssembly path/to/previous/UMapx.dll `
  -CurrentAssembly sources/bin/Release/netstandard2.0/UMapx.dll `
  -OutputFile artifacts/complex-comparison.jsonl `
  -Domain complex -Sizes 32,128 -Names SVD,EVD-SPD,EVD,Schur,QZ,GEVD
```

`EVD-SPD` uses an exactly symmetric or Hermitian positive definite input.
Domain mode checks reconstruction or eigenvector residuals and orthogonality
before timing, and reports both errors in the JSON.

## Measurement method

- Both versions receive the same seeded input within a domain. Generation,
  conversion, validation and reflection are outside the timed region. Real and
  complex inputs have different spectra, so their timing ratio is not a direct
  measure of scalar arithmetic cost.
- The timed operation includes construction and retrieval of primary factors
  for legacy APIs, or `Decompose` for static APIs. Compiled adapters retain the
  results in a small reference array included in the allocation measurement.
- Warmup lasts at least 250 ms and three calls. Seven subsequent batches target
  80 ms each. JSON records the median time, individual samples and allocated
  bytes per call across all threads. Allocation is not peak live memory.
- Garbage collection during a batch is included; forced collection between
  batches is excluded. Versions run sequentially, with tiered compilation and
  ReadyToRun disabled to avoid JIT tier transitions.
- SVD, Polar and GSVD use an iteration limit of 100. In the default real mode,
  Power and NMF use 100 iterations and NMF uses rank `min(columns, 8)`. Legacy
  random initialization is left unchanged.

API differences can change the work being measured: current GramSchmidt also
returns R, LU/LDU use pivoting, and Power also computes a Rayleigh quotient.
Default-mode residuals cover only the operations implemented in `Residual`;
a missing residual is not a correctness check. Use the
[test suite](../README.md) for broader validation.

Timings depend on the runtime, machine, matrix shape and spectrum. Inspect sample
variation and repeat uncertain comparisons before claiming a speedup or
regression. Keep raw results under the ignored `artifacts` directory.
