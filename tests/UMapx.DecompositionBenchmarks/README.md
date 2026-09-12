# Decomposition performance comparison

This standalone .NET 8 runner compares built UMapx assemblies, including the
constructor/property API of 7.5.1.5 and the static tuple API of 8.0.0.1.
No reference to the current library is compiled into the runner. Each operation
loads one assembly in a separate process; versions never benchmark concurrently.

Build the library versions in **Release**, then run from the repository root
with PowerShell 7:

```powershell
dotnet build tests/UMapx.DecompositionBenchmarks -c Release
./tests/UMapx.DecompositionBenchmarks/Compare.ps1 `
  -PreviousAssembly path/to/previous/UMapx.dll `
  -CurrentAssembly sources/bin/Release/netstandard2.0/UMapx.dll `
  -OutputFile artifacts/decomposition-comparison.jsonl
```

Use a new output filename for each run. `-Sizes 256 -Names QR,Cholesky,SVD`
selects a subset; `-Sizes 32 -Rows 512 -Names QR,SVD` tests a tall matrix.
`Lanczos-Full` selects the explicit second reorthogonalization pass. Check any
error records before interpreting a ratio. Earlier releases can fail to produce
finite factors on otherwise valid inputs.

## What is measured

- Seeded real matrices, with symmetric positive definite inputs where required;
  generation and reconstruction checks are outside the timed region.
- The full public operation: construction and retrieval of the primary factors
  for old versions, or `Decompose` for current versions. A compiled delegate
  adapts both APIs; reflection is outside the timed region. Both adapters create
  a small array of references to retain results.
- At least 250 ms and three calls of warmup, then seven measurement batches
  targeting 80 ms each. The JSON contains the median, all samples, and allocated
  bytes per operation across all threads. Garbage collection during a batch is
  included; forced collections between batches are excluded.
- Tiered compilation and ReadyToRun are disabled in each child process to avoid
  JIT tier transitions. These are controlled comparisons; actual application
  timings depend on runtime settings, CPU, sizes, and input spectra.
- SVD/GSVD have an explicit limit of 100; Power/NMF perform 100 iterations, with
  NMF rank 8. The inputs use the same seed for both versions. Legacy methods
  with internal random initialization retain their original behavior.

GramSchmidt now also produces R; LU/LDU use pivoting; Power additionally computes
a Rayleigh quotient. These are comparisons of public operations, not identical
arithmetic. The default Power iteration count changed from 10 to 100 between
releases; this runner deliberately supplies the same count to both.

Reconstruction checks cover the factorizations implemented in the runner's
`Residual` method and are reported independently of timing. The runner is not
a substitute for the correctness suite. In particular, the old default Lanczos
mode can produce inaccurate factors on the clustered spectrum used here.

## Regression coverage

`DecompositionPerformanceRepairTests` verifies large real kernels, rectangular
factors, independently scaled and singular QZ pencils, complementary GSVD
subspaces, and adaptive Lanczos reorthogonalization. NMF has a deterministic
allocation test: iteration count must not increase workspace allocations.
Wall-clock thresholds are deliberately kept out of unit tests.

Debug builds once again have optimization enabled, matching 7.5.1.5. For
diagnostic builds that require unoptimized stepping, explicitly pass
`-p:Optimize=false`; use Release for production performance comparisons.
