# Selected imaging repairs — September 11, 2026

The requested items 1, 2, 3, and 5 are repaired: lookup-table normalization,
bitmap stride, histogram median, and ColorTransfer. **All 21 assigned failures
now pass.** The complete Windows Release run has **16,099 cases: 16,071 passed,
28 failed, none skipped**. All 15,622 previous test IDs remain; no passing case
regressed, and all **477 added cases pass**.

| Requested item | Resolved failing cases |
| --- | ---: |
| 1. Correction tables and pixel equations | 8 |
| 2. Padded and negative bitmap stride | 8 |
| 3. Histogram median | 3 |
| 5. ColorTransfer statistics and constant channels | 2 |
| **Total** | **21** |

The remaining B11 failures are tensor alpha (4), depth histogram overflow (2),
constant-series Figure rendering (15), rectangle IoU overflow (1), and depth
Merge coordinates (1): **23 cases**. B12 retains its five video-parser failures.
Their exact test IDs and outcomes are unchanged. Counts refer to test cases,
not independent root causes. B11 remains partially open.

## Correction tables

[Intensity.cs](../sources/Imaging/Intensity.cs) samples its 14 one-dimensional
function tables at `i / (length - 1)`, including both unit-interval endpoints.
The old `i / length` grid disagreed with Correction's multiplication by 255 and
darkened neutral byte values. Empty tables stay empty; a singleton evaluates
the function at zero. Existing local two-dimensional masks and the Quantize
palette were not changed.

Neutral Shift bypasses the logarithm/exponential round trip, and Contrast uses
the algebraically equivalent `x + contrast * (x - 0.5)` to preserve the identity
exactly. Inversion samples the complementary index directly; subtraction of
rounded normalized values could otherwise lose one byte during truncation.
Pixel saturation/truncation conventions are retained.

Tests check all 14 tables against independent double-precision functions at
lengths 0, 1, 2, 3, 17, and 256. Seven neutral filters preserve every possible
byte through five applications. Additional checks cover independent nonneutral
pixel equations and exact byte inversion, including its two-pass identity.

## Signed stride

[RGBFilter.cs](../sources/Imaging/RGBFilter.cs),
[Grayscale.cs](../sources/Imaging/Grayscale.cs),
[TransparencyCorrection.cs](../sources/Imaging/TransparencyCorrection.cs), and
[ErrorDiffusionDithering.cs](../sources/Imaging/ErrorDiffusionDithering.cs) restart
each logical row at `Scan0 + (long)y * Stride`. Pixel loops no longer traverse
padding or assume increasing physical row addresses. Diffusion's neighbor
offsets use the same signed stride as the main raster.

Guarded buffers verify every nonpixel byte as well as expected pixels for both
stride signs, zero/4/20-byte padding, singleton and thin images, and rectangular
images. RGB, grayscale, transparency, and Floyd-Steinberg results are checked
against independent scalar equations. All 11 predefined diffusion kernels also
produce the same pixels for packed and padded/reversed physical row layouts.

## Histogram median

[Statistics.Median](../sources/Imaging/Statistics.cs) uses the one-based rank
`(population + 1) / 2`. This fixes odd populations and retains the lower median
for even populations. Empty populations return zero. Counts must be
nonnegative; negative counts now produce an argument exception.

Totals and cumulative ranks use 64-bit integers. References sort independent
sample lists for odd/even populations and repeated values, and include
histograms whose population exceeds `Int32.MaxValue`.

## ColorTransfer

[ColorTransfer.cs](../sources/Imaging/ColorTransfer.cs) computes the mean and
population standard deviation over every pixel of each color plane using
Welford's centered update in double precision. The old deviation of column
deviations measured spatial variation in contrast instead of overall contrast.

The established `Factor` and `Inverted` behavior is retained. With destination
statistics `mu_t, sigma_t` and reference statistics `mu_s, sigma_s`, output is
`mu_s + gain * (value - mu_t)`, where:

- `Inverted=false`: `gain = (sigma_t / sigma_s) * (1 + Factor)`;
- `Inverted=true`: `gain = (sigma_s / sigma_t) / (1 + Factor)`.

These are the existing complementary modes; their boolean meanings were not
swapped. For the documented factor domain [0, 10], if either deviation is zero,
gain is defined as zero and the output plane becomes the reference mean. This
preserves identical constant images and gives finite results when a constant
channel cannot provide a usable deviation ratio.

The independent RGB reference computes population variance from pairwise
pixel differences: `sum(i < j, (x_i - x_j)^2) / N^2`. Tests include distinct
means/contrasts, constant and nonconstant channel combinations, both modes,
three factors, unequal image sizes, reference reshaping and replication, and
reference-image immutability. Same-image checks also exercise singleton/thin
images and all four supported color spaces, retaining their established
conversion/quantization tolerances.

## Evidence

- [New tests](../tests/UMapx.Tests/ImagingRepairTests.cs): 477 passing cases.
- [Full summary](audit-imaging-selected/summary.json) and
  [test-ID comparison](audit-imaging-selected/comparison.json).
- [Remaining failures](audit-imaging-selected/failures.md) and
  [28 unique block assignments](audit-imaging-selected/repair-blocks.json).
- [Source inventory and hashes](audit-imaging-selected/source-inventory.json).
- Raw final evidence: `artifacts/math-audit/imaging-selected-repair/final-documented/`.
  TRX and coverage hashes are recorded in the summary; raw artifacts are ignored
  by Git. The preceding snapshot remains under `docs/audit-b07-b10/`.

The focused run passed all 554 selected and related cases. The complete audit
records **83.75% line coverage and 77.45% branch coverage**, including failing
tests. These are execution metrics, not a mathematical correctness guarantee.
Original tests and tolerances were retained. No public signatures changed.
All new numerical helpers live in primary class files with English XML comments.

Run the complete audit from the repository root:

```powershell
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/imaging-selected-repair/final-documented
```

The command intentionally exits with status 1 while the remaining 28 failures
stay enabled. The solution, project files, and NuGet configuration were not changed.
