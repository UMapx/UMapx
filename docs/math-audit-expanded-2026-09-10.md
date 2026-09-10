# UMapx expanded mathematical audit — September 10, 2026

This is the historical baseline. The subsequent [special-function repair report](special-functions-repair-2026-09-10.md)
records 423 resolved cases and the current 450 remaining failures.

The audit targets library version **7.5.1.5**, source commit
`0a1e916` (the preceding audit was committed there). No production algorithm was
changed in this expansion. The solution now includes a small subprocess probe
for tests that would otherwise hang the test runner.

**8,408 test cases: 7,535 passed, 873 failed, none skipped.** Two complete runs
produced the same totals and coverage. The build succeeds; the test command
correctly exits with status 1. The failures assert the intended mathematical
results, not the current incorrect output.

This is an audit across the whole source tree, **not a certificate of mathematical
correctness or exhaustive test coverage**. All 426 C# source files are inventoried.
Of 383 files with instrumented executable code, 382 were executed by tests; the
exception is the console diagnostic helper `Core/Debugger.cs`. Running code in a
file does not establish that every method, branch, parameter range, or formula
in it has been verified.

## Evidence and reproducibility

- [Test commands, categories, reference generation, and tolerances](../tests/UMapx.Tests/README.md).
- [Per-file inventory](audit/source-inventory.md), [CSV](audit/source-inventory.csv), and
  [JSON including exact uncovered lines and source hashes](audit/source-inventory.json).
- [All failing test families](audit/failures.md) and
  [individual failing cases, messages, test locations, and IDs](audit/failures.json).
- [Unexecuted methods](audit/uncovered-methods.json), including constructors,
  accessors, compiler-generated methods, and algorithm branches.
- [Machine-readable run summary and evidence hashes](audit/summary.json).
- [Original report: groups 01–31](math-audit-2026-09-10.md). Its 840-case counts
  describe the earlier stage, not this expanded run.

The final raw evidence is under `artifacts/math-audit/final/`: `full-audit.trx`
and a coverlet result subdirectory. Raw artifacts are ignored by Git. The compact
evidence above is checked into the working tree and can be regenerated with
`tools/Run-MathAudit.ps1`. Numeric assertion messages use invariant culture.

### Coverage by source area

Coverage is **28,437 / 34,443 executable lines (82.56%)** and
**9,989 / 13,422 branches (74.42%)**. Coverage includes failing tests and contract
checks; it is an execution metric, not the fraction of mathematics proved correct.
The probe processes are not included in the parent process's coverlet data.

| Area | Source files | Covered / executable lines | Main evidence |
| --- | ---: | ---: | --- |
| Analysis | 15 | 1,021 / 1,116 | Quadrature, finite differences, ODEs, interpolation, fitting, roots, Pade coefficients |
| Colorspace | 17 | 773 / 887 | RGB grids, conversion round trips, independent D65/XYZ/LAB values |
| Core | 27 | 6,708 / 9,059 | Scalar references, exact integer arithmetic, matrix products, statistics, geometry, containers |
| Decomposition | 25 | 2,044 / 2,224 | Reconstruction, orthogonality, structural constraints, generalized eigenproblems, Penrose equations |
| Distance | 19 | 312 / 332 | Independent real/complex distance equations and Boolean contingency counts |
| Distribution | 56 | 1,732 / 2,228 | Independent densities, CDFs, moments, entropy, support, mode/median properties |
| Imaging | 121 | 6,030 / 7,880 | Pixel equations, padding/stride guards, tensor/channel order, neutral settings, image invariants |
| Response | 3 | 230 / 236 | Causal difference equations, transfer polynomials, stability |
| Transform | 42 | 2,842 / 2,970 | Independent transform coefficients, direct/fast agreement, inversion, filter properties |
| Video | 19 | 366 / 757 | In-memory MJPEG parsing, synthetic frames, configuration/lifecycle contracts |
| Visualization | 12 | 575 / 787 | Coordinate maps, curve rendering, constant/singular input handling |
| Wavelet | 37 | 3,868 / 3,934 | Analysis/synthesis banks, impulses, reconstruction, continuous functions |
| Window | 33 | 1,936 / 2,033 | Independent formulas, even/odd sizes, explicit-size overloads |

## Repair register

Groups are work items, **not a count of independent defects**. One cause can fail
many overloads and downstream algorithms. Conversely, a group may require several
formula corrections. The original groups 01–31 all remain reproducible; this
section adds groups 32–88. Reference disagreements and API-contract questions are
explicitly distinguished from localized implementation mistakes.

### Probability distributions

The numeric values below come from independent densities integrated with mpmath,
not from re-evaluating the corresponding UMapx property formula.
`DistributionReferenceTests` contains the full constructor parameters and expected
values in its embedded JSON; `DistributionShapeAuditTests` adds mode/median checks.

| ID | Source | Reproducible finding and repair direction |
| --- | --- | --- |
| 32 | [BetaPrime.Excess](../sources/Distribution/BetaPrime.cs#L157) | `(2,8)` gives `9.2`, expected `13.6`. Correct the fourth central moment / excess formula. |
| 33 | [BirnbaumSaunders](../sources/Distribution/BirnbaumSaunders.cs#L123) | `(0,1,1)` gives median `2.618034` instead of `1`, mode `-1` instead of `0.353209964`, and entropy `1.225791` instead of `1.322062571`. The median is `mu+beta`. The stationary-point cubic has constant **-1**, not +1. The entropy also needs a separate correction. |
| 34 | [FisherSnedecor.Excess](../sources/Distribution/FisherSnedecor.cs#L165) | `(4,12)` gives `3.4`, expected `26.142857143`. Recompute the central fourth moment with the correct degrees of freedom. |
| 35 | [FisherZ.Function](../sources/Distribution/FisherZ.cs#L199) | `(4,12).Function(30)` gives `NaN`; the density is approximately `2.7605e-152`. Zero is an acceptable binary32 underflow result here; `NaN` is not. Evaluate the exponential/power ratio in log space. |
| 36 | [FoldedNormal.ComputeMode](../sources/Distribution/FoldedNormal.cs#L176) | `(.75,1.25)` gives mode `0.683053851`; the correct mode is `0`. For `abs(mu)<=sigma`, the density decreases from zero. The current Lambert-W expression does not solve the stationarity equation. |
| 37 | [Gompertz.Variance](../sources/Distribution/Gompertz.cs#L82) | `(.75,1.25)` gives `0.454489022`, expected `0.149104776`. Correct the second moment, then subtract the squared mean. |
| 38 | [Kumaraswamy.Entropy](../sources/Distribution/Kumaraswamy.cs#L146) | `(2,3)` gives `-2.029632330`, expected `-0.208426136`. Re-derive `-E[log f(X)]`, including the digamma and shape terms. Negative differential entropy itself is valid. |
| 39 | [Levy.Entropy](../sources/Distribution/Levy.cs#L129) | `(.5,1.5)` gives `2.121671200`, expected `3.729947910`. Correct the additive constants and scale term; entropy must not depend on the location parameter. |
| 40 | [PowerNormal](../sources/Distribution/PowerNormal.cs), [PowerLognormal](../sources/Distribution/PowerLognormal.cs) | The code uses `Phi(x)^p`, whereas the NIST definition cited by these classes is `1-Phi(-x)^p`. Decide on the documented law and make both the PDF and CDF consistent with it. This is a distribution-definition mismatch, not a rounding error. |
| 41 | [Trapezoidal moments](../sources/Distribution/Trapezoidal.cs#L111) | `(0,1,2,4)` gives mean `1.5` and variance `3.5`; expected `1.8` and `0.726666667`. Include the middle segment's proper weight in the moments and subtract the squared mean from the second raw moment. |
| 42 | [TukeyLambda moments](../sources/Distribution/TukeyLambda.cs#L50) | `lambda=-2` reports mean `0` although the mean does not exist. At `-.75` and `-2`, variance is not finite, but the implementation returns an invalid finite negative value or `NaN`. Check existence before evaluating beta/gamma formulas; a symmetric principal value is not an expectation. |
| 43 | [Wigner.Entropy](../sources/Distribution/Wigner.cs#L153) | Radius `2` gives `1.644729972`, expected `1.337877066`. The additive difference is `1-log(2)`; correct the entropy constant. |

The failing Gamma, Erlang, ChiSquare, GeneralizedNormal, Nakagami-related special
functions, inverse-gamma laws, and Poisson CDF checks must be re-run after original
groups **01–02** are fixed. Several are consequences of the shared incomplete-gamma
implementation. Beta PDF/CDF failures at `(20,30)` inherit original group **13**.
Poisson mode/entropy checks at rate 100 additionally encounter group **19**.
Generic median tests use the library CDF and therefore detect consistency failures;
they do not by themselves localize the defect to the median getter.

### Arithmetic, matrix operations, distances, and windows

| ID | Source / test suite | Reproducible finding and repair direction |
| --- | --- | --- |
| 44 | [SokalSneath](../sources/Distance/SokalSneath.cs), `DistanceAuditTests` | For mismatch count `m` and shared-true count `tt`, use `2*m/(tt+2*m)`. The implementation uses a different denominator. Three contingency-table cases fail. |
| 45 | [BartlettHann](../sources/Window/BartlettHann.cs), `WindowAuditTests` | The cosine term has the wrong sign. The independent formula is `.62-.48*abs(t-.5)-.38*cos(2*pi*t)`. |
| 46 | [Confined](../sources/Window/Confined.cs), `WindowAuditTests` | `GetWindow(n)` depends on the object's stored frame size even when sigma is held equal. Pass the explicit size through the internal Gaussian helper. This is separate from original group 23's integer centering. |
| 47 | [Maths.IsPrime / Pollard](../sources/Core/Maths.cs#L1178), `NumberTheoryAuditTests` | `IsPrime(1)` fails to terminate for both integer widths. Guard values below 2 before entering factor search. Tests kill the probe process after five seconds. |
| 48 | [Maths.IsPrime / Itf / Radical](../sources/Core/Maths.cs#L1534), `NumberTheoryAuditTests` | `IsPrime(25)` and `IsPrime(169)` return true. A failed Pollard split is not a primality proof. `Itf(n,true)` can return composite factors; `Radical(100)` gives `50` instead of `10`, also corrupting totients. Recursively split and validate prime factors. Default `Itf(n,false)` is allowed to return composite factors by its contract. |
| 49 | [Maths.Gcd / Euclidean](../sources/Core/Maths.cs), `NumberTheoryAuditTests` | `Gcd(18,-6)` returns `-6`. `Euclidean(-3,4)` does not return gcd `1` and valid Bezout coefficients. Normalize the gcd sign and keep quotient/remainder conventions consistent for negative inputs. |
| 50 | [Maths.Lcm](../sources/Core/Maths.cs), `NumberTheoryAuditTests` | `Lcm(50000,50000)` gives `35899`; `Lcm(12345,54321)` gives `223530922` instead of `223530915`. Avoid overflow before division and avoid floating-point arithmetic in exact integer operations. |
| 51 | [Maths.Coprime](../sources/Core/Maths.cs), `NumberTheoryAuditTests` | Searching from 2 for `n=6` returns a non-coprime value. The successful candidate is incremented before being returned. |
| 52 | [Maths.Etf](../sources/Core/Maths.cs), `NumberTheoryAuditTests` | `Etf(3)` gives `1`, expected `2`, even with an unambiguous prime factor. Floating-point multiplication followed by truncation loses one; use exact integer updates. Group 48 is an additional source of totient errors. |
| 53 | [Maths.NumLength](../sources/Core/Maths.cs#L1964), `NumberTheoryAuditTests` | `Decimal2Base(0)` and `ModPow(2,0,7)` throw `OverflowException` via `log(0)`. Zero needs a defined one-digit representation; exponent zero must yield the modular identity. |
| 54 | [Maths.Base2Decimal / Vector2Numeral](../sources/Core/Maths.cs#L1915), `NumberTheoryAuditTests` | A base-conversion round trip changes `123456789` to `123456787`. Accumulate digits in integer arithmetic rather than using binary32 powers and products. |
| 55 | [Maths.Numeral2Vector / Vector2Numeral](../sources/Core/Maths.cs#L1935), `NumberTheoryAuditTests` | `123456789` round-trips as `987654321`. Encoding and decoding use opposite digit orders. Correct the order separately from group 54's precision issue. |
| 67 | [Matrice diagonal Dot overloads](../sources/Core/Matrice.cs), `MatrixAuditTests` | Left diagonal multiplication on a `3x5` matrix indexes the three-element vector by column, leading to incorrect scaling or an exception. Left multiplication must scale row `i` by `v[i]`. Eight mixed overload cases fail. |
| 68 | [Matrice.Shift](../sources/Core/Matrice.cs), `MatrixAuditTests` | Rectangular matrices throw because row/column moduli are swapped. Test both `3x5` and `5x3`, not just squares. |
| 69 | [Matrice.Merge](../sources/Core/Matrice.cs#L8385), [DepthTransform.Merge](../sources/Imaging/DepthTransform.cs), matrix/geometry suites | Inserting a `2x2` block at `(1,1)` does not copy the full region. The loops confuse an extent with an absolute end coordinate. Tests use a constant patch to isolate placement from bicubic interpolation. |
| 70 | [LinealgOptions.MeanFilter](../sources/Core/LinealgOptions.cs), `MatrixAuditTests` | A three-sample local mean of an interior impulse omits a neighbor's expected contribution. Correct initialization and the order of sliding-window removal/addition; both real and complex overloads fail. |
| 71 | [LinealgOptions.MorphologySortFilter](../sources/Core/LinealgOptions.cs#L3004), `MatrixAuditTests` | The rank helper returns a zero-based rank, but the caller subtracts one again. On `[1,2,3,4,5]` with window size 3, the center median is `2` instead of `3`, and dilation is `3` instead of `4`. Fix both vector and matrix callers. |
| 72 | [Matrice.Resize internals](../sources/Core/LinealgOptions.cs), [Imaging.Resize](../sources/Imaging/Resize.cs), `KernelAndContainerAuditTests`, `ImagingAuditTests` | Bicubic resampling to the same dimensions changes the samples. Review pixel-center mapping: the current half-pixel offset is not an identity at scale 1. Nearest-neighbor and bilinear identity cases pass. |
| 73 | [Matrice.Rotate](../sources/Core/LinealgOptions.cs), `MatrixFilterAuditTests` | A 180-degree nearest-neighbor rotation is not exact index reversal. Floating-point coordinates just below an integer are truncated to the previous pixel. Real and complex cases fail. |
| 74 | [Matrice.Compute mesh overloads](../sources/Core/Matrice.cs#L10064), `MatrixStructureAuditTests` | All four real/complex combinations evaluate `f(x[i],y[i])` rather than `f(x[i],y[j])`. Rectangular inputs can also throw. |
| 75 | [Matrice.Swap](../sources/Core/Matrice.cs#L9422), `MatrixStructureAuditTests` | Row swaps iterate over the row count, and column swaps over the column count. Rectangular inputs leave entries untouched or throw; 12 real/complex/direction cases fail. |

### Decompositions, wavelets, and transforms

| ID | Source / test suite | Finding and repair direction |
| --- | --- | --- |
| 56 | [Schur](../sources/Decomposition/Schur.cs), `DecompositionAuditTests` | The zero matrix gives nonfinite factors at size 2 and hangs at sizes 3/5. Handle exact-zero deflation before normalizing or dividing; a relative threshold of zero must still recognize exact zero. |
| 57 | [SVD.P](../sources/Decomposition/SVD.cs#L92), `DecompositionAuditTests` | Zero and rank-deficient matrices produce invalid pseudoinverses through `1/0`. Set reciprocal singular values to zero below an appropriate rank threshold, then check all four Penrose equations. |
| 58 | [SVD.svdcmp](../sources/Decomposition/SVD.cs#L116), `DecompositionAuditTests` | The `4x4` rank-one matrix `A[i,j]=(i+1)*(j+1)` fails reconstruction even with 100 iterations. This precedes pseudoinverse construction. Investigate QR/deflation and convergence handling separately from group 57; the exact internal fault is not yet isolated. |
| 59 | [GEVD initialization](../sources/Decomposition/GEVD.cs#L48), `DecompositionAuditTests` | Finite nonsingular matrix pairs return nonfinite eigenvectors. The accumulated eigenvector matrix `Z` starts as zero rather than identity, then normalization divides by zero. QZ equation checks on the same pairs pass. |
| 60 | [LaplacianPyramidTransform.Backward](../sources/Transform/LaplacianPyramidTransform.cs#L152), `TransformAuditTests` | Vector overloads set `nlev=pyramid.Length` and access `pyramid[nlev]`. Start from the final valid index, as the matrix overload does. Laplacian/local-Laplacian vector filters inherit the exception. |
| 61 | [HankelTransform](../sources/Transform/HankelTransform.cs#L56), `AdditionalTransformAuditTests` | Order 5, size 2: matrix entry `(0,0)` is about `-2.116e10`, expected `0.524632933`. The unbracketed zero search uses the inaccurate Bessel J implementation, binary32 arithmetic, and a `1e-16` stopping tolerance. Validate ordered positive roots and convergence after repairing Bessel J. Smaller discrepancies also occur at lower orders. |
| 62 | [BilateralGridFilter](../sources/Transform/BilateralGridFilter.cs), `TransformAuditTests` | A constant complex input `.375+.125i` becomes `.3952847+0i`: magnitude is returned and phase disappears. This is a **contract question** for the complex overload, with an enabled phase-preservation test. Either preserve complex samples or explicitly document magnitude-only filtering. |
| 63 | [LocalLaplacianFilter](../sources/Transform/LocalLaplacianFilter.cs), `TransformAuditTests` | Constant real input `.375` becomes `.3924532` with five intensity samples. This is a measured **approximation/DC-bias issue**. Review interpolation of the coarsest band; define and test an error budget if the bias is intentional. |
| 64 | [ThresholdFilter](../sources/Transform/ThresholdFilter.cs), `TransformAuditTests` | Complex matrix Over/Under thresholding clears an entire element where the vector overload clears individual components. This is a reproducible **overload-contract inconsistency**; choose and apply one rule. |
| 65 | [Meyer](../sources/Wavelet/Meyer.cs), `WaveletAuditTests` | Removable singularities are not handled: `Wavelet(.5)=4/pi`, `Scaling(+/-.75)=2/(3*pi)` must be finite. Evaluate these limits explicitly or use stable local expansions. |
| 86 | [Special.FactorialDown / Binomial](../sources/Core/Special.cs), reference suite | For nonnegative integer `n` and integer `k>n`, falling factorials/binomial coefficients should be zero; gamma ratios at poles return `NaN`. Handle integer identities before applying gamma ratios. |
| 87 | [ConeShape.Distribution](../sources/Distribution/ConeShape.cs), `AdditionalTransformAuditTests` | Negative `tau` yields a negative pulse height. The inverse transform of the even sinc kernel requires width `abs(tau)` and height proportional to `1/abs(tau)`. |
| 88 | [Special](../sources/Core/Special.cs), `SpecialFunctionReferenceTests` | **Remaining accuracy/branch work:** 358 reference cases fail in the enlarged special-function grid. This includes original groups 01–14 and 28–29, their dependent functions, and additional discrepancies. The full list is retained rather than treating every failed point as a new root cause. Repair foundational functions first, then re-evaluate Bessel I/Y/J/K, Struve H/L, Fresnel, Dawson, Faddeeva, generalized erf, and hypergeometric accuracy. |

For generalized complex erf, the reference uses the analytic integral continuation.
The current lower-gamma composition also introduces branch changes through `x^n`.
The complex API needs an explicit continuation convention before all such
differences can be classified as bugs. They remain visible in the reference suite.

Original group **05** is now exercised across all available reconstruction banks.
There are 52 failing Bior/CDF reconstruction cases; this is not evidence of 52
unrelated defects. Legendre and the explicitly non-perfect-reconstruction prototype
banks are not forced through an orthogonal perfect-reconstruction requirement.

### Images, geometry, rendering, and video

| ID | Source / test suite | Finding and repair direction |
| --- | --- | --- |
| 66 | [AHSL](../sources/Colorspace/AHSL.cs), `ColorSpaceAuditTests` | RGB `(64,0,128)` loses chroma and becomes gray `(64,64,64)`. Saturation collapses when red equals the channel average. Derive chroma from all channels. |
| 76 | [Intensity / Correction](../sources/Imaging/Intensity.cs), `ImagingAuditTests` | Neutral gamma/brightness/contrast/linear/levels/shift/transparency settings can reduce an 8-bit channel by one; white becomes 254. Align lookup-table sample coordinates with the scale used when applying the table. A separate gamma pixel-equation test detects the same normalization problem. |
| 77 | [RGBFilter](../sources/Imaging/RGBFilter.cs), [Grayscale](../sources/Imaging/Grayscale.cs), [TransparencyCorrection](../sources/Imaging/TransparencyCorrection.cs), [ErrorDiffusionDithering](../sources/Imaging/ErrorDiffusionDithering.cs) | Padded and negative-stride `BitmapData` cases produce wrong pixels or write into padding. Address each row through `Scan0+y*Stride`; do not assume contiguous forward rows. Guarded managed buffers keep the tests inside allocated memory. |
| 78 | [Statistics.Median](../sources/Imaging/Statistics.cs#L471), `ImagingAuditTests` | Odd populations select an earlier observation; one sample at intensity 20 yields median 0. Correct the cumulative-count target and its comparison. |
| 79 | [TensorMatrix](../sources/Imaging/TensorMatrix.cs#L256), `ImagingAuditTests` | RGB byte/float tensor reconstruction leaves alpha zero, making the resulting ARGB bitmap transparent. Set opaque alpha when the input has no alpha channel. |
| 80 | [Rectangles.IoU](../sources/Imaging/Rectangles.cs#L306), `GeometryAndRenderingAuditTests` | A `50000x50000` rectangle compared with itself gives about `.56004649`, expected `1`. Widen operands before area multiplication. |
| 81 | [Figure](../sources/Visualization/Figure.cs), `GeometryAndRenderingAuditTests` | Automatic rendering of constant `y=2` series throws `OverflowException`. Expand zero-width/height data ranges before calculating tick spacing and coordinate scales. All 15 tested constant-series style combinations fail. |
| 82 | [Boundary.FromResponse](../sources/Video/Internal/Boundary.cs), `VideoAuditTests` | MIME boundary parsing includes trailing parameters in the boundary value and rejects case variants of media types/attribute names. Parse parameter tokens independently and compare their names case-insensitively. |
| 83 | [MJPEGStreamParser](../sources/Video/Internal/MJPEGStreamParser.cs#L211), `VideoAuditTests` | Chunk sizes 1, 2, and 3 can produce a negative search start and throw before enough header bytes arrive. Bound the retained search position to the valid buffer prefix. |
| 84 | [ColorTransfer](../sources/Imaging/ColorTransfer.cs), `RemainingImagingAuditTests` | Transferring an identical constant RGB image `(77,99,123)` produces black. The variance ratio becomes `0/0`; additionally, `StnDev().StnDev()` is not the standard deviation of all pixels. Define the zero-variance case and compute channel moments over the image. |
| 85 | [DepthTransform.Equalize](../sources/Imaging/DepthTransform.cs#L500), `UtilityContractAuditTests` | A constant `256x256` depth image maps to zero instead of 65535. Histogram bins are `ushort` and wrap at 65536. Use counters wide enough for the image population. |

## Recommended repair order

1. Fix nontermination and basic indexing: 47, 56, 59, 60, 67–75, 77, 83, 85.
2. Repair shared arithmetic and special functions: original 01–04, 06–14, 17,
   26–30, plus 48–55 and 86. Re-run dependent distribution/transform suites before
   counting residual failures as independent causes.
3. Repair reconstruction and distribution formulas: original 05 and 18–25/31,
   plus 32–46, 57–58, 65–66, 76, 78–81, 84, 87.
4. Resolve accuracy and contract decisions in 40, 61–64, 82, 88. Set explicit
   accuracy domains rather than increasing tolerances until the suite passes.

Every fix should retain its regression test and add a nearby parameter/shape case
where necessary. Run the complete audit after shared mathematical primitives
change. Tests marked as contract/invariant checks supplement independent references;
they must not replace them.

## Limits and remaining work

- **6,006 executable lines and 3,433 branches remain unexecuted.** The exact list
  is exported, not hidden behind a single percentage. The 1,320 zero-hit method
  entries include accessors, constructors, compiler-generated functions and
  alternate implementations; they are not 1,320 distinct mathematical algorithms.
- Many image filters are checked for constant preservation, neutral settings,
  channel behavior, or stride handling rather than against a complete independent
  implementation on arbitrary images. In particular, artistic/nonlinear filters,
  all color-space parameter combinations, unusual pixel formats, large images,
  and every boundary policy need more reference cases.
- Internal alternative morphology implementations, general random-matrix
  generators, parsing/formatting variants, diagnostic printing, and some setters
  and invalid-argument paths are not exhaustively exercised. Public mathematical
  paths were prioritized over coverage of simple accessors.
- Probability reference data samples selected parameter sets. Unimplemented
  getters (`NotSupportedException` or explicitly undefined `NaN` results) are
  recorded in the distribution catalog and contract tests; they are not counted
  as verified numerical answers. Parameter-boundary and diverging-moment cases
  are still incomplete. Approximate medians require an explicit accuracy budget.
- Ill-conditioned, very large, adversarial, and all repeated-eigenvalue matrix
  cases are not exhausted. Some production decompositions use internal random
  initialization without a public seed; the tests allow explicit residual budgets.
- Reference fixtures exclude unrepresentably large finite values, poles, and
  many branch-cut limits. Absolute error allowances do not establish good relative
  accuracy near zero. One negative-real-part complex LogGamma fixture remains
  excluded because its unwrapped branch is unspecified.
- Live HTTP streams, cameras, and actual screen capture were not invoked. Video
  coverage comprises parser tests, generated in-memory frames, and unstarted
  configuration objects. No performance, race-condition, resource-leak, unsafe
  memory proof, or complete platform-compatibility audit is claimed.
- The complete bitmap/rendering run requires Windows and the .NET 8 runtime.
  `[SupportedOSPlatform]` documents this requirement; it does not skip tests on
  other operating systems. The numeric-only test filter is documented in the README.

Several initial suspicions were deliberately rejected after checking conventions:
complex two-dimensional transform matrices use conjugation on the right; Weyl–
Heisenberg inversion requires a nondegenerate translated-window system; real
arccotangent uses its documented positive-real branch; morphology arguments in
`Matrice` are halved before reaching the internal radius; and a two-level diffusion
palette changes at 129 in this API. The final tests reflect these conventions.
The perspective helper's corner-order comments should also be clarified: the
public warp path uses top-left, top-right, bottom-left, bottom-right. No numerical
perspective defect is reported for that documentation discrepancy.

## Independent definitions used

- [mpmath 1.3.0 documentation](https://mpmath.org/doc/1.3.0/) for high-precision
  special-function evaluation, quadrature, and positive Bessel zeros.
- [NIST DLMF incomplete gamma definitions](https://dlmf.nist.gov/8.2),
  [continued fractions](https://dlmf.nist.gov/8.9),
  [Bessel integral representations](https://dlmf.nist.gov/10.9), and
  [large-argument asymptotics](https://dlmf.nist.gov/10.17).
- [NIST power normal](https://www.itl.nist.gov/div898/handbook/eda/section3/eda366d.htm)
  and [power lognormal](https://www.itl.nist.gov/div898/handbook/eda/section3/eda366e.htm)
  definitions for group 40.
- [SciPy Sokal–Sneath definition](https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/generated.scipy.spatial.distance.sokalsneath.html)
  for group 44.
- [NIST Fourier transforms](https://dlmf.nist.gov/1.14) for the even sinc/rectangle
  pair in group 87.
- [RFC 2045](https://www.rfc-editor.org/rfc/rfc2045) for media-type and parameter
  token handling in group 82.

The supplied tests and exported evidence make these findings reproducible. They
also make the remaining gaps explicit; absence of further errors cannot yet be
confirmed.
