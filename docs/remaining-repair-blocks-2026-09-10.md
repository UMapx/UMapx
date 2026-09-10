# Remaining repair blocks — September 10, 2026

After the [matrix and distribution repairs](matrix-distribution-repair-2026-09-10.md),
**198 failures remain in eight open blocks**. B01–B04 are closed. The latest step
resolved all 113 B03/B04 failures and added 1,982 passing cases, with no regressions
or removed tests. Every remaining failing test ID belongs to exactly one block.
The complete verified run has 14,112 cases, with 13,914 passed and no skipped cases.

**Recommended next: B05 — approximation, interpolation, and local array filters (15 failures).**
Its interpolation and resampling fixes support subsequent transform and imaging work.
No failures outside B03/B04 changed in the latest run. The decomposition failures
remain reproducible after the matrix repairs.

The table is a suggested work order. Blocks without dependencies can be taken earlier.
The IsPrime(1) termination defect is fixed. Zero-matrix Schur remains a termination
defect and should be the first subtask in B06. Keep the subprocess timeouts enabled.

| Block | Scope | Assigned failing cases | Recommended after |
| --- | --- | ---: | --- |
| [B01](#b01) | Real and complex arithmetic | 0 (31 resolved) | Closed |
| [B02](#b02) | Exact integer arithmetic and number theory | 0 (108 resolved) | Closed |
| [B03](#b03) | Matrix indexing, complex statistics, and distances | 0 (46 resolved) | Closed |
| [B04](#b04) | Probability distributions | 0 (67 resolved) | Closed |
| [B05](#b05) | Approximation, interpolation, and local array filters | 15 | B03 |
| [B06](#b06) | Matrix decompositions | 15 | B01, B03 |
| [B07](#b07) | Wavelets | 57 | B01 |
| [B08](#b08) | Window functions | 18 | B01 |
| [B09](#b09) | Transforms and response filters | 23 | B01, B05 |
| [B10](#b10) | Color spaces | 21 | Independent |
| [B11](#b11) | Images, geometry, and rendering | 44 | B03, B05, B10 |
| [B12](#b12) | Video stream parsing | 5 | Independent |
| **Remaining** | | **198** | |

Counts represent failed test cases, not independent bugs. Several parameter sets and original
regressions can expose the same cause. Fixing one block may also change another block's count.
Contract/approximation questions are identified below rather than presented as proven formula errors.

## Completion rule for each block

1. Reproduce and localize its assigned failures; retain the original case IDs as evidence.
2. Correct the algorithms and add independent references or identities for the newly fixed edge cases.
3. Run the focused methods and their related passing tests. Do not hide failures or loosen tolerances just to pass.
4. Run the complete audit and compare IDs with the previous snapshot; require no passing tests to regress.
5. Refresh the remaining-failure ownership and counts before selecting the next block.

For a contract question, record the selected mathematical/API definition and justify any test change explicitly.
A block is complete only when its assigned failures have been resolved or separately and explicitly adjudicated.

<a id="b01"></a>
## B01. Real and complex arithmetic — closed (31 cases resolved)

The following list records the repaired baseline failures. See the [repair report](arithmetic-repair-2026-09-10.md) for additional boundary tests.

Sources: [Maths.cs](../sources/Core/Maths.cs), [Complex32.cs](../sources/Core/Complex32.cs).

- Complex division under extreme scaling: 2 failing cases.
- Real Tanh, Ctanh, and Asinh: 12 cases including the original regressions.
- Complex Actan/Acosh branches and a negative real base with a complex exponent: 12 cases.
- Quadratic equations with a zero linear coefficient and real cube roots in the cubic solver: 5 cases.

**Validation:** Check principal branches, scaling invariance, finite limiting values, polynomial residuals, and Vieta identities; retain all passing elementary-function and special-function checks.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [CoreScalarAuditTests.ComplexElementaryFunctionsRespectPrincipalValues](../tests/UMapx.Tests/CoreScalarAuditTests.cs) | 9 |
| [CoreScalarAuditTests.RealElementaryFunctionsAgreeWithDoublePrecision](../tests/UMapx.Tests/CoreScalarAuditTests.cs) | 10 |
| [MathematicalRegressionTests.AsinhAvoidsCancellationOnNegativeArguments](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.ComplexAcoshUsesPrincipalBranchOnNegativeRealAxis](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.ComplexArccotangentAgreesWithPositiveRealBranch](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.ComplexDivisionIsInvariantToScaling](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 |
| [MathematicalRegressionTests.CubicUsesRealCubeRootsForNegativeRadicands](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 |
| [MathematicalRegressionTests.NegativeRealBaseSupportsComplexExponent](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.QuadraticHandlesZeroLinearCoefficient](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 3 |
| [MathematicalRegressionTests.TanhAvoidsInfinityDividedByInfinity](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |

</details>

<a id="b02"></a>
## B02. Exact integer arithmetic and number theory — closed (108 cases resolved)

The following list records the repaired baseline failures. All related old and new integer checks pass.

Sources: [Maths.cs](../sources/Core/Maths.cs).

- First fix termination of IsPrime(1) for both integer widths (2 cases); keep the process timeouts enabled.
- Then repair primality/factorization and the dependent radical/totient calculations.
- Normalize GCD/Bezout signs; remove premature overflow and float arithmetic from LCM and modular powers.
- Repair coprime search, digit order, integer base conversion, and zero/exponent-zero handling.

**Validation:** Use BigInteger references, factor-product and primality checks, Bezout identities, and both integer widths. The 89 parameterized IntegerOperations cases cover several operations, not 89 separate defects.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [NumberTheoryAuditTests.BaseConversionsPreserveExactIntegerValues](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 3 |
| [NumberTheoryAuditTests.CompositeNumbersAreNotDeclaredPrimeWhenPollardFailsToSplitThem](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 |
| [NumberTheoryAuditTests.CoprimeSearchActuallyReturnsACoprime](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 |
| [NumberTheoryAuditTests.DecimalDigitVectorsRoundTripIndependentlyOfBaseConversion](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 3 |
| [NumberTheoryAuditTests.IntegerGcdLcmAndBezoutAgreeWithExactArithmetic](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 5 |
| [NumberTheoryAuditTests.IntegerOperationsIndependentlyMatchExactArithmetic](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 89 |
| [NumberTheoryAuditTests.ModularExponentiationMatchesBigInteger](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 1 |
| [NumberTheoryAuditTests.OneIsNotPrimeAndTheCheckTerminates](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 |
| [NumberTheoryAuditTests.PrimeSieveFactorizationTotientAndRadicalAgreeWithIntegerDefinitions](../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 1 |

</details>

<a id="b03"></a>
## B03. Matrix indexing, complex statistics, and distances — closed (46 cases resolved)

The following list records the repaired baseline failures. The
[repair report](matrix-distribution-repair-2026-09-10.md) also records 316 added
matrix/distance cases, including extreme shifts and complex covariance scaling.

Sources: [Matrice.cs](../sources/Core/Matrice.cs), [LinealgOptions.cs](../sources/Core/LinealgOptions.cs), [SokalSneath.cs](../sources/Distance/SokalSneath.cs).

- Rectangular indexing and placement: diagonal multiplication 8, Shift 4, Merge 6, mesh evaluation 8, Swap 12; total 38 cases.
- Hermitian variance, covariance, and norms: 5 cases including the two original regressions.
- Sokal-Sneath denominator: 3 independent contingency-table checks.

**Validation:** Use non-square arrays in both orientations, all real/complex overloads, independent index maps, and conjugate inner products. Verify covariance symmetry and nonnegative norms.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [DistanceAuditTests.BooleanDistancesUseContingencyCounts](../tests/UMapx.Tests/DistanceAuditTests.cs) | 3 |
| [MathematicalRegressionTests.ComplexVarianceUsesSquaredMagnitudes](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.ComplexVectorModulusCannotCancelNonzeroComponents](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MatrixAuditTests.ComplexStatisticsRespectHermitianInnerProducts](../tests/UMapx.Tests/MatrixAuditTests.cs) | 3 |
| [MatrixAuditTests.DiagonalMultiplicationScalesTheDocumentedRowsOrColumns](../tests/UMapx.Tests/MatrixAuditTests.cs) | 8 |
| [MatrixAuditTests.RectangularArrayOperationsMatchIndexDefinitions](../tests/UMapx.Tests/MatrixAuditTests.cs) | 10 |
| [MatrixStructureAuditTests.MeshEvaluationUsesBothIndependentCoordinates](../tests/UMapx.Tests/MatrixStructureAuditTests.cs) | 8 |
| [MatrixStructureAuditTests.SwappingRowsAndColumnsPermutesEveryEntry](../tests/UMapx.Tests/MatrixStructureAuditTests.cs) | 12 |

</details>

<a id="b04"></a>
## B04. Probability distributions — closed (67 cases resolved)

The following list records the repaired baseline failures. The
[repair report](matrix-distribution-repair-2026-09-10.md) adds 1,666 passing
distribution cases, including 1,611 independent high-precision references.

Sources: [Distribution](../sources/Distribution).

- PowerNormal and PowerLognormal: 29 reference failures; align the documented law, PDF, and CDF.
- Binomial and Poisson: 9 failures across boundaries, mass, entropy, modes, and medians.
- Birnbaum-Saunders: 10 failures across median, mode, and entropy.
- Other moments/entropy: 16; FoldedNormal mode and ChiSquare median checks: 3.

**Validation:** Use independent PDF integrals, CDF/PDF differentiation, normalization, parameter-dependent existence of moments, and both median inequalities. Recheck each distribution after shared scalar fixes.

**Resolved contracts:** PowerNormal/PowerLognormal follow NIST's survival-power laws;
ChiSquare median now uses CDF inversion. Discrete medians select the lower median.
The 67 resolved cases include 50 reference checks, 11 shape checks, and 6 original
regressions. Existing assertions and tolerances were retained.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [DistributionReferenceTests.MatchesIndependentProbabilityReference](../tests/UMapx.Tests/DistributionReferenceTests.cs) | 50 |
| [DistributionShapeAuditTests.BirnbaumSaundersModeIsThePositiveStationaryPoint](../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 2 |
| [DistributionShapeAuditTests.FoldedNormalWithLocationSmallerThanScaleHasModeZero](../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 1 |
| [DistributionShapeAuditTests.ModesAndMediansSatisfyTheirProbabilityDefinitions](../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 8 |
| [MathematicalRegressionTests.BinomialCertainSuccessHasUnitMass](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 |
| [MathematicalRegressionTests.BinomialMedianSatisfiesBothHalfProbabilityInequalities](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.InverseChiSquareEntropyAgreesWithInverseGamma](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.PoissonMassIsFiniteAtItsModeForLambda100](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.PoissonMedianCanBeOneBelowLambdaOne](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |

</details>

<a id="b05"></a>
## B05. Approximation, interpolation, and local array filters — 15 cases

Sources: [Pade.cs](../sources/Analysis/Pade.cs), [Interpolation.cs](../sources/Analysis/Interpolation.cs), [LinealgOptions.cs](../sources/Core/LinealgOptions.cs), [Resize.cs](../sources/Imaging/Resize.cs).

- Pade denominator degree greater than numerator degree: 2 cases.
- Bilinear grid-edge interpolation: 2 cases.
- Bicubic resizing at scale one: 3 cases (real array, complex array, bitmap).
- Local mean: 2; morphology ranks: 4; exact half-turn rotation: 2.

**Validation:** Compare Taylor coefficients, exact affine interpolation, impulse responses, sorted edge neighborhoods, and pixel/sample identities. Check endpoints and rectangular shapes.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [AnalysisAuditTests.PadeCoefficientsMatchTheTaylorSeriesThroughTheRequestedOrder](../tests/UMapx.Tests/AnalysisAuditTests.cs) | 1 |
| [ImagingAuditTests.BitmapResizingToTheSameDimensionsIsAnIdentity](../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 |
| [KernelAndContainerAuditTests.ResizingAnArrayToItsCurrentSizePreservesSamples](../tests/UMapx.Tests/KernelAndContainerAuditTests.cs) | 2 |
| [MathematicalRegressionTests.BilinearInterpolationIsLinearAlongGridEdges](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 |
| [MathematicalRegressionTests.PadeSupportsDenominatorDegreeLargerThanNumeratorDegree](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MatrixAuditTests.LocalMeanHasTheCorrectInteriorImpulseResponse](../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 |
| [MatrixAuditTests.MorphologyMatchesSortingOfReplicatedEdgeNeighborhoods](../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 |
| [MatrixAuditTests.MorphologyUsesZeroBasedRanksInAThreeSampleWindow](../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 |
| [MatrixFilterAuditTests.MatrixRotationsAtExactHalfTurnsMatchIndexReversal](../tests/UMapx.Tests/MatrixFilterAuditTests.cs) | 2 |

</details>

<a id="b06"></a>
## B06. Matrix decompositions — 15 cases

Sources: [Schur.cs](../sources/Decomposition/Schur.cs), [SVD.cs](../sources/Decomposition/SVD.cs), [GEVD.cs](../sources/Decomposition/GEVD.cs).

- Handle Schur zero-matrix deflation and termination first: 3 cases.
- Separate SVD reconstruction from pseudoinverse singular-value thresholding: 9 cases.
- Repair GEVD eigenvector accumulation/normalization: 3 cases.

**Validation:** Check reconstruction, orthogonality/unitarity, all four Penrose equations, generalized eigen-equations, and bounded termination on zero/rank-deficient inputs.

**Open point:** The exact internal cause of the rank-one SVD reconstruction failure is not fully localized. Fixing reciprocal zero singular values alone is insufficient to declare the SVD block complete.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [DecompositionAuditTests.GeneralizedFactorizationsSatisfyBothMatrixEquations](../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 |
| [DecompositionAuditTests.RectangularAndRankDeficientMatricesRetainTheirInformation](../tests/UMapx.Tests/DecompositionAuditTests.cs) | 9 |
| [DecompositionAuditTests.SchurDecompositionOfTheZeroMatrixTerminates](../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 |

</details>

<a id="b07"></a>
## B07. Wavelets — 57 cases

Sources: [Wavelet](../sources/Wavelet).

- Bior/CDF reconstruction banks: 54 cases including 2 original Bior13 regressions.
- Meyer removable singularities: 3 cases.

**Validation:** Check synthesis/analysis coefficient pairing, normalization, phase/alignment, impulses, deterministic signals, and analytical limits. Retain the distinction between reconstruction banks and prototype wavelets that do not promise perfect reconstruction.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [MathematicalRegressionTests.BiorthogonalWaveletReconstructsImpulse](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 |
| [WaveletAuditTests.MeyerWaveletAndScalingHaveFiniteRemovableSingularities](../tests/UMapx.Tests/WaveletAuditTests.cs) | 3 |
| [WaveletAuditTests.ReconstructionBanksRecoverImpulsesAndDeterministicSignals](../tests/UMapx.Tests/WaveletAuditTests.cs) | 52 |

</details>

<a id="b08"></a>
## B08. Window functions — 18 cases

Sources: [BartlettHann.cs](../sources/Window/BartlettHann.cs), [Normal.cs](../sources/Window/Normal.cs), [Confined.cs](../sources/Window/Confined.cs).

- Independent window formulas: 11 cases, including the Bartlett-Hann sign and Gaussian centering.
- Explicit frame size versus stored frame size: 5 cases.
- Original even Normal/Confined symmetry regressions: 2 cases.

**Validation:** Check odd/even lengths, endpoint values, symmetry, independent sample formulas, and independence from stored size when an explicit size is supplied.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [MathematicalRegressionTests.EvenConfinedWindowIsSymmetric](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.EvenNormalWindowIsSymmetric](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [WindowAuditTests.ExplicitFrameSizeDoesNotDependOnStoredFrameSize](../tests/UMapx.Tests/WindowAuditTests.cs) | 5 |
| [WindowAuditTests.SamplesAgreeWithIndependentWindowFormulas](../tests/UMapx.Tests/WindowAuditTests.cs) | 11 |

</details>

<a id="b09"></a>
## B09. Transforms and response filters — 23 cases

Sources: [Transform](../sources/Transform), [IIR.cs](../sources/Response/IIR.cs), [ConeShape.cs](../sources/Distribution/ConeShape.cs).

- Hankel positive-zero bracketing and normalization: 4 cases.
- Laplacian pyramid/vector filter reconstruction and constant preservation: 8 cases.
- BilateralGrid and LocalLaplacian constant/phase behavior: 4 cases.
- Complex threshold overload consistency: 2 cases.
- IIR characteristic polynomial/stability: 3 cases; ConeShape sign for negative tau: 2 cases.

**Validation:** Check ordered positive Bessel zeros, independent transform matrices, reconstruction, DC/phase preservation, threshold equality, and poles of the actual difference equation.

**Open point:** BilateralGrid complex semantics, LocalLaplacian approximation/DC bias, and component-versus-magnitude thresholding need explicit API decisions. These are six failing expectations, not six already-localized algebraic mistakes. The order-5 Hankel search currently converges toward the origin.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [AdditionalTransformAuditTests.HankelMatrixMatchesHighPrecisionBesselZeros](../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 4 |
| [AdditionalTransformAuditTests.TimeFrequencyKernelsMatchTheirFourierPair](../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 2 |
| [MathematicalRegressionTests.IirStabilityUsesTheSamePolynomialAsReaction](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [ResponseAuditTests.StabilityMatchesThePoleOfTheActualDifferenceEquation](../tests/UMapx.Tests/ResponseAuditTests.cs) | 2 |
| [TransformAuditTests.LaplacianPyramidsReconstructVectors](../tests/UMapx.Tests/TransformAuditTests.cs) | 6 |
| [TransformAuditTests.SmoothingAndDetailFiltersPreserveConstants](../tests/UMapx.Tests/TransformAuditTests.cs) | 6 |
| [TransformAuditTests.ThresholdFiltersRespectEqualityAndSignedComponents](../tests/UMapx.Tests/TransformAuditTests.cs) | 2 |

</details>

<a id="b10"></a>
## B10. Color spaces — 21 cases

Sources: [XYZ.cs](../sources/Colorspace/XYZ.cs), [LAB.cs](../sources/Colorspace/LAB.cs), [RYB.cs](../sources/Colorspace/RYB.cs), [AHSL.cs](../sources/Colorspace/AHSL.cs).

- XYZ/LAB white point, range, and round trips: 15 cases including the original XYZ regression.
- RYB chroma/white preservation: 4 cases; AHSL chroma: 2 cases.

**Validation:** Use D65 reference coordinates, primary/neutral colors, gamut boundaries, and round trips with the existing quantization budgets. Repair XYZ before judging downstream LAB differences.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [ColorSpaceAuditTests.CieXyzAndLabMatchTheD65ReferenceWhite](../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 6 |
| [ColorSpaceAuditTests.ColorConversionsPreserveColorsWithinTheirQuantizationBudget](../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 13 |
| [MathematicalRegressionTests.RybConversionPreservesWhite](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |
| [MathematicalRegressionTests.XyzCanRepresentItsD65WhitePoint](../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 |

</details>

<a id="b11"></a>
## B11. Images, geometry, and rendering — 44 cases

Sources: [Imaging](../sources/Imaging), [Figure.cs](../sources/Visualization/Figure.cs).

- Lookup-table normalization and pixel equations: 8 cases.
- Padded/negative bitmap stride: 8; histogram median: 3; tensor alpha/channel reconstruction: 4.
- Constant-image ColorTransfer: 2; depth histogram counter overflow: 2.
- Figure rendering of constant series: 15; rectangle area overflow: 1; depth Merge coordinates: 1.

**Validation:** Use guarded buffers for stride/padding, exact neutral pixels, independent histogram ranks, opaque alpha for RGB inputs, zero-variance limits, wide area/count arithmetic, and finite rendering ranges.

**Open point:** The bicubic bitmap resize failure belongs to B05 and is not counted here. ColorTransfer should be rechecked after matrix statistics and color-space repairs.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [GeometryAndRenderingAuditTests.FiguresRenderFiniteConstantAndDiscontinuousSeries](../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 15 |
| [GeometryAndRenderingAuditTests.RectangleOverlapUsesGeometricAreaWithoutIntegerOverflow](../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 |
| [GeometryAndRenderingAuditTests.RectangularDepthTransformsPreserveCoordinateMeaning](../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 |
| [ImagingAuditTests.HistogramMedianSelectsTheMiddleObservationForOddCounts](../tests/UMapx.Tests/ImagingAuditTests.cs) | 3 |
| [ImagingAuditTests.NeutralFilterSettingsPreservePixels](../tests/UMapx.Tests/ImagingAuditTests.cs) | 7 |
| [ImagingAuditTests.PixelFiltersHonorBitmapDataStrideAndLeavePaddingUntouched](../tests/UMapx.Tests/ImagingAuditTests.cs) | 8 |
| [ImagingAuditTests.PixelwiseFiltersMatchIndependentChannelEquations](../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 |
| [ImagingAuditTests.TensorConversionsMatchPixelChannelOrder](../tests/UMapx.Tests/ImagingAuditTests.cs) | 4 |
| [RemainingImagingAuditTests.ConstantColorTransferPreservesIdenticalImages](../tests/UMapx.Tests/RemainingImagingAuditTests.cs) | 2 |
| [UtilityContractAuditTests.DepthHistogramEqualizationCountsMoreThan65535PixelsWithoutOverflow](../tests/UMapx.Tests/UtilityContractAuditTests.cs) | 2 |

</details>

<a id="b12"></a>
## B12. Video stream parsing — 5 cases

Sources: [Boundary.cs](../sources/Video/Internal/Boundary.cs), [MJPEGStreamParser.cs](../sources/Video/Internal/MJPEGStreamParser.cs).

- MIME boundary parameters/case handling: 2 cases.
- MJPEG headers split across one-, two-, and three-byte transport chunks: 3 cases.

**Validation:** Use in-memory streams with arbitrary split points, bounded buffer indices, quoted/trailing parameters, and case variants. No live camera or network is needed.

**Open point:** This is a parsing/contract block, not a numerical mathematics defect.

<details>
<summary>Assigned failing test families</summary>

| Test family | Cases |
| --- | ---: |
| [VideoAuditTests.MjpegFramesSurviveArbitraryTransportChunkBoundaries](../tests/UMapx.Tests/VideoAuditTests.cs) | 3 |
| [VideoAuditTests.MultipartBoundaryParsingHonorsMimeParameters](../tests/UMapx.Tests/VideoAuditTests.cs) | 2 |

</details>

## Evidence and commands

- [Current block metadata, focused filters, and all 198 unique test-ID assignments](audit-matrix-distribution/repair-blocks.json).
- [Current individual failures](audit-matrix-distribution/failures.json) and [run summary](audit-matrix-distribution/summary.json).
- [Matrix and distribution repair results](matrix-distribution-repair-2026-09-10.md).
- [Previous 311-case planning snapshot](audit-arithmetic/repair-blocks.json).
- [Arithmetic and number-theory repair results](arithmetic-repair-2026-09-10.md).
- [Previous 450-case planning snapshot](audit-special-functions/repair-blocks.json).
- [Special-function repair results](special-functions-repair-2026-09-10.md).
- [Historical repair register](math-audit-expanded-2026-09-10.md); its superseded special-function findings are not pending blocks.

Run a block's focused test families from the repository root:

```powershell
$plan = Get-Content -Raw docs/audit-matrix-distribution/repair-blocks.json | ConvertFrom-Json
$block = $plan.blocks | Where-Object id -eq 'B05'
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter $block.focused_test_filter
```

The focused filter includes all parameters of the assigned test methods, not just failing rows.
It is a starting point, not a substitute for related passing tests and the complete audit.

For B01 also run the complete CoreScalarAuditTests, MathematicalIdentityTests, and SpecialFunction suites
because changes to common arithmetic can affect previously repaired functions:

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~CoreScalarAuditTests|FullyQualifiedName~MathematicalIdentityTests|FullyQualifiedName~SpecialFunction"
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/next-block
```

The full suite currently exits with 1 because all 198 remaining expectations stay enabled.
The historical snapshots are retained; the current JSON matches the verified matrix/distribution repair run.
