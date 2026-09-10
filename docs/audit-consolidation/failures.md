# Failed audit checks

Generated from the recorded TRX. These are failing test families, not independent root causes.
The explanatory repair register is in [the expanded report](../math-audit-expanded-2026-09-10.md).

| Test family | Failed cases | Example |
| --- | ---: | --- |
| [AdditionalTransformAuditTests.HankelMatrixMatchesHighPrecisionBesselZeros](../../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 4 | Expected 0.00045426352213434666; actual NaN; tolerance 0.000500227. |
| [AdditionalTransformAuditTests.TimeFrequencyKernelsMatchTheirFourierPair](../../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 2 | Expected 1.248456843919346; actual -1.2484568357467651; tolerance 2.69691E-05. |
| [ColorSpaceAuditTests.CieXyzAndLabMatchTheD65ReferenceWhite](../../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 6 | Expected -107.85373425232731; actual -113.36367034912109; tolerance 0.00515707. |
| [ColorSpaceAuditTests.ColorConversionsPreserveColorsWithinTheirQuantizationBudget](../../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 13 | Expected RGB(64, 0, 128); actual RGB(64, 64, 64); channel tolerance 2. |
| [DecompositionAuditTests.GeneralizedFactorizationsSatisfyBothMatrixEquations](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 | GEVD returned nonfinite eigenvectors for a finite nonsingular matrix pair. |
| [DecompositionAuditTests.RectangularAndRankDeficientMatricesRetainTheirInformation](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 9 | Expected 0; actual NaN; tolerance 0.005. |
| [DecompositionAuditTests.SchurDecompositionOfTheZeroMatrixTerminates](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 | Assert.Equal() Failure: Strings differ |
| [GeometryAndRenderingAuditTests.FiguresRenderFiniteConstantAndDiscontinuousSeries](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 15 | System.OverflowException : Overflow error. |
| [GeometryAndRenderingAuditTests.RectangleOverlapUsesGeometricAreaWithoutIntegerOverflow](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Expected 1; actual 0.5600464940071106; tolerance 2.2E-05. |
| [GeometryAndRenderingAuditTests.RectangularDepthTransformsPreserveCoordinateMeaning](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.HistogramMedianSelectsTheMiddleObservationForOddCounts](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 3 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.NeutralFilterSettingsPreservePixels](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 7 | Expected RGBA 47,11,89,57; actual 46,10,88,57; tolerance 0. |
| [ImagingAuditTests.PixelFiltersHonorBitmapDataStrideAndLeavePaddingUntouched](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 8 | Expected RGBA 0,0,255,60; actual 0,0,0,60; tolerance 0. |
| [ImagingAuditTests.PixelwiseFiltersMatchIndependentChannelEquations](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 | Expected RGBA 216,11,140,125; actual 214,11,138,125; tolerance 1. |
| [ImagingAuditTests.TensorConversionsMatchPixelChannelOrder](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 4 | Assert.InRange() Failure: Value not in range |
| [MathematicalRegressionTests.BiorthogonalWaveletReconstructsImpulse](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected 1; actual 1.0156248807907104; tolerance 0.0001. |
| [MathematicalRegressionTests.EvenConfinedWindowIsSymmetric](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.13257220387458801; actual -0.16181044280529022; tolerance 0.0001. |
| [MathematicalRegressionTests.EvenNormalWindowIsSymmetric](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.36787945032119751; actual 0.16901330649852753; tolerance 0.0001. |
| [MathematicalRegressionTests.IirStabilityUsesTheSamePolynomialAsReaction](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Assert.False() Failure |
| [MathematicalRegressionTests.RybConversionPreservesWhite](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [MathematicalRegressionTests.XyzCanRepresentItsD65WhitePoint](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1.089; actual 1; tolerance 2.378E-05. |
| [RemainingImagingAuditTests.ConstantColorTransferPreservesIdenticalImages](../../tests/UMapx.Tests/RemainingImagingAuditTests.cs) | 2 | Expected RGBA 77,99,123,255; actual 0,0,0,255; tolerance 1. |
| [ResponseAuditTests.StabilityMatchesThePoleOfTheActualDifferenceEquation](../../tests/UMapx.Tests/ResponseAuditTests.cs) | 2 | Assert.Equal() Failure: Values differ |
| [TransformAuditTests.LaplacianPyramidsReconstructVectors](../../tests/UMapx.Tests/TransformAuditTests.cs) | 6 | System.IndexOutOfRangeException : Index was outside the bounds of the array. |
| [TransformAuditTests.SmoothingAndDetailFiltersPreserveConstants](../../tests/UMapx.Tests/TransformAuditTests.cs) | 6 | Assert.All() Failure: 16 out of 16 items in the collection did not pass. |
| [TransformAuditTests.ThresholdFiltersRespectEqualityAndSignedComponents](../../tests/UMapx.Tests/TransformAuditTests.cs) | 2 | Expected <0; -1>; actual (0, 0); error 1. |
| [UtilityContractAuditTests.DepthHistogramEqualizationCountsMoreThan65535PixelsWithoutOverflow](../../tests/UMapx.Tests/UtilityContractAuditTests.cs) | 2 | Assert.Equal() Failure: Values differ |
| [VideoAuditTests.MjpegFramesSurviveArbitraryTransportChunkBoundaries](../../tests/UMapx.Tests/VideoAuditTests.cs) | 3 | System.ArgumentOutOfRangeException : Index was out of range. Must be non-negative and less than or equal to the size of the collection. (Parameter 'startIndex') |
| [VideoAuditTests.MultipartBoundaryParsingHonorsMimeParameters](../../tests/UMapx.Tests/VideoAuditTests.cs) | 2 | System.ArgumentException : Invalid content type |
| [WaveletAuditTests.MeyerWaveletAndScalingHaveFiniteRemovableSingularities](../../tests/UMapx.Tests/WaveletAuditTests.cs) | 3 | Expected 0.21220659078919379; actual 0; tolerance 1.42441E-05. |
| [WaveletAuditTests.ReconstructionBanksRecoverImpulsesAndDeterministicSignals](../../tests/UMapx.Tests/WaveletAuditTests.cs) | 52 | Expected 0.15000000596046448; actual 0.12346269190311432; tolerance 0.0005. |
| [WindowAuditTests.ExplicitFrameSizeDoesNotDependOnStoredFrameSize](../../tests/UMapx.Tests/WindowAuditTests.cs) | 5 | Expected 0.99660146236419678; actual 2.8826443667639978E-06; tolerance 2.9932E-05. |
| [WindowAuditTests.SamplesAgreeWithIndependentWindowFormulas](../../tests/UMapx.Tests/WindowAuditTests.cs) | 11 | Expected 0; actual 0.75999999046325684; tolerance 5E-05. |
