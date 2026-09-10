# Failed audit checks

Generated from the recorded TRX. These are failing test families, not independent root causes.
The explanatory repair register is in [the expanded report](../math-audit-expanded-2026-09-10.md).

| Test family | Failed cases | Example |
| --- | ---: | --- |
| [AdditionalTransformAuditTests.HankelMatrixMatchesHighPrecisionBesselZeros](../../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 11 | Expected -0.32625831763452789; actual -0.32743236422538757; tolerance 0.000663129. |
| [AdditionalTransformAuditTests.TimeFrequencyKernelsMatchTheirFourierPair](../../tests/UMapx.Tests/AdditionalTransformAuditTests.cs) | 2 | Expected 1.248456843919346; actual -1.2484568357467651; tolerance 2.69691E-05. |
| [AnalysisAuditTests.PadeCoefficientsMatchTheTaylorSeriesThroughTheRequestedOrder](../../tests/UMapx.Tests/AnalysisAuditTests.cs) | 1 | System.IndexOutOfRangeException : Index was outside the bounds of the array. |
| [ColorSpaceAuditTests.CieXyzAndLabMatchTheD65ReferenceWhite](../../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 6 | Expected -107.85373425232731; actual -113.36367034912109; tolerance 0.00515707. |
| [ColorSpaceAuditTests.ColorConversionsPreserveColorsWithinTheirQuantizationBudget](../../tests/UMapx.Tests/ColorSpaceAuditTests.cs) | 13 | Expected RGB(64, 0, 128); actual RGB(64, 64, 64); channel tolerance 2. |
| [CoreScalarAuditTests.ComplexElementaryFunctionsRespectPrincipalValues](../../tests/UMapx.Tests/CoreScalarAuditTests.cs) | 9 | Expected <3.6907476555303655; -1.6206926493437195>; actual (-3.69070601, 1.62068748); error 8.06179. |
| [CoreScalarAuditTests.RealElementaryFunctionsAgreeWithDoublePrecision](../../tests/UMapx.Tests/CoreScalarAuditTests.cs) | 10 | Expected -5.2983423656105888; actual -5.2988667488098145; tolerance 0.00016195. |
| [DecompositionAuditTests.GeneralizedFactorizationsSatisfyBothMatrixEquations](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 | GEVD returned nonfinite eigenvectors for a finite nonsingular matrix pair. |
| [DecompositionAuditTests.RectangularAndRankDeficientMatricesRetainTheirInformation](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 9 | Expected 0; actual NaN; tolerance 0.005. |
| [DecompositionAuditTests.SchurDecompositionOfTheZeroMatrixTerminates](../../tests/UMapx.Tests/DecompositionAuditTests.cs) | 3 | Assert.Equal() Failure: Strings differ |
| [DistanceAuditTests.BooleanDistancesUseContingencyCounts](../../tests/UMapx.Tests/DistanceAuditTests.cs) | 3 | Expected 0.70588235294117652; actual 0.75; tolerance 1.61176E-05. |
| [DistributionReferenceTests.MatchesIndependentProbabilityReference](../../tests/UMapx.Tests/DistributionReferenceTests.cs) | 91 | Expected 0; actual NaN; tolerance 0.0003. |
| [DistributionShapeAuditTests.BirnbaumSaundersModeIsThePositiveStationaryPoint](../../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 2 | Expected 0.35320996419932449; actual -1; tolerance 3.70642E-05. |
| [DistributionShapeAuditTests.FoldedNormalWithLocationSmallerThanScaleHasModeZero](../../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 1 | Expected 0; actual 0.68305385112762451; tolerance 2E-06. |
| [DistributionShapeAuditTests.ModesAndMediansSatisfyTheirProbabilityDefinitions](../../tests/UMapx.Tests/DistributionShapeAuditTests.cs) | 12 | Assert.False() Failure |
| [GeometryAndRenderingAuditTests.FiguresRenderFiniteConstantAndDiscontinuousSeries](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 15 | System.OverflowException : Overflow error. |
| [GeometryAndRenderingAuditTests.RectangleOverlapUsesGeometricAreaWithoutIntegerOverflow](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Expected 1; actual 0.5600464940071106; tolerance 2.2E-05. |
| [GeometryAndRenderingAuditTests.RectangularDepthTransformsPreserveCoordinateMeaning](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.BitmapResizingToTheSameDimensionsIsAnIdentity](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 | Assert.InRange() Failure: Value not in range |
| [ImagingAuditTests.HistogramMedianSelectsTheMiddleObservationForOddCounts](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 3 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.NeutralFilterSettingsPreservePixels](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 7 | Expected RGBA 47,11,89,57; actual 46,10,88,57; tolerance 0. |
| [ImagingAuditTests.PixelFiltersHonorBitmapDataStrideAndLeavePaddingUntouched](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 8 | Expected RGBA 0,0,255,60; actual 0,0,0,60; tolerance 0. |
| [ImagingAuditTests.PixelwiseFiltersMatchIndependentChannelEquations](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 | Expected RGBA 216,11,140,125; actual 214,11,138,125; tolerance 1. |
| [ImagingAuditTests.TensorConversionsMatchPixelChannelOrder](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 4 | Assert.InRange() Failure: Value not in range |
| [KernelAndContainerAuditTests.ResizingAnArrayToItsCurrentSizePreservesSamples](../../tests/UMapx.Tests/KernelAndContainerAuditTests.cs) | 2 | Expected 0.93000000715255737; actual 1.0312892198562622; tolerance 3.86E-05. |
| [MathematicalRegressionTests.AsinhAvoidsCancellationOnNegativeArguments](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected -9.9034875550361274; actual -Infinity; tolerance 0.00020007. |
| [MathematicalRegressionTests.BesselJAsymptoticAccuracyDependsOnOrder](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.14853180559607407; actual 0.46753552556037903; tolerance 4.97064E-06. |
| [MathematicalRegressionTests.BesselJAtImaginaryArgumentSatisfiesConnectionToI](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <-145831809975.96713; 0>; actual (-7.20047899E+11, 1.10999201E+12); error 1.24972E+12. |
| [MathematicalRegressionTests.BesselYIncludesBothExponentialTerms](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.088256964215676956; actual 0.32845768332481384; tolerance 3.76514E-06. |
| [MathematicalRegressionTests.BetaAvoidsIntermediateGammaOverflow](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1.7681885473062026E-15; actual NaN; tolerance 4.53638E-20. |
| [MathematicalRegressionTests.BilinearInterpolationIsLinearAlongGridEdges](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected 0.5; actual 0; tolerance 1.2E-05. |
| [MathematicalRegressionTests.BinomialCertainSuccessHasUnitMass](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected 1; actual NaN; tolerance 2.2E-05. |
| [MathematicalRegressionTests.BinomialMedianSatisfiesBothHalfProbabilityInequalities](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1; actual 2; tolerance 2.2E-05. |
| [MathematicalRegressionTests.BiorthogonalWaveletReconstructsImpulse](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected 1; actual 1.0156248807907104; tolerance 0.0001. |
| [MathematicalRegressionTests.ChebyshevPolynomialIsDefinedOutsideUnitInterval](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 7; actual NaN; tolerance 0.000142. |
| [MathematicalRegressionTests.ChebyshevUHasFiniteEndpointValue](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 3; actual NaN; tolerance 6.2E-05. |
| [MathematicalRegressionTests.ComplexAcoshUsesPrincipalBranchOnNegativeRealAxis](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <1.3169578969248166; 3.141592653589793>; actual (-1.31695783, 3.14159274); error 2.63392. |
| [MathematicalRegressionTests.ComplexArccotangentAgreesWithPositiveRealBranch](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <0.7853981633974483; 0>; actual (-0.785398185, 0); error 1.5708. |
| [MathematicalRegressionTests.ComplexBesselKResolvesOscillatoryIntegral](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <3.015242162785077; -31.062215312887222>; actual (140853712, -48188600); error 1.48869E+08. |
| [MathematicalRegressionTests.ComplexDivisionIsInvariantToScaling](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected <1; 0>; actual (NaN, NaN); error NaN. |
| [MathematicalRegressionTests.ComplexErfAgreesWithEntireFunctionReference](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <-1.0035022433130363; 0.004740903031294336>; actual (-0.889444888, 0.677245617); error 0.682108. |
| [MathematicalRegressionTests.ComplexVarianceUsesSquaredMagnitudes](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <2; 0>; actual (-2, 0); error 4. |
| [MathematicalRegressionTests.ComplexVectorModulusCannotCancelNonzeroComponents](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <1.4142135623730951; 0>; actual (0, 0); error 1.41421. |
| [MathematicalRegressionTests.CubicUsesRealCubeRootsForNegativeRadicands](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 2 | Expected <0; 0>; actual (NaN, NaN); error NaN. |
| [MathematicalRegressionTests.EvenConfinedWindowIsSymmetric](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.13257220387458801; actual -0.16181044280529022; tolerance 0.0001. |
| [MathematicalRegressionTests.EvenNormalWindowIsSymmetric](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.36787945032119751; actual 0.16901330649852753; tolerance 0.0001. |
| [MathematicalRegressionTests.GammaContinuedFractionMustUseNextDenominator](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.19914827347145578; actual 0.17923341691493988; tolerance 5.98297E-06. |
| [MathematicalRegressionTests.GammaSeriesMustIncludeGammaNormalization](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.55950671493478765; actual 13.428162574768066; tolerance 1.31901E-05. |
| [MathematicalRegressionTests.HypergeometricSeriesAvoidsIntermediatePochhammerOverflow](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 21.789416887313024; actual NaN; tolerance 0.000437788. |
| [MathematicalRegressionTests.IirStabilityUsesTheSamePolynomialAsReaction](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Assert.False() Failure |
| [MathematicalRegressionTests.InverseChiSquareEntropyAgreesWithInverseGamma](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1.4612841492431206; actual 0.30685296654701233; tolerance 3.12257E-05. |
| [MathematicalRegressionTests.LambertWReturnsRequestedPrincipalBranch](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <0.21453028412564626; 0.3510955496497862>; actual (-1.87301111, -2.73923755); error 3.72934. |
| [MathematicalRegressionTests.NegativeRealBaseSupportsComplexExponent](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected <0; 1>; actual (NaN, NaN); error NaN. |
| [MathematicalRegressionTests.PadeSupportsDenominatorDegreeLargerThanNumeratorDegree](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | System.IndexOutOfRangeException : Index was outside the bounds of the array. |
| [MathematicalRegressionTests.PoissonMassIsFiniteAtItsModeForLambda100](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0.039860996809147099; actual NaN; tolerance 2.79722E-06. |
| [MathematicalRegressionTests.PoissonMedianCanBeOneBelowLambdaOne](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1; actual 0; tolerance 2.2E-05. |
| [MathematicalRegressionTests.PositiveOrderPochhammerOfZeroIsZero](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 0; actual 1; tolerance 2E-06. |
| [MathematicalRegressionTests.QuadraticHandlesZeroLinearCoefficient](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 3 | Expected <0; 0>; actual (-1, -0); error 1. |
| [MathematicalRegressionTests.RybConversionPreservesWhite](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [MathematicalRegressionTests.TanhAvoidsInfinityDividedByInfinity](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1; actual NaN; tolerance 2.2E-05. |
| [MathematicalRegressionTests.XyzCanRepresentItsD65WhitePoint](../../tests/UMapx.Tests/MathematicalRegressionTests.cs) | 1 | Expected 1.089; actual 1; tolerance 2.378E-05. |
| [MatrixAuditTests.ComplexStatisticsRespectHermitianInnerProducts](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 3 | Expected <14.33333333333333; 0>; actual (-5, -5.33333349); error 20.0555. |
| [MatrixAuditTests.DiagonalMultiplicationScalesTheDocumentedRowsOrColumns](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 8 | System.IndexOutOfRangeException : Index was outside the bounds of the array. |
| [MatrixAuditTests.LocalMeanHasTheCorrectInteriorImpulseResponse](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 | Expected 1; actual 0; tolerance 4E-05. |
| [MatrixAuditTests.MorphologyMatchesSortingOfReplicatedEdgeNeighborhoods](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 | Expected 0.94680589437484741; actual 0.86018317937850952; tolerance 2.09361E-05. |
| [MatrixAuditTests.MorphologyUsesZeroBasedRanksInAThreeSampleWindow](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 2 | Expected 4; actual 3; tolerance 8.2E-05. |
| [MatrixAuditTests.RectangularArrayOperationsMatchIndexDefinitions](../../tests/UMapx.Tests/MatrixAuditTests.cs) | 10 | Expected 2; actual 2.2578125; tolerance 6E-05. |
| [MatrixFilterAuditTests.MatrixRotationsAtExactHalfTurnsMatchIndexReversal](../../tests/UMapx.Tests/MatrixFilterAuditTests.cs) | 2 | Expected 2.2699999809265137; actual 2.0999999046325684; tolerance 0.000654. |
| [MatrixStructureAuditTests.MeshEvaluationUsesBothIndependentCoordinates](../../tests/UMapx.Tests/MatrixStructureAuditTests.cs) | 8 | Expected 3.4699999690055847; actual 3.25; tolerance 8.94E-05. |
| [MatrixStructureAuditTests.SwappingRowsAndColumnsPermutesEveryEntry](../../tests/UMapx.Tests/MatrixStructureAuditTests.cs) | 12 | System.IndexOutOfRangeException : Index was outside the bounds of the array. |
| [NumberTheoryAuditTests.BaseConversionsPreserveExactIntegerValues](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 3 | System.OverflowException : Arithmetic operation resulted in an overflow. |
| [NumberTheoryAuditTests.CompositeNumbersAreNotDeclaredPrimeWhenPollardFailsToSplitThem](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 | Assert.False() Failure |
| [NumberTheoryAuditTests.CoprimeSearchActuallyReturnsACoprime](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 | Assert.Equal() Failure: Values differ |
| [NumberTheoryAuditTests.DecimalDigitVectorsRoundTripIndependentlyOfBaseConversion](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 3 | System.OverflowException : Arithmetic operation resulted in an overflow. |
| [NumberTheoryAuditTests.IntegerGcdLcmAndBezoutAgreeWithExactArithmetic](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 5 | Assert.Equal() Failure: Values differ |
| [NumberTheoryAuditTests.IntegerOperationsIndependentlyMatchExactArithmetic](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 89 | Assert.Equal() Failure: Values differ |
| [NumberTheoryAuditTests.ModularExponentiationMatchesBigInteger](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 1 | System.OverflowException : Arithmetic operation resulted in an overflow. |
| [NumberTheoryAuditTests.OneIsNotPrimeAndTheCheckTerminates](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 2 | IsPrimeInt(1) did not terminate within five seconds. |
| [NumberTheoryAuditTests.PrimeSieveFactorizationTotientAndRadicalAgreeWithIntegerDefinitions](../../tests/UMapx.Tests/NumberTheoryAuditTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [RemainingImagingAuditTests.ConstantColorTransferPreservesIdenticalImages](../../tests/UMapx.Tests/RemainingImagingAuditTests.cs) | 2 | Expected RGBA 77,99,123,255; actual 0,0,0,255; tolerance 1. |
| [ResponseAuditTests.StabilityMatchesThePoleOfTheActualDifferenceEquation](../../tests/UMapx.Tests/ResponseAuditTests.cs) | 2 | Assert.Equal() Failure: Values differ |
| [SpecialFunctionReferenceTests.AgreesWithHighPrecisionReference](../../tests/UMapx.Tests/SpecialFunctionReferenceTests.cs) | 358 | Expected <1.7681885473062026E-15; 0>; actual (NaN, 0); error NaN. |
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
