# Failed audit checks

Generated from the recorded TRX. These are failing test families, not independent root causes.
The explanatory repair register is in [the expanded report](../math-audit-expanded-2026-09-10.md).

| Test family | Failed cases | Example |
| --- | ---: | --- |
| [GeometryAndRenderingAuditTests.FiguresRenderFiniteConstantAndDiscontinuousSeries](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 15 | System.OverflowException : Overflow error. |
| [GeometryAndRenderingAuditTests.RectangleOverlapUsesGeometricAreaWithoutIntegerOverflow](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Expected 1; actual 0.5600464940071106; tolerance 2.2E-05. |
| [GeometryAndRenderingAuditTests.RectangularDepthTransformsPreserveCoordinateMeaning](../../tests/UMapx.Tests/GeometryAndRenderingAuditTests.cs) | 1 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.HistogramMedianSelectsTheMiddleObservationForOddCounts](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 3 | Assert.Equal() Failure: Values differ |
| [ImagingAuditTests.NeutralFilterSettingsPreservePixels](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 7 | Expected RGBA 47,11,89,57; actual 46,10,88,57; tolerance 0. |
| [ImagingAuditTests.PixelFiltersHonorBitmapDataStrideAndLeavePaddingUntouched](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 8 | Expected RGBA 0,0,255,60; actual 0,0,0,60; tolerance 0. |
| [ImagingAuditTests.PixelwiseFiltersMatchIndependentChannelEquations](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 1 | Expected RGBA 216,11,140,125; actual 214,11,138,125; tolerance 1. |
| [ImagingAuditTests.TensorConversionsMatchPixelChannelOrder](../../tests/UMapx.Tests/ImagingAuditTests.cs) | 4 | Assert.InRange() Failure: Value not in range |
| [RemainingImagingAuditTests.ConstantColorTransferPreservesIdenticalImages](../../tests/UMapx.Tests/RemainingImagingAuditTests.cs) | 2 | Expected RGBA 77,99,123,255; actual 0,0,0,255; tolerance 1. |
| [UtilityContractAuditTests.DepthHistogramEqualizationCountsMoreThan65535PixelsWithoutOverflow](../../tests/UMapx.Tests/UtilityContractAuditTests.cs) | 2 | Assert.Equal() Failure: Values differ |
| [VideoAuditTests.MjpegFramesSurviveArbitraryTransportChunkBoundaries](../../tests/UMapx.Tests/VideoAuditTests.cs) | 3 | System.ArgumentOutOfRangeException : Index was out of range. Must be non-negative and less than or equal to the size of the collection. (Parameter 'startIndex') |
| [VideoAuditTests.MultipartBoundaryParsingHonorsMimeParameters](../../tests/UMapx.Tests/VideoAuditTests.cs) | 2 | System.ArgumentException : Invalid content type |
