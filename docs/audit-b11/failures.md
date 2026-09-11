# Failed audit checks

Generated from the recorded TRX. These are failing test families, not independent root causes.
The explanatory repair register is in [the expanded report](../math-audit-expanded-2026-09-10.md).

| Test family | Failed cases | Example |
| --- | ---: | --- |
| [VideoAuditTests.MjpegFramesSurviveArbitraryTransportChunkBoundaries](../../tests/UMapx.Tests/VideoAuditTests.cs) | 3 | System.ArgumentOutOfRangeException : Index was out of range. Must be non-negative and less than or equal to the size of the collection. (Parameter 'startIndex') |
| [VideoAuditTests.MultipartBoundaryParsingHonorsMimeParameters](../../tests/UMapx.Tests/VideoAuditTests.cs) | 2 | System.ArgumentException : Invalid content type |
