# Matrix and probability-distribution repairs — September 10, 2026

**B03 and B04 are closed: all 113 assigned failures now pass.** The complete
Windows Release audit contains **14,112 cases: 13,914 passed, 198 failed, none
skipped**. All 12,130 previous test IDs remain. There are **no regressions, no
removed tests, and 1,982 added passing cases**. Existing assertions, fixtures,
and tolerances were retained.

| Scope | Previous failures | Current failures | Resolved |
| --- | ---: | ---: | ---: |
| B03: matrix indexing, complex statistics, distances | 46 | 0 | 46 |
| B04: probability distributions | 67 | 0 | 67 |
| B05–B12 | 198 | 198 | 0 |
| **Total** | **311** | **198** | **113** |

These are test-case counts, not counts of independent defects. The resolved ID
set is exactly the union of B03/B04 assignments in the arithmetic snapshot.
No additional downstream failures disappeared. All 3,749 dedicated
special-function checks and the 1,626 arithmetic/number-theory repair cases
continue to pass. The remaining decomposition failures are still in B06.

## Matrix changes

Production changes are in [Matrice.cs](../sources/Core/Matrice.cs), including its
private statistics helpers, and
[SokalSneath.cs](../sources/Distance/SokalSneath.cs).

- All four left-diagonal products now implement `diag(v) * A`, scaling rows and
  requiring a diagonal of the correct length. Right-diagonal products retain
  their column scaling. Inverse diagonal products preserve the existing convention
  that a zero diagonal entry produces a zero row/column.
- Matrix shifts use the correct dimension for each axis. Matrix and vector
  shifts widen the index subtraction before wrapping, including `int.MinValue`.
- The four mesh overloads evaluate `f(x[i], y[j])`. Row/column swaps traverse
  the correct axis lengths in rectangular matrices.
- Merge clips the requested patch to the destination intersection, including
  negative offsets, and leaves the inputs intact. A patch already at the requested
  size is copied directly. Actual resizing still uses the existing interpolation
  code; its known defects remain assigned to B05.
- Complex variance and norms use squared magnitudes. Covariance uses the
  centered Hermitian product with the first column conjugated. Double accumulation
  prevents premature float overflow, and standard deviation takes the square
  root before converting to float. Covariance pairs are computed together to
  preserve conjugate symmetry.
- Sokal–Sneath uses `2*(tf+ft)/(tt+2*(tf+ft))`, with widened denominator arithmetic.
  This is the [SciPy documented contingency formula](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.sokalsneath.html).
  The library's existing zero result for two all-zero vectors is retained.

## Distribution changes

Seventeen distribution classes were corrected. The internal
[DistributionNumerics.cs](../sources/Distribution/DistributionNumerics.cs) helpers
provide stable logarithmic calculations and transformed quadrature.
[Special.cs](../sources/Core/Special.cs)
exposes the already audited double kernels internally, avoiding intermediate
rounding through the public float APIs. The underlying special-function algorithms
were not changed.

| Distribution | Corrected behavior |
| --- | --- |
| PowerNormal, PowerLognormal | Consistent survival-power PDF/CDF; stable tails for powers below one and large arguments. |
| Binomial | Point masses at `p=0/1`, zero trials, log-space PMF, double-precision incomplete-beta CDF, lower median by CDF search. Integer support checks avoid rounding an unrepresentable trial count onto a neighboring float input. |
| Poisson | Log-space PMF, incomplete-gamma CDF, lower-median search, exact integer mode test, entropy from a probability recurrence or large-rate expansion. |
| Birnbaum–Saunders | Median equals location plus scale; mode is the positive stationary root; entropy is integrated through the normal-variable transformation. |
| BetaPrime, FisherSnedecor | Correct excess-kurtosis formula and its parameter dependence. |
| InverseChiSquare, Kumaraswamy, Levy, Wigner | Correct differential entropy, including scale factors and digamma terms. |
| FisherZ | Stable log-beta PDF and incomplete-beta CDF, including large finite arguments. |
| Gompertz | Mean and variance from transformed quadrature; centered variance avoids subtracting nearly equal second moments. |
| FoldedNormal | Mode zero when `abs(mu) <= sigma`; otherwise solves the density stationarity equation. |
| ChiSquare | Median by incomplete-gamma CDF inversion. |
| Trapezoidal | Mean and variance from the weighted component moments of both power ramps and the middle region; zero-width pieces and translations are supported. |
| TukeyLambda | Moment-existence thresholds; stable variance near zero and the exact logistic limit at zero. |

PowerNormal follows `F(x)=1-Phi(-x)^p`, as specified by
[NIST](https://www.itl.nist.gov/div898/handbook/eda/section3/eda366d.htm).
PowerLognormal applies the same law to `log(x)/sigma` with the PDF Jacobian
`1/(x*sigma)`, following [NIST's power-lognormal definition](https://www.itl.nist.gov/div898/handbook/eda/section3/eda366e.htm).
For integer powers this describes the minimum of independent normal/lognormal
variables. The previous CDF based on `Phi(x)^p` described a different law.

## Contracts and compatibility

Public method/property signatures, parameter defaults, and result types are
unchanged. Corrected numerical results intentionally differ from the old values.

| API | Explicit convention |
| --- | --- |
| `Matrice.Var(v)`, `Cov(v)` for complex values | Sample variance `sum(abs(v-mean)^2)/(n-1)`, returned as a real `Complex32`. Fewer than two observations give NaN in the real component. |
| Complex `Var(x,y)` and `StnDev(x,y)` | Squared differences normalized by `n-1`, and their square root. These overloads retain their existing difference-statistic meaning; they are not cross-covariance. |
| `Matrice.Cov(matrix)` | Columns are variables; `C[j,k] = sum(conj(x[j]-mean[j])*(x[k]-mean[k]))/(n-1)`. |
| Complex `Abs` | Euclidean norm, or squared norm when requested; matrix overload returns row norms. Imaginary component is zero. |
| `Swap` | Horizontal swaps rows; Vertical swaps columns; Both performs both swaps. Existing axis meanings are retained. |
| `Merge` | Returns a destination clone with the overlapping patch inserted. Negative requested sizes are rejected. |
| Binomial/Poisson median | Lower median, rounded to binary32. For symmetric binomial laws the lower median is selected exactly; the upper endpoint of a median interval can also be a valid mathematical median. |
| ChiSquare median | Numerically inverts the CDF; no longer returns the Wilson–Hilferty approximation. |
| Differential entropy | Natural logarithms in the repaired entropy getters. Negative differential entropy is valid. |

Unsupported getters remain explicit and are not counted as implemented by this
repair. Undefined/divergent moments retain explicit NaN/infinity results.

## Validation and independent references

The new tests add **1,982 passing cases**:

- [MatrixRepairTests](../tests/UMapx.Tests/MatrixRepairTests.cs): **316 cases**.
  All real/complex diagonal combinations, rectangular/empty shapes, zero diagonal
  entries, extreme shifts, clipped and empty patches, nonseparable meshes,
  swap involutions, Hermitian moments, covariance outer products/positivity,
  and exact boolean contingency counts. Matrix expectations use explicit index
  maps and independent `System.Numerics.Complex` double arithmetic.
- [DistributionRepairTests](../tests/UMapx.Tests/DistributionRepairTests.cs):
  **1,666 cases**, including **1,611 mpmath references at 70 decimal digits**
  and 55 mode, scaling, support, NaN, and endpoint checks. The
  [generator](../tests/UMapx.Tests/Data/generate_distribution_repair_reference.py)
  evaluates the exact binary32 inputs. Moments and entropies include direct
  PDF/survival integration, discrete probability summation, and high-precision
  CDF inversion independent of the production algorithms.

New finite distribution references use `2*float.Epsilon + 2e-5*abs(expected)`;
small nonzero probabilities cannot pass merely by returning zero. Discrete
medians compare exactly, and NaN/infinity expectations are explicit. Additional
mode checks use the stationary equation and independently evaluated densities.
Existing test budgets were not loosened.

The extended checks also reproduced and fixed vector-shift integer overflow,
premature overflow in complex deviation/covariance, and a rounded integer-support
comparison in a degenerate binomial law. The symmetric binomial median convention
was made explicit. The fixture generator was rerun and produced byte-identical
JSON (SHA256 `d9699e2b273fe7f7d78fd2b6f26f98f9c944be5b20e5f5e70cef89f61facbc91`).

## Reproduction and evidence

```powershell
dotnet test tests/UMapx.Tests -c Release --no-restore -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~MatrixRepairTests|FullyQualifiedName~DistributionRepairTests|FullyQualifiedName~DistributionReferenceTests|FullyQualifiedName~DistributionShapeAuditTests|FullyQualifiedName~MatrixStructureAuditTests|FullyQualifiedName~DistanceAuditTests"
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/matrix-distribution-repair/verification
```

The recorded complete run is `artifacts/math-audit/matrix-distribution-repair/final`.
It exits with **1** because the other 198 failures remain enabled.

- [Run summary and source/input hashes](audit-matrix-distribution/summary.json).
- [Every-test-ID comparison with the 12,130-case baseline](audit-matrix-distribution/comparison.json).
- [Remaining failures](audit-matrix-distribution/failures.json) and [failing families](audit-matrix-distribution/failures.md).
- [Source inventory and coverage](audit-matrix-distribution/source-inventory.md) and [uncovered methods](audit-matrix-distribution/uncovered-methods.json).
- [Current block ownership](audit-matrix-distribution/repair-blocks.json) and [remaining work order](remaining-repair-blocks-2026-09-10.md).
- [Historical B03/B04 assignments](audit-arithmetic/repair-blocks.json).

Coverage is **83.29% of executable lines** (28,129/33,769) and **76.42% of branches**
(10,613/13,886). All 436 production source hashes match this audit snapshot.
The 198 remaining IDs each have exactly one block owner. Source, test, tool, and
report text was checked for Cyrillic characters.

## Limits and next work

Closing these blocks resolves their known audit failures; it does not prove
correctness over every input. High-precision fixtures cover selected parameter
ranges, including very small probabilities, not every representable float or
every extreme combination. Numerical entropy/moment getters now perform
quadrature; their accuracy and runtime have not been benchmarked over the entire
parameter domain. Poisson's large-rate entropy uses an asymptotic expansion.

The remaining eight block counts are unchanged: B05 15, B06 15, B07 57, B08 18,
B09 23, B10 21, B11 44, B12 5. **B05 is the recommended next block** because its
interpolation, resampling, and local-filter fixes support subsequent work.
The known zero-matrix Schur termination failure remains the first subtask for B06;
its subprocess timeout stays enabled.
