# Special-function repairs — September 10, 2026

This is the historical special-function snapshot. The subsequent
[arithmetic and number-theory repair report](arithmetic-repair-2026-09-10.md)
records 139 more resolved failures and the current 311 remaining failures.

All failures in the dedicated special-function suites are resolved. The final
Windows Release run has **10,504 cases: 10,054 passed, 450 failed, none skipped**.
The remaining failures stay enabled, and the complete test command exits with 1.

Every one of the previous 8,408 test IDs is still present. Comparing those IDs
shows **423 failures resolved, no passing tests regressed, and no tests removed**.
All **2,096 added cases pass**. The original acceptance tolerances were retained.

| Existing test suite | Failures before | Failures after | Resolved |
| --- | ---: | ---: | ---: |
| Special-function references | 358 | 0 | 358 |
| Mathematical regressions | 43 | 30 | 13 |
| Distribution references | 91 | 50 | 41 |
| Distribution shape checks | 15 | 11 | 4 |
| Additional transforms | 13 | 6 | 7 |
| All other suites | 353 | 353 | 0 |
| **Total on the original cases** | **873** | **450** | **423** |

The dedicated special-function suites contain **3,749 passing cases**: 1,654
original reference records, 2,073 additional references, and 22 boundary/identity
cases. The formerly excluded complex LogGamma record is now included because
its continuation is explicitly defined. The 13 repaired special-function
counterexamples in `MathematicalRegressionTests` also pass.

## What changed

Production changes are confined to `sources/Core/Special.cs` and six new partial
class files in the same directory. Public method signatures and return types
are preserved. Intermediate calculations use `double` and `System.Numerics.Complex`;
conversion to public single-precision types happens at the boundary.

- **Gamma and beta:** restore incomplete-gamma normalization and the next
  denominator in the continued fraction; evaluate complementary tails directly;
  retain tiny upper tails for shapes approaching zero; use logarithmic gamma
  ratios for beta and a symmetric continued fraction for real incomplete beta.
  Gamma, log-gamma, digamma, trigamma, and zeta use double-precision kernels.
- **Factorials and polynomials:** evaluate integer rising/falling factorials and
  binomial coefficients without gamma poles; retain factorials that fit in the
  public `double` result. Chebyshev endpoint/outside-interval values and negative
  integer orders are supported. Orthogonal-polynomial recurrences are iterative.
  Euler polynomials no longer divide by zero at `x = 0.5`.
- **Hypergeometric functions:** accumulate terms through ratios, factor removable
  polynomial zeros near one, and continue 2F1/1F1/0F1 using their differential
  equations when the series at zero is unsuitable. Existing NaN sentinels for
  absent hypergeometric parameters remain supported.
- **Error and integral functions:** replace the real-only approximation previously
  used for complex Erf. Faddeeva uses a decaying integral in the upper half-plane
  and reflection; large-argument expansions stop near their smallest term.
  Erfc and Q retain small tails, and inverse Q avoids forming a rounded `1-2*p`.
  Dawson, Erfi, Fresnel, Gerf, Ei, Ci, Si, Li, and Owen use the corrected kernels.
  Fresnel S has a separate series near zero.
- **Lambert W:** use branch-aware initial values and refinement, including near
  the branch point and at large arguments. The square super-root returns its
  value at one and rejects non-real results in its real overload.
- **Bessel and Struve:** use series, order-zero/order-one asymptotics, and stable
  recurrence directions. Complex Y uses its connection to J and K to avoid
  amplification of a dominant solution. Complex phase is retained. Struve H/L
  use a smooth integral after an endpoint substitution, including their integer
  order continuation.
- **Other special functions:** stable Erlang blocking recurrence, logistic and
  smoothed Heaviside functions, Gudermannian limits, exact integer Fibonacci/Lucas
  evaluation, dyadic Rademacher zeros, and the harmonic number at zero.

## Contracts and compatibility

- Real `LogGamma(x)` means `log(abs(Gamma(x)))`. Complex LogGamma is the analytic
  continuation cut along the negative real axis, using the upper side on the cut;
  its imaginary part is not reduced modulo `2*pi`.
- `Gerf(x,n)` means `n!/sqrt(pi)` times the integral from zero to `x` of
  `exp(-t^n)`, with its entire continuation for complex `x` and integer `n >= 0`.
- Fresnel functions retain the library's unnormalized integrands `cos(t*t)` and
  `sin(t*t)`. H/L remain Struve functions, not Hankel functions.
- Real incomplete beta accepts finite positive shapes and `0 <= x <= 1`.
  Out-of-domain values now return NaN instead of being clamped to endpoints.
  Inverse Erf/Q likewise return NaN outside their real domains and the correct
  infinities at the endpoints. Real 2F1 returns NaN when continuation is non-real.
- Integer factorial identities, negative integer Bessel/Struve orders, and
  selected removable singularities now return their mathematical values.
- Fibonacci returns exact Int32 results for `-46 <= n <= 46`; Lucas does so for
  `-44 <= n <= 44`. Outside those ranges they throw `ArgumentOutOfRangeException`
  instead of returning an overflowed or rounded integer.

## Independent validation and limits

The new JSON references are generated with **mpmath 1.3.0 at 80 decimal digits**,
after rounding inputs to binary32. The generator is checked in beside the data.
It uses mpmath's 1F1 representation for incomplete gamma with a negative-real-part
complex endpoint, avoiding an endpoint-ordering recursion problem in mpmath 1.3.0.

The additional grid includes Bessel/Struve orders from -5 through 50, real
arguments through 100, complex arguments in both half-planes, both sides of
algorithm switches, gamma shapes through 500, large beta parameters, Lambert W
branches, tails down to subnormal floats, and polynomial degrees through 50.
References for new cases use a relative budget of `2e-5`, with an absolute budget
of `1e-7` normally and `1e-44` for nonzero reference magnitudes below `1e-4`.
Boundary tests assert exact identities, endpoint values, invalid domains, and
overflow behavior separately. Existing reference tolerances were not loosened.

Special-function execution coverage is **1,232 / 1,297 lines (94.99%)**. Across
the whole library it is **27,991 / 33,720 lines (83.01%)** and
**10,272 / 13,574 branches (75.67%)**. These are execution measures, not a
percentage of mathematically proven algorithms. The source denominator changed
because obsolete implementations were replaced.

This validation does not establish correctness for every float/complex input,
every singularity, arbitrary orders, or arbitrarily ill-conditioned parameter
combinations. Finite reference grids exclude poles, non-real results for real
APIs, and results outside the relevant output range. Adaptive/asymptotic kernels
and their extreme-argument limits remain candidates for further specialized
accuracy and performance work. The uncovered lines are retained in the inventory.

## Remaining failures

The [repair-block plan](remaining-repair-blocks-2026-09-10.md) tracks the original
12 blocks; B01/B02 are now closed and ten blocks remain open.

The historical [failure list](audit-special-functions/failures.md) and
[individual records](audit-special-functions/failures.json) contain all 450
remaining failures. These are failed test cases, not 450 independent causes.

The distribution suites retain 61 failures, including PowerNormal/PowerLognormal
CDFs and moments/entropy formulas in several other distributions. Their
special-function dependencies are fixed, but those failures do not disappear.
The original regression suite also retains separate distribution counterexamples.

Four Hankel matrix checks at order 5 still fail. Its Newton search starts near
`10.2101765`, steps to `6.1687536`, then to `-4.5240078`, and approaches the zero
at the origin instead of the first positive zero. The resulting normalization
degenerates. Repairing `HankelTransform.BesselZeroJ` requires a bracketed positive
root search; it is separate from this special-function change. Two other
additional-transform failures concern the sign of the cone-shaped kernel.

Other outstanding areas include number theory, matrix operations and
decompositions, wavelets, windows, response filters, geometry, imaging, and
video/contracts. The [previous repair register](math-audit-expanded-2026-09-10.md)
remains the historical baseline; use the current individual failure list to
determine which counterexamples still reproduce.

## Evidence and reproduction

- [Final summary, source digest, and raw-evidence hashes](audit-special-functions/summary.json).
- [Comparison by test ID, including all 423 resolved records](audit-special-functions/comparison.json).
- [Source inventory](audit-special-functions/source-inventory.md) and
  [JSON with source hashes and uncovered lines](audit-special-functions/source-inventory.json).
- [Unexecuted methods](audit-special-functions/uncovered-methods.json).
- Historical evidence in `docs/audit` is preserved.

The production source is identified by the file hashes and source digest, not
solely by the base commit: these repairs were tested in the working tree.
Raw final evidence is in `artifacts/math-audit/special-repair/final/` (Git-ignored).

```powershell
dotnet test tests/UMapx.Tests -c Release -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~SpecialFunction"
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/special-repair/final
python -X utf8 tools/compare_math_audits.py --baseline artifacts/math-audit/final/full-audit.trx --current artifacts/math-audit/special-repair/final/full-audit.trx --output docs/audit-special-functions/comparison.json
```

The complete audit is expected to exit with 1 until the remaining failures are
fixed. The dedicated special-function run exits with 0.

## Mathematical references

The implementations use standard identities and representations documented in
NIST's Digital Library of Mathematical Functions:

- [Incomplete gamma series](https://dlmf.nist.gov/8.7) and
  [continued fractions](https://dlmf.nist.gov/8.9).
- [Error-function integral representations](https://dlmf.nist.gov/7.7).
- [Bessel power series](https://dlmf.nist.gov/10.8),
  [large-argument expansions](https://dlmf.nist.gov/10.17), and
  [modified Bessel series](https://dlmf.nist.gov/10.31).
- [Struve integral representations](https://dlmf.nist.gov/11.5).
- [Gauss hypergeometric differential equation](https://dlmf.nist.gov/15.10) and
  [Kummer's equation](https://dlmf.nist.gov/13.2).

Reference generation uses the independent
[mpmath special-function implementation](https://mpmath.org/doc/1.3.0/functions/index.html).
