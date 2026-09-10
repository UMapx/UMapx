# UMapx Mathematical Audit — September 10, 2026

The library working copy was audited against commit
`9d623dcdede20411f220dee976e62cf4c8cf5728`, project version 7.5.1.5.
The algorithms were not modified. The `UMapx.Tests` project was added to the
existing `sources/UMapx.sln`, together with reproducible reference data.

**The absence of mathematical errors cannot be confirmed: concrete
counterexamples were found. This report describes 31 issue groups covered by
43 regression tests.** Some groups combine related errors across multiple
overloads. Additional special function comparisons reveal accuracy discrepancies;
the number of failing tests must not be treated as the number of independent defects.

## Reproducible Test Results

```powershell
dotnet test sources/UMapx.sln -c Release -p:GeneratePackageOnBuild=false
```

| Suite | Total | Passed | Failed |
| --- | ---: | ---: | ---: |
| Mathematical identities (`Identity`) | 48 | 48 | 0 |
| Counterexamples (`Regression`) | 43 | 0 | 43 |
| Special function references (`Reference`) | 749 | 567 | 182 |
| **Total** | **840** | **615** | **225** |

The build succeeds. The test exit code is 1 because of failed assertions.
The tests are enabled and expect mathematically correct results. A repeat run
after finalizing the project produced the same counts. The environment was
Windows with .NET SDK 10.0.401; the tests target .NET 8, while the library
continues to target `netstandard2.0`.

The methodology, test filters, and reference generation instructions are in the
[test project README](../tests/UMapx.Tests/README.md).
The complete local test report is at
`artifacts/math-audit/test-results/math-audit-final.trx`.
The `artifacts` directory is ignored by Git; a new report can be generated with
`--logger trx`.

## Counterexamples and Causes

Values below are rounded for readability. Exact constants and tolerances are in
[MathematicalRegressionTests.cs](../tests/UMapx.Tests/MathematicalRegressionTests.cs).
Line numbers refer to the original commit. Here, `i` is the imaginary unit, and
function argument order follows the UMapx API, for example `J(x, n)`.

### 01. Missing normalization in the lower incomplete gamma function

[Special.cs, line 1931](../sources/Core/Special.cs#L1931); complex version: line 1953.

`Special.GammaP(5, 5)` returns **13.4281625748**, whereas the expected value is
**0.559506714935**. For positive arguments, the regularized function must lie
in `[0,1]`. `LowerRegGammaSeries` computes the lower incomplete gamma function
without dividing by `Gamma(s)`: the example is wrong by a factor of `Gamma(5) = 24`.
The normalization should be included in the exponential prefactor through
`-LogGamma(s)`. The definition is given in [NIST DLMF §8.2](https://dlmf.nist.gov/8.2).

The error propagates through `GammaP`, `GammaQ`, and the incomplete gamma functions
to distributions that use them, including Gamma, Erlang, and ChiSquare.
Checking only `P + Q = 1` does not detect the shared normalization error.

### 02. Incorrect denominator shift in the upper gamma continued fraction

[Special.cs, line 1982](../sources/Core/Special.cs#L1982); complex version: line 2015.

`Special.GammaQ(2, 3)` returns **0.179233416915**, whereas the correct value is
`4*exp(-3)` = **0.199148273471**. After initializing `b0 = x + 1 - s`, the
first iteration reuses `b0` instead of using `b1 = b0 + 2`. The increment of `b`
is placed at the end of the iteration. This is a separate error in the continued
fraction branch and would remain after fixing issue 01.
See [NIST DLMF §8.9](https://dlmf.nist.gov/8.9).

### 03. Quadratic equations with a zero linear coefficient

[Maths.cs, line 2081](../sources/Core/Maths.cs#L2081).

`Maths.Quadratic(1, 0, -1)` returns zero and `NaN`, although the roots are `-1`
and `1`. The stable formula uses `Math.Sign(b)`, which makes `q` zero when
`b = 0`, followed by division in `c/q`. The tests also cover `x² = 0` and
`x² + 1 = 0`. A zero `q` needs separate handling, and the sign choice must be
correct when `b = 0`. `BiQuadratic` also depends on this function.

### 04. Negative real cube roots are computed using fractional powers

[Maths.cs, line 2054](../sources/Core/Maths.cs#L2054); also line 2062.

`Maths.Cubic(0, 0, -1)`, corresponding to `x³ - 1 = 0`, returns `NaN` for
every root. The case `x³ - 3x - 2`, whose roots are `2, -1, -1`, also fails.
`Math.Pow(negative, 1/3)` does not implement the real cube root.
A sign-preserving real cube root is needed. The tests substitute the returned
roots into the original polynomial.

### 05. Biorthogonal wavelets use incorrect reconstruction filters

[WaveletPacket.cs, lines 141–146](../sources/Wavelet/WaveletPacket.cs#L141).

For `Bior13`, one level, and input `[1,0,0,0,0,0,0,0]`, the result of
`Backward(Forward(x))` with `normalized = false` is approximately
`[1.015625,-0.015625,0,-0.125,-0.015625,0.015625,0,0.125]`.
The normalized variant also fails. The discrepancy is not a uniform scale factor.

`Create(scaling, wavelet)` constructs synthesis filters by reversing the same
analysis filters. Biorthogonal Bior/CDF banks require the corresponding dual
reconstruction filters. The four separate filters are described in the
[PyWavelets documentation](https://pywavelets.readthedocs.io/en/latest/ref/wavelets.html).
This finding does not automatically apply to orthogonal Haar/Daubechies banks:
their tested variants reconstruct the signal.

### 06. Incorrect sign in complex arccotangent

[Maths.cs, line 909](../sources/Core/Maths.cs#L909).

`Maths.Actan(1 + 0i)` returns **−π/4**, whereas **π/4** is expected, consistent
with the real overload on the positive axis. The sign before the logarithm does
not match the numerator and denominator order. This also breaks inversion with
cotangent; the difference is not merely an additive integer multiple of `π`.

### 07. Complex Acosh does not follow the principal branch

[Maths.cs, line 1124](../sources/Core/Maths.cs#L1124).

`Maths.Acosh(-2 + 0i)` returns **−1.316957831 + πi** instead of the principal
value **+1.316957897 + πi**. With the standard principal square root,
`Log(z + Sqrt(z*z - 1))` does not select the principal branch throughout the
plane. The root must be chosen consistently, for example using
`sqrt(z-1)*sqrt(z+1)` with appropriate branch cut handling.
The test explicitly requires the principal branch described in
[NIST DLMF §4.37](https://dlmf.nist.gov/4.37).

### 08. The complex erf approximation fails for supported argument types

[Special.cs, line 2491](../sources/Core/Special.cs#L2491).

`Erf(-2 + 0.5i)` returns **−0.889444888 + 0.677245617i**, whereas the reference
is **−1.003502243 + 0.004740903i**. A real rational approximation is applied
directly to complex arguments without appropriate domain splitting and
reflection; the asymptotic branch also needs to account for the argument sector.
Increasing the iteration count alone is insufficient. Related `Erfc` and `Erfi`
functions are affected; even the real value `Erfi(2)` is approximately `18.42167`
instead of `18.56480`.

### 09. Lambert W converges to a different branch

[Special.cs, line 781](../sources/Core/Special.cs#L781).

`LambertW(0.1f + 0.5i, 0)` returns **−1.873011112 − 2.739237547i**, whereas
the expected value is **0.214530284 + 0.351095550i**. For this argument, the
initial approximation `L - Log(L)` sends the iterations away from the principal
branch. A small residual `w*exp(w) - z` does not establish that the requested
branch was selected. The initial approximation must account for the argument
region and branch number.

### 10. Missing terms in the Bessel Y integral formula

[Special.cs, line 3516](../sources/Core/Special.cs#L3516); complex version: line 3555.

`Special.Y(1f, 0)` returns **0.328457683**, whereas the reference is
**0.088256964**. For integer `n`, the second integrand should contain
`exp(-x*sinh(t)) * (exp(n*t) + (-1)^n*exp(-n*t))`.
The implementation uses only `exp(-x*sinh(t) - n*t)`. Even when `n = 0`, a
factor of 2 is missing. See [NIST DLMF 10.9.7](https://dlmf.nist.gov/10.9.E7).

### 11. Bessel J asymptotics: applicability and lost complex phase

[Special.cs, line 3440](../sources/Core/Special.cs#L3440); complex version: line 3482.

`J(21f, 10)` returns **0.467535526** instead of **0.148531806**. The threshold
`|x| >= 20` does not account for the function order or the size of omitted terms.
The complex version also uses `|z|` in the amplitude and inverse powers where
`z` must be retained. As a result, `J(30i, 10)` has a huge imaginary component,
although its relation to `I₁₀(30)` requires a real answer: **−1.458318100e11**.
See the dependence of the coefficients on the order and complex argument in
[NIST DLMF §10.17](https://dlmf.nist.gov/10.17).

### 12. Complex Bessel K quadrature does not resolve oscillations

[Special.cs, line 3695](../sources/Core/Special.cs#L3695).

`K(1 + 5i, 10)` returns **140853712 − 48188600i**, whereas the expected value
is **3.015242163 − 31.062215313i**. Fixed quadrature with 16 nodes does not
provide sufficient accuracy for the oscillatory integral with a large
`cosh(n*t)` factor. The argument is in the right half-plane, where the integral
representation applies; the problem is its numerical evaluation.
Quadrature error control or another stable representation is needed.

### 13. Beta overflows in intermediate gamma functions

[Special.cs, line 2904](../sources/Core/Special.cs#L2904).

`Beta(20f, 30f)` returns `NaN`, although the result **1.768188547306e−15** is
representable in `float`. Direct evaluation of `Gamma(a)*Gamma(b)/Gamma(a+b)`
overflows before cancellation. Computation should use logarithms through the
existing `LogBeta`, with appropriate handling of the parameter domain.
The error affects Beta consumers, including the regularized incomplete beta function.

### 14. The hypergeometric series overflows inside its convergence disk

[Special.cs, line 2667](../sources/Core/Special.cs#L2667).

`Hypergeom(2f, 3f, 4f, 0.9f)` returns `NaN`, whereas the expected value is
**21.7894168873**. This argument is inside `|z| < 1`, so missing analytic
continuation does not explain the problem. Factorials and Pochhammer symbols
are computed separately and overflow before division. Series terms should be
updated using the ratio of consecutive terms, with convergence monitoring.

### 15. Pade approximation accesses a negative index

[Pade.cs, line 93](../sources/Analysis/Pade.cs#L93); complex version: line 142.

`new Pade(1, 3).Compute([1,1,1/2,1/6,1/24])` throws
`IndexOutOfRangeException`. This is a valid `[1/3]` approximation of `exp(x)`:
the numerator is `[1,1/4]` and the denominator is `[1,-3/4,1/4,-1/24]`.
The system matrix accesses `taylorCoeffs[m+i-j]` with an index that can become
negative. The corresponding coefficient of the formal series should be zero.

### 16. Bilinear interpolation loses linearity at the x boundaries

[Interpolation.cs, line 166](../sources/Analysis/Interpolation.cs#L166).

For a grid with `x = y = [0,1]` and `z[i,j] = x[i]+y[j]`, the value at
`(0,0.5)` is **0** instead of **0.5**; at `(1,0.5)`, it is **1** instead of
**1.5**. The early return selects the lower node in `y` and skips interpolation
along the boundary. Even if values outside the grid are clamped, linear
interpolation along the other coordinate must be preserved on the boundary itself.

### 17. Complex variance and norm use squares without conjugation

[Matrice.cs, line 5057](../sources/Core/Matrice.cs#L5057); norm: line 8634.

The sample `Var([i,-i])` is **−2**, whereas the usual complex variance
`sum(abs(z-mean)^2)/(n-1)` is **2**. `Abs([1,i])` is **0**, although the
Euclidean norm is **sqrt(2)**. The implementation uses `(z-mean)^2` and `z^2`
instead of squared magnitudes. This is not rounding error: nonzero components
cancel out. If a non-conjugated pseudovariance is intended, it should be exposed
and named explicitly; the current result is not an ordinary variance or norm.

### 18. Binomial probability becomes NaN when p = 1

[Binomial.cs, line 204](../sources/Distribution/Binomial.cs#L204).

`new Binomial(5,1).Function(5)` and `new Binomial(0,1).Function(0)` return
`NaN` instead of **1**. The log probability calculation produces `0*log(0)` in
the `(n-k)*log(1-p)` term. A similar safeguard for `k*log(p)` already exists.
Degenerate distributions and zero multipliers need correct handling.

### 19. The Poisson PMF overflows for moderate parameters

[Poisson.cs, line 144](../sources/Distribution/Poisson.cs#L144).

`new Poisson(100).Function(100)` returns `NaN`, whereas the reference is
**0.0398609968091**. Direct evaluation of `exp(-lambda)*lambda^k/k!` overflows
the power and factorial, even though the probability itself is moderate.
Suitable alternatives include the log probability using `LogGamma(k+1)` or a
stable recurrence.

### 20. The binomial median is computed using a mode formula

[Binomial.cs, line 149](../sources/Distribution/Binomial.cs#L149).

`new Binomial(2,0.7f).Median` returns **2**, whereas the correct median is **1**.
For exact `p = 0.7`, the probabilities are `P(X <= 1) = 0.51` and
`P(X >= 2) = 0.49`; rounding `p` to `float` does not change the conclusion.
The formula `floor((n+1)*p)` does not define the median. Both probability
inequalities that characterize a median must be checked.

### 21. The Poisson median is incorrect between log(2) and 1

[Poisson.cs, line 105](../sources/Distribution/Poisson.cs#L105).

`new Poisson(0.9f).Median` returns **0**, whereas the median is **1**:
`P(X=0) = exp(-0.9) ≈ 0.40657 < 0.5`. The branch returning zero for all
`lambda < 1` covers too wide a range. An asymptotic median estimate needs
verification against probabilities, particularly for small parameters.

### 22. IIR.Stability checks a different characteristic polynomial

[IIR.cs, line 310](../sources/Response/IIR.cs#L310).

For `b = [1]` and `a = [0.5,-0.75]`, `Stability` returns **true**.
However, `Reaction` divides by `1-a[0]`, so the pole is **−1.5** and the
impulse response begins `[2,-3,4.5,-6.75,...]`: its amplitude grows.

The stability check constructs coefficients `[1,-a[0],-a[1],...]`, whereas
the denominator of the implemented recurrence corresponds to
`[1-a[0],-a[1],...]`. The polynomial must agree with `Reaction`, including
coefficient order and normalization.

### 23. Even-length Normal and Confined windows are asymmetric

[Normal.cs, line 69](../sources/Window/Normal.cs#L69) and
[Confined.cs, line 72](../sources/Window/Confined.cs#L72).

At length 8, the differences between the endpoint values are approximately
**0.198866** and **0.294383**, respectively, instead of zero.
In `float a = (frameSize - 1) / 2`, division is performed using integers,
placing the center at 3 instead of 3.5. Division by `2f` is needed.
Very short windows also need explicitly defined limiting values.

### 24. XYZ clips a valid white-point coordinate

[XYZ.cs, line 36](../sources/Colorspace/XYZ.cs#L36); setter: line 77.

The declared white point `(0.9505,1,1.089)` actually stores `Z = 1`.
Converting RGB white through XYZ returns **(255,254,244)** instead of
**(255,255,255)**. XYZ coordinates are not individually bounded by one,
even with the normalization used here. The incorrect restriction should be
removed, and the white point, matrices, and dependent LAB conversions checked
for consistency.

### 25. RGB ↔ RYB conversion loses white

[RYB.cs, line 202](../sources/Colorspace/RYB.cs#L202); inverse conversion: line 265.

`RYB.FromRGB(255,255,255).ToRGB` returns **(0,0,0)**. After removing the
achromatic component, all components are zero and normalization computes
`0/0`; the safeguards are commented out. A zero chromatic remainder needs
handling in both conversion directions.

### 26. Unstable formulas for real hyperbolic functions

[Maths.cs, line 959](../sources/Core/Maths.cs#L959); `Asinh`: line 995.

`Tanh(100f)` returns `NaN` instead of **1** because the `Sinh/Cosh` ratio
produces `Infinity/Infinity`. `Asinh(-10000f)` returns **−Infinity** instead
of **−9.903487555** because of catastrophic cancellation in
`x + sqrt(x*x+1)`. Numerically stable formulas must account for the sign and
magnitude of the argument.

### 27. Complex division is unstable under scaling

[Complex32.cs, line 250](../sources/Core/Complex32.cs#L250).

For `z = 1e20f + 1e20f*i`, the expression `z/z` returns `NaN`, although the
answer is **1**. The test with scale `1e-20f` also fails. Directly squaring
components causes overflow or loss of precision in the subnormal range before
the final division. A scaled complex division algorithm is needed.

### 28. Chebyshev polynomials fail at valid arguments

[Special.cs, line 33](../sources/Core/Special.cs#L33); `ChebyshevU`: line 53.

`ChebyshevU(1f,2)` returns `NaN`, although **U₂(1) = 3**: the trigonometric
formula produces `0/0`. `ChebyshevT(2f,2)` also returns `NaN`, although
**T₂(2) = 7**. The polynomial is defined outside `[-1,1]`, where the real
`acos` is undefined. Evaluating the polynomials by recurrence removes these
restrictions.

### 29. Rising factorial with a zero base

[Special.cs, line 2180](../sources/Core/Special.cs#L2180).

`FactorialUp(0,2)` returns **1**, although `(0)₂ = 0*1 = 0`.
The `n == 0` branch is incorrectly treated as equivalent to zero order.
The value 1 is required for order `k = 0`; for a positive integer order and
a zero base, the result is zero.

### 30. Negative real base with a complex exponent

[Maths.cs, line 594](../sources/Core/Maths.cs#L594).

`Maths.Pow(-1f, new Complex32(0.5f,0))` returns `NaN` instead of the principal
value **i**. The overload uses the real logarithm of a negative number.
The base must be converted to a complex representation before taking the
logarithm, using a branch consistent with the fully complex overload.

### 31. Incorrect digamma coefficient in InverseChiSquare entropy

[InverseChiSquare.cs, line 153](../sources/Distribution/InverseChiSquare.cs#L153).

`new InverseChiSquare(2).Entropy` returns **0.306852967**, whereas the correct
differential entropy is **1.461284149**. Using the equivalent InverseGamma
distribution with `alpha = v/2` and `beta = 1/2`, the formula is
`alpha + log(beta) + LogGamma(alpha) - (1+alpha)*DiGamma(alpha)`.
The implementation uses `1-alpha` instead of `1+alpha`. At `v = 2`, the
reference reduces to `1 - log(2) + 2*EulerGamma`.

## Additional Special Function Comparisons

Reference values were computed independently with mpmath 1.3.0 at 60 decimal
digits of precision. Inputs were first rounded to `float`; the tests use a
mixed tolerance of `2e-4 + 2e-4*abs(expected)`. For complex values, the error
is measured by its magnitude. This tolerance is an audit criterion, not an
accuracy guarantee found in the library.

| Function | Passed | Failed |
| --- | ---: | ---: |
| Gamma / LogGamma / DiGamma / TriGamma | 19 / 15 / 19 / 19 | 0 / 0 / 0 / 0 |
| Ci / Si / Ei / Li / Zeta | 16 / 19 / 19 / 15 / 18 | 0 / 0 / 0 / 0 / 0 |
| Erf / Erfc / Erfi | 17 / 17 / 5 | 2 / 2 / 11 |
| Fresnelc / Fresnels | 20 / 21 | 2 / 1 |
| J / Y / I / K | 43 / 12 / 45 / 41 | 22 / 53 / 20 / 24 |
| Struve H / L | 56 / 45 | 9 / 20 |
| LambertW | 24 | 2 |
| Beta / BetaIncomplete / BetaIncompleteRegularized | 3 / 18 / 15 | 1 / 2 / 5 |
| Hypergeom, both overloads | 26 | 6 |

In particular, discrepancies in `I`, `H`, `L`, and some `K` cases require
further analysis of the approximation regions. The 31 issue groups above are
not a complete breakdown of all 182 failing reference comparisons by cause.
Fresnel errors remain near the switching point `x = 6`: for example,
`Fresnels(6)` returns approximately `0.639870465` instead of `0.638459189`,
even after accounting for the UMapx normalization convention.

## Confirmed Cases and Audit Limitations

The maintained test suite successfully checks:

- FFT against an independent DFT at lengths 1, 2, 3, 4, 5, 7, 8, 15, 16, 17,
  and 31; inversion of a rectangular 3×5 complex FFT along both directions.
- Inversion of direct and fast cosine/sine/Hartley/Walsh–Hadamard/Chebyshev
  transforms at selected powers of two.
- Matrix reconstruction for SVD, QR, LU, LDU, LQ, QL, RQ, Polar, Hessenberg,
  Schur, EVD, Cholesky, LDL, and UDL on small deterministic matrices;
  rectangular SVD and QR, and orthonormality of QR columns.
- Exact integration of a cubic polynomial by Simpson's method with even and
  odd subinterval counts; reproduction of a quadratic polynomial by
  Lagrange/Newton/Barycentric interpolation.
- Signal reconstruction by Haar/D2/D4 at three levels, the FIR impulse
  response, and the definitions of four basic distances.

The exploratory runs also checked other decompositions, wavelets,
distributions, differentiation, and ODEs. Individual passing examples do not
replace maintained tests or provide a general guarantee.

The structural review covered a library of 426 C# files. Numerical checks
focused on the mathematical core, analysis, decompositions, transforms,
wavelets, distributions, filters, windows, and selected color conversions.
**This audit did not verify all 426 files line by line, exhaustively test all
parameters and conditioning, cover every complex branch cut, or fully audit
image processing, video, and visualization algorithms.** Passing tests confirm
the tested cases, not correctness throughout the entire domain.

The following observations were deliberately excluded from the list of proven
formula errors:

- Unpivoted LU does not reconstruct `[[0,1],[1,0]]`, but that matrix does not
  admit a standard unpivoted LU with a unit diagonal in L. The contract for
  a zero pivot must specify either pivoting or explicit rejection. Testing
  ordinary random matrices alone cannot establish this behavior.
- The Chord method may fail to converge for particular functions and initial
  values. Without examining convergence conditions, that is not proof of an
  incorrect formula.
- Legendre filters are explicitly nonorthogonal. Failure of a simple inverse
  transform to reproduce the input is not automatically equivalent to the
  Bior/CDF defect.
- UMapx Fresnel functions use `cos(t²)` and `sin(t²)`. The difference from
  the normalized mpmath convention is not itself an error.
- One complex `LogGamma` example is excluded from the reference tests because
  the branch contract is unspecified; the discrepancy is not hidden by a
  larger tolerance.

## Suggested Fix Order

First address both incomplete gamma branches, the polynomial root solvers,
Bior/CDF reconstruction filters, and the IIR stability polynomial. Then
address basic complex arithmetic and the branches and approximations of
special functions, which affect many consumers. Follow with the remaining
local formulas, distribution boundaries, interpolation, windows, and colors.

After fixing each cause, run the corresponding regressions and dependent
reference tests while preserving the correct mathematical expectations.
Increasing tolerances to match existing errors or disabling failing cases
does not establish mathematical quality.
