# Arithmetic and number-theory repairs — September 10, 2026

**B01 and B02 are closed: all 139 assigned failures now pass.** The complete
Windows Release run contains **12,130 cases: 11,819 passed, 311 failed, none skipped**.
All 10,504 previous test IDs remain. There are **no regressions, no removed tests,
and 1,626 added passing cases**. Existing test assertions and tolerances were retained.

| Scope | Previous failures | Current failures | Resolved |
| --- | ---: | ---: | ---: |
| B01: real and complex arithmetic | 31 | 0 | 31 |
| B02: exact integer arithmetic and number theory | 108 | 0 | 108 |
| Other blocks | 311 | 311 | 0 |
| **Total** | **450** | **311** | **139** |

These are test-case counts, not counts of independent defects. The resolved ID
set is exactly the union of B01/B02 assignments in the previous planning snapshot.
No additional downstream failures disappeared in this run. All 3,749 dedicated
special-function checks continue to pass.

## Arithmetic changes

Production changes are in [Maths.cs](../sources/Core/Maths.cs),
[Complex32.cs](../sources/Core/Complex32.cs), and the private helpers in
[Maths.Arithmetic.cs](../sources/Core/Maths.Arithmetic.cs). Public declarations,
parameter types/defaults, and result types are unchanged.

- Complex division and multiplication use double-precision intermediate products.
  This fixes `z/z` at tiny/large scales and avoids intermediate `Infinity-Infinity`
  when the resulting component is representable. Magnitude, squared magnitude,
  and logarithm likewise avoid premature float overflow or underflow.
- Real hyperbolic functions use stable evaluation, small-argument limits, explicit
  real domains, and odd symmetry. Large negative Asinh no longer suffers
  subtractive cancellation; Tanh/Ctanh approach their finite limits. Reciprocal
  functions retain representable subnormal results.
- Complex Acosh uses the product of the two principal square roots and a stable
  `log(1+z)` near one. Asinh handles the left half-plane and the negative imaginary
  cut. Actan uses a stable principal reciprocal-atan formula, including very large
  arguments. A negative real base accepts a complex exponent via the complex logarithm.
- Complex Tanh/Ctanh/Sech/Cosch evaluate their intermediate hyperbolic factors in
  double precision. They remain accurate near poles and zero and saturate only
  where the omitted term is below binary32's representable range.
- Quadratic roots handle zero linear/constant coefficients, complex conjugate
  roots, and coefficient scaling without unstable subtraction. The cubic solver
  uses real cube roots, double intermediates, exact zero roots, and real-root
  refinement. When one root is much larger than the others, it recovers the
  smaller pair from their sum and product instead of subtracting large Cardano terms.

## Number-theory changes

- Primality is independent of factorization. Miller–Rabin uses the first twelve
  prime bases, sufficient throughout the signed 64-bit domain by the bound in
  [Sorenson and Webster, Theorem 1.1](https://arxiv.org/html/1509.00864v1#S1.Thm1).
  Inputs below two return false immediately, including the former `IsPrime(1)` hang.
- Factorization retries failed Pollard–Brent cycles, uses overflow-safe modular
  multiplication, and recursively splits composite factors. Totients and radicals
  use the resulting distinct primes with exact integer arithmetic.
- GCD and extended Euclid return a nonnegative GCD with consistent Bezout
  coefficients. Intermediate extended-Euclid arithmetic uses BigInteger.
  LCM divides before multiplying and detects unrepresentable results.
- Both modular-exponentiation algorithms handle exponent zero, signed moduli,
  negative bases, and products exceeding Int64. Modular inverses and remainders
  also handle the minimum signed integer correctly.
- Coprime search returns the first valid candidate at or above its starting
  bound. Impossible searches fail explicitly rather than wrapping around.
- Base conversion, decimal digit vectors, and digit lengths use integer
  division and checked accumulation. Zero has one digit. Floating-point powers
  no longer corrupt large decimal values.
- Segmented sieve bounds remain Int64 until an actual output prime is emitted;
  the exclusive upper bound can therefore represent `int.MaxValue + 1`.

## Contracts and compatibility

Corrected results intentionally differ from the old erroneous values. The
following conventions are explicit in the repair tests:

| API | Defined behavior |
| --- | --- |
| `Gcd`, `Euclidean` | Nonnegative GCD, including zero for `(0,0)`. A positive GCD outside the return type throws `OverflowException`. |
| `Lcm` | Nonnegative; zero if either argument is zero. Unrepresentable results throw `OverflowException`. |
| `Mod`, `ModPow`, `ModInv` | Use the magnitude of a nonzero modulus. Remainders/inverses are normalized. Modulus zero throws `DivideByZeroException`. |
| `ModPow` | Exponent must be nonnegative; otherwise `ArgumentOutOfRangeException`. Exponent zero returns `1 mod abs(p)`, including zero for modulus one. |
| `ModInv` | Keeps the existing zero sentinel when no inverse exists. |
| `Coprime` | `increment` is the inclusive starting value, not a stride. Exhausting the representable search range throws `OverflowException`. |
| `Itf` | Positive integers only. Default includes prime multiplicities; `onlyPrimes=true` returns distinct primes. Results are sorted. `Itf(1)` is empty. |
| `Pollard` | A proper divisor for a composite; the input itself for a prime or one. |
| `Etf`, `Radical` | Positive integers only; both return one at one. Nonpositive factorization/totient/radical inputs throw `ArgumentOutOfRangeException`. |
| `Decimal2Base`, `Base2Decimal` | Least-significant digit first, preserving the existing conversion pair's order. Negative source integers are rejected. |
| `Numeral2Vector`, `Vector2Numeral` | Most-significant digit first, matching normal decimal notation and the existing vector-to-numeral contract. |
| Digit arrays | Radix at least two, valid nonnegative digits; null arrays are rejected. Empty arrays decode to zero; output overflow throws `OverflowException`. |
| `NumLength` | Counts magnitude digits, including negative inputs and `long.MinValue`; zero has one digit. |
| `Quadratic` | Requires a nonzero leading coefficient; zero throws `ArgumentOutOfRangeException`. |
| Real inverse hyperbolics | Undefined real-domain inputs return NaN; poles return signed infinity where appropriate. In particular, Acosch at signed zero has its corresponding infinite limit. |
| Complex `Actan` | Principal `atan(1/z)`, with `pi/2` at zero. On its imaginary cuts the real part has the sign of `Im(1/z)`. The real overload retains its `(0,pi)` convention. |

## Validation and independent references

[ArithmeticRepairTests](../tests/UMapx.Tests/ArithmeticRepairTests.cs) contains
**1,486 passing cases**, including **1,029 independent mpmath references at 100
decimal digits**, exact binary32 input rounding, real/complex extreme scales,
principal branches, polynomial residuals, and all Vieta identities.

[NumberTheoryRepairTests](../tests/UMapx.Tests/NumberTheoryRepairTests.cs) contains
**140 passing cases**, including multiple assertions/inputs per case:

- The Cartesian product of signed boundaries, both integer widths, and 300 seeded
  Int64 samples with BigInteger modular/Bezout references.
- Exhaustive elementary definitions through 1,024, pseudoprimes, and an independent
  Lucas–Lehmer certificate for the Mersenne prime `2^61-1`.
- Large semiprimes and a prime square, including
  `2147483647 * 4294967291`, with prime factors checked by independent trial division.
- Subprocess deadlines for large factorization; the existing `IsPrime(1)` deadlines
  remain enabled. Subprocess execution is not included in parent coverage totals.
- Independent digit accumulation/order checks, base-power boundaries, invalid
  inputs, and a boolean-sieve reference crossing the first segment boundary.

The [fixture generator](../tests/UMapx.Tests/Data/generate_arithmetic_reference.py)
and [embedded data](../tests/UMapx.Tests/Data/arithmetic-repair.json) are checked in.
Complex reference acceptance is `4*float.Epsilon + 3e-6*abs(expected)`; real
boundary tests use `2*float.Epsilon + 2e-6*abs(expected)` with explicit NaN/infinity
checks. Polynomial residuals are scaled by absolute polynomial terms; Vieta sums
use the magnitudes of the rounded roots to account for cancellation. No original
test budgets were relaxed.

Branch definitions follow [NIST DLMF inverse hyperbolic functions](https://dlmf.nist.gov/4.37)
and [inverse trigonometric functions](https://dlmf.nist.gov/4.23).

## Reproduction and evidence

```powershell
dotnet test tests/UMapx.Tests -c Release --no-restore -p:GeneratePackageOnBuild=false --filter "FullyQualifiedName~ArithmeticRepairTests|FullyQualifiedName~NumberTheoryRepairTests|FullyQualifiedName~CoreScalarAuditTests|FullyQualifiedName~NumberTheoryAuditTests|FullyQualifiedName~SpecialFunction"
./tools/Run-MathAudit.ps1 -NoRestore -ResultsDirectory artifacts/math-audit/arithmetic-repair/verification
```

The complete command exits with **1** while the other 311 failures remain enabled.
The recorded final run is `artifacts/math-audit/arithmetic-repair/final-confirmed`.

- [Run summary and source/input hashes](audit-arithmetic/summary.json).
- [Every-test-ID comparison with the previous 10,504-case run](audit-arithmetic/comparison.json).
- [Remaining failures](audit-arithmetic/failures.json) and [failed method families](audit-arithmetic/failures.md).
- [Source inventory and coverage](audit-arithmetic/source-inventory.md) and [uncovered methods](audit-arithmetic/uncovered-methods.json).
- [Current block ownership](audit-arithmetic/repair-blocks.json) and [remaining work order](remaining-repair-blocks-2026-09-10.md).
- [Unchanged historical block assignments](audit-special-functions/repair-blocks.json).

Coverage is **83.06% of executable lines** (28,043/33,759) and **75.95% of branches**
(10,431/13,734). All 433 production source hashes match the final audit snapshot.

## Limits and next work

Passing tests do not prove correctness over every input. The suite exercises
finite complex inputs, selected cuts/poles, and polynomial coefficient scales;
it is not an exhaustive binary32 test or a general complex-infinity contract.
The full `Sieve(int.MaxValue)` allocation was not run; segment transitions were
tested through 2,097,167 and the upper-bound arithmetic was reviewed. Pollard
factorization has bounded rho attempts and an exact trial-division fallback;
this is a correctness safeguard, not a constant-time performance guarantee.

The next recommended block is **B03: matrix indexing, complex statistics, and
distances (46 failures)**. B04 contains 67 distribution failures. The counts in
all ten remaining blocks are unchanged by these repairs.
