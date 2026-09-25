using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class NumberTheoryTests
{
    [Theory]
    [InlineData(12, 18)]
    [InlineData(-12, 18)]
    [InlineData(12, -18)]
    [InlineData(18, -6)]
    [InlineData(12345, 54321)]
    [InlineData(50000, 50000)]
    public void IntegerGcdLcmAndBezoutAgreeWithExactArithmetic(int a, int b)
    {
        long gcd = (long)BigInteger.GreatestCommonDivisor(a, b);
        long lcm = Math.Abs((long)a / gcd * b);
        Assert.Equal(gcd, Maths.Gcd(a, b));
        Assert.Equal(gcd, Maths.Gcd((long)a, b));
        Assert.Equal(lcm, Maths.Lcm(a, b));
        Assert.Equal(lcm, Maths.Lcm((long)a, b));
        var v = Maths.Euclidean(a, b);
        Assert.Equal(gcd, v[0]);
        Assert.Equal(gcd, (long)a * v[1] + (long)b * v[2]);
        var w = Maths.Euclidean((long)a, b);
        Assert.Equal(gcd, w[0]);
        Assert.Equal(gcd, a * w[1] + b * w[2]);
    }

    [Theory]
    [InlineData(2, 0, 7)]
    [InlineData(2, 10, 17)]
    [InlineData(17, 13, 41)]
    [InlineData(123456789, 13, 1000000007)]
    public void ModularExponentiationMatchesBigInteger(int a, int exponent, int modulus)
    {
        long expected = (long)BigInteger.ModPow(a, exponent, modulus);
        foreach (bool modified in new[]
        {
            false,
            true
        }

        )
        {
            Assert.Equal(expected, Maths.ModPow(a, exponent, modulus, modified));
            Assert.Equal(expected, Maths.ModPow((long)a, exponent, modulus, modified));
        }
    }

    [Fact]
    public void ModularArithmeticAndCoprimeSearchSatisfyTheirDefinitions()
    {
        foreach (int a in new[]
        {
            -17,
            -3,
            1,
            5,
            17
        }

        )
            foreach (int n in new[]
            {
                3,
                7,
                11
            }

            )
            {
                int expected = (a % n + n) % n;
                Assert.Equal(expected, Maths.Mod(a, n));
                Assert.Equal(expected, Maths.Mod((long)a, n));
                NumericAssert.Close(expected, Maths.Mod((float)a, n));
                if (BigInteger.GreatestCommonDivisor(a, n) == 1)
                {
                    Assert.Equal(1, Maths.Mod(a * Maths.ModInv(a, n), n));
                    Assert.Equal(1, Maths.Mod((long)a * Maths.ModInv((long)a, n), n));
                }
            }
    }

    [Theory]
    [InlineData(6)]
    [InlineData(12)]
    [InlineData(35)]
    public void CoprimeSearchActuallyReturnsACoprime(int n)
    {
        Assert.Equal(BigInteger.One, BigInteger.GreatestCommonDivisor(n, Maths.Coprime(n, 2)));
        Assert.Equal(BigInteger.One, BigInteger.GreatestCommonDivisor(n, Maths.Coprime((long)n, 2)));
    }

    [Fact]
    public void PrimeSieveFactorizationTotientAndRadicalAgreeWithIntegerDefinitions()
    {
        bool Prime(int n) => n >= 2 && !Enumerable.Range(2, Math.Max(0, (int)Math.Sqrt(n) - 1)).Any(d => n % d == 0);
        var primes = Enumerable.Range(0, 101).Where(Prime).ToArray();
        Assert.Equal(primes, Maths.Sieve(100));
        foreach (int n in Enumerable.Range(2, 99))
        {
            Assert.Equal(Prime(n), Maths.IsPrime(n));
            Assert.Equal(Prime(n), Maths.IsPrime((long)n));
            int totient = Enumerable.Range(1, n).Count(k => BigInteger.GreatestCommonDivisor(k, n) == 1);
            Assert.Equal(totient, Maths.Etf(n));
            Assert.Equal(totient, Maths.Etf((long)n));
            int radical = primes.Where(p => n % p == 0).Aggregate(1, (a, b) => a * b);
            Assert.Equal(radical, Maths.Radical(n));
            Assert.Equal(radical, Maths.Radical((long)n));
        }

        foreach (int n in new[]
        {
            12,
            30,
            72,
            121
        }

        )
        {
            var factors = Maths.Itf(n);
            Assert.Equal(n, factors.Aggregate(1, (a, b) => a * b));
            var longs = Maths.Itf((long)n);
            Assert.Equal(n, longs.Aggregate(1L, (a, b) => a * b));
            Assert.Equal(primes.Where(p => n % p == 0), Maths.Itf(n, true).OrderBy(p => p));
            Assert.Equal(primes.Where(p => n % p == 0).Select(p => (long)p), Maths.Itf((long)n, true).OrderBy(p => p));
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void OneIsNotPrime(bool wide)
    {
        Assert.False(wide ? Maths.IsPrime(1L) : Maths.IsPrime(1));
    }

    [Theory]
    [InlineData(0L)]
    [InlineData(1L)]
    [InlineData(123456789L)]
    [InlineData(987654321012345L)]
    public void BaseConversionsPreserveExactIntegerValues(long value)
    {
        foreach (int radix in new[]
        {
            2,
            3,
            8,
            10,
            16,
            36
        }

        )
        {
            var digits = Maths.Decimal2Base(value, radix);
            Assert.All(digits, d => Assert.InRange(d, 0, radix - 1));
            Assert.Equal(value, Maths.Base2Decimal(digits, radix));
        }

        Assert.Equal(value, Maths.Vector2Numeral(Maths.Numeral2Vector(value)));
    }

    // Independent cases ensure one defect cannot conceal another operation or overload.
    public static IEnumerable<object[]> IndependentCases()
    {
        foreach (string operation in new[]
        {
            "Gcd",
            "Lcm",
            "Euclidean",
            "ModPow",
            "Coprime",
            "IsPrime",
            "Etf",
            "Radical",
            "Itf",
            "Digits"
        }

        )
            foreach (int n in new[]
            {
                2,
                3,
                6,
                12,
                35,
                72,
                97,
                100,
                121,
                50000
            }

            )
                foreach (bool wide in new[]
                {
                    false,
                    true
                }

                )
                    yield return new object[]
                    {
                        operation,
                        n,
                        wide
                    };
    }

    [Theory]
    [MemberData(nameof(IndependentCases))]
    public void IntegerOperationsIndependentlyMatchExactArithmetic(string operation, int n, bool wide)
    {
        bool Prime(int p) => p >= 2 && !Enumerable.Range(2, Math.Max(0, (int)Math.Sqrt(p) - 1)).Any(d => p % d == 0);
        var primes = Enumerable.Range(2, n - 1).Where(p => n % p == 0 && Prime(p)).Select(p => (long)p).ToArray();
        switch (operation)
        {
            case "Gcd":
                Assert.Equal(n, wide ? Maths.Gcd(-(long)n, n) : Maths.Gcd(-n, n));
                break;
            case "Lcm":
                Assert.Equal(n, wide ? Maths.Lcm((long)n, n) : Maths.Lcm(n, n));
                break;
            case "Euclidean":
                var v = wide ? Maths.Euclidean(-(long)n, n + 1) : Maths.Euclidean(-n, n + 1).Select(x => (long)x).ToArray();
                Assert.Equal(1, v[0]);
                Assert.Equal(1, -n * v[1] + (n + 1) * v[2]);
                break;
            case "ModPow":
                Assert.Equal(1, wide ? Maths.ModPow((long)n, 0, 101) : Maths.ModPow(n, 0, 101));
                break;
            case "Coprime":
                Assert.Equal(BigInteger.One, BigInteger.GreatestCommonDivisor(n, wide ? Maths.Coprime((long)n, 2) : Maths.Coprime(n, 2)));
                break;
            case "IsPrime":
                Assert.Equal(Prime(n), wide ? Maths.IsPrime((long)n) : Maths.IsPrime(n));
                break;
            case "Etf":
                Assert.Equal(Enumerable.Range(1, n).Count(k => BigInteger.GreatestCommonDivisor(k, n) == 1), wide ? Maths.Etf((long)n) : Maths.Etf(n));
                break;
            case "Radical":
                Assert.Equal(primes.Aggregate(1L, (a, b) => a * b), wide ? Maths.Radical((long)n) : Maths.Radical(n));
                break;
            case "Itf":
                var f = wide ? Maths.Itf((long)n) : Maths.Itf(n).Select(p => (long)p).ToArray();
                Assert.Equal(n, f.Aggregate(1L, (a, b) => a * b));
                Assert.Equal(primes, (wide ? Maths.Itf((long)n, true) : Maths.Itf(n, true).Select(p => (long)p)).OrderBy(p => p));
                break;
            case "Digits":
                Assert.Equal(n, wide ? Maths.Vector2Numeral(Maths.Numeral2Vector((long)n)) : Maths.Vector2Numeral(Maths.Numeral2Vector(n)));
                break;
        }
    }

    [Theory]
    [InlineData(0L)]
    [InlineData(1L)]
    [InlineData(123456789L)]
    [InlineData(987654321012345L)]
    public void DecimalDigitVectorsRoundTripIndependentlyOfBaseConversion(long value) => Assert.Equal(value, Maths.Vector2Numeral(Maths.Numeral2Vector(value)));
    [Theory]
    [InlineData(25)]
    [InlineData(49)]
    [InlineData(169)]
    [InlineData(341)]
    [InlineData(561)]
    public void CompositeNumbersAreNotDeclaredPrimeWhenPollardFailsToSplitThem(int n)
    {
        Assert.False(Maths.IsPrime(n));
        Assert.False(Maths.IsPrime((long)n));
    }

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(49)]
    [InlineData(2097152)]
    [InlineData(2097153)]
    [InlineData(2097167)]
    public void SegmentedSieveAgreesWithAnIndependentBooleanSieve(int limit)
    {
        var composite = new bool[limit + 1];
        for (int divisor = 2; divisor <= limit / divisor; divisor++)
            if (!composite[divisor])
                for (int multiple = divisor * divisor; multiple <= limit; multiple += divisor)
                    composite[multiple] = true;
        var expected = Enumerable.Range(0, limit + 1).Where(n => n >= 2 && !composite[n]);
        Assert.Equal(expected, Maths.Sieve(limit));
    }

    public static IEnumerable<object[]> SignedPairs()
    {
        long[] values =
        {
            long.MinValue,
            long.MinValue + 1,
            int.MinValue,
            -50000,
            -1,
            0,
            1,
            50000,
            int.MaxValue,
            long.MaxValue
        };
        foreach (long a in values)
            foreach (long b in values)
                yield return new object[]
                {
                    a,
                    b
                };
    }

    [Theory, MemberData(nameof(SignedPairs))]
    public void SignedBoundariesPreserveGcdLcmBezoutAndRemainders(long a, long b)
    {
        BigInteger gcd = BigInteger.GreatestCommonDivisor(a, b);
        BigInteger lcm = gcd == 0 ? 0 : BigInteger.Abs((BigInteger)a / gcd * b);
        if (gcd <= long.MaxValue)
        {
            Assert.Equal((long)gcd, Maths.Gcd(a, b));
            var bezout = Maths.Euclidean(a, b);
            Assert.Equal((long)gcd, bezout[0]);
            Assert.Equal(gcd, (BigInteger)a * bezout[1] + (BigInteger)b * bezout[2]);
        }
        else
        {
            Assert.Throws<OverflowException>(() => Maths.Gcd(a, b));
            Assert.Throws<OverflowException>(() => Maths.Euclidean(a, b));
        }

        if (lcm <= long.MaxValue)
            Assert.Equal((long)lcm, Maths.Lcm(a, b));
        else
            Assert.Throws<OverflowException>(() => Maths.Lcm(a, b));
        if (b != 0)
        {
            BigInteger modulus = BigInteger.Abs(b);
            Assert.Equal((long)(((BigInteger)a % modulus + modulus) % modulus), Maths.Mod(a, b));
            long inverse = Maths.ModInv(a, b);
            if (gcd == 1)
                Assert.Equal(BigInteger.One % modulus, (((BigInteger)a * inverse) % modulus + modulus) % modulus);
            else
                Assert.Equal(0, inverse);
        }

        if (a < int.MinValue || a > int.MaxValue || b < int.MinValue || b > int.MaxValue)
            return;
        int x = (int)a, y = (int)b;
        if (gcd <= int.MaxValue)
        {
            Assert.Equal((int)gcd, Maths.Gcd(x, y));
            var bezout = Maths.Euclidean(x, y);
            Assert.Equal(gcd, (BigInteger)x * bezout[1] + (BigInteger)y * bezout[2]);
        }
        else
        {
            Assert.Throws<OverflowException>(() => Maths.Gcd(x, y));
            Assert.Throws<OverflowException>(() => Maths.Euclidean(x, y));
        }

        if (lcm <= int.MaxValue)
            Assert.Equal((int)lcm, Maths.Lcm(x, y));
        else
            Assert.Throws<OverflowException>(() => Maths.Lcm(x, y));
        if (y != 0)
        {
            Assert.Equal(Maths.Mod(a, b), Maths.Mod(x, y));
            Assert.Equal(Maths.ModInv(a, b), Maths.ModInv(x, y));
        }
    }

    [Theory]
    [InlineData(17)]
    [InlineData(731)]
    [InlineData(2026)]
    public void ModularPowersAndBezoutMatchBigIntegerAcrossRandomInt64Values(int seed)
    {
        var random = new Random(seed);
        for (int i = 0; i < 100; i++)
        {
            long a = random.NextInt64(long.MinValue, long.MaxValue);
            long b = random.NextInt64(long.MinValue, long.MaxValue);
            if (b == 0)
                b = 1;
            var v = Maths.Euclidean(a, b);
            Assert.Equal(BigInteger.GreatestCommonDivisor(a, b), (BigInteger)a * v[1] + (BigInteger)b * v[2]);
            foreach (long exponent in new[]
            {
                0,
                1,
                random.NextInt64(2, long.MaxValue)
            }

            )
            {
                BigInteger modulus = BigInteger.Abs(b);
                BigInteger expected = BigInteger.ModPow(((BigInteger)a % modulus + modulus) % modulus, exponent, modulus);
                foreach (bool modified in new[]
                {
                    false,
                    true
                }

                )
                    Assert.Equal((long)expected, Maths.ModPow(a, exponent, b, modified));
            }
        }

        foreach (bool modified in new[]
        {
            false,
            true
        }

        )
        {
            Assert.Equal(0, Maths.ModPow(long.MinValue, 1, long.MinValue, modified));
            Assert.Equal(1, Maths.ModPow(long.MinValue, 0, long.MinValue, modified));
            Assert.Equal(long.MaxValue, Maths.ModPow(-1, long.MaxValue, long.MinValue, modified));
            Assert.Equal(0, Maths.ModPow(0, 0, -1, modified));
            Assert.Equal(0, Maths.ModPow(0L, 0, 1, modified));
        }
    }

    [Fact]
    public void SmallIntegersMatchExhaustiveDefinitions()
    {
        for (int n = -32; n <= 1024; n++)
        {
            bool prime = TrialPrime(n);
            Assert.Equal(prime, Maths.IsPrime(n));
            Assert.Equal(prime, Maths.IsPrime((long)n));
            if (n < 1)
                continue;
            int[] divisors = Enumerable.Range(2, n - 1).Where(p => n % p == 0 && TrialPrime(p)).ToArray();
            Assert.Equal(divisors, Maths.Itf(n, true));
            Assert.Equal(divisors.Select(p => (long)p), Maths.Itf((long)n, true));
            Assert.Equal(n, Maths.Itf(n).Aggregate(1, (a, b) => a * b));
            Assert.Equal(n, Maths.Itf((long)n).Aggregate(1L, (a, b) => a * b));
            int totient = Enumerable.Range(1, n).Count(k => BigInteger.GreatestCommonDivisor(k, n) == 1);
            Assert.Equal(totient, Maths.Etf(n));
            Assert.Equal(totient, Maths.Etf((long)n));
            int radical = divisors.Aggregate(1, (a, b) => a * b);
            Assert.Equal(radical, Maths.Radical(n));
            Assert.Equal(radical, Maths.Radical((long)n));
            int factor = Maths.Pollard(n);
            Assert.Equal(0, n % factor);
            if (!prime && n > 1)
                Assert.InRange(factor, 2, n - 1);
        }
    }

    [Theory]
    [InlineData(25L)]
    [InlineData(169L)]
    [InlineData(561L)]
    [InlineData(3215031751L)]
    [InlineData(341550071728321L)]
    [InlineData(3825123056546413051L)]
    [InlineData(long.MaxValue)]
    public void StrongPseudoprimesAndPrimeSquaresAreComposite(long n) => Assert.False(Maths.IsPrime(n));
    [Theory]
    [InlineData(long.MinValue)]
    [InlineData(int.MinValue)]
    [InlineData(-1L)]
    [InlineData(0L)]
    [InlineData(1L)]
    public void NonpositiveValuesAndOneAreNotPrime(long n)
    {
        Assert.False(Maths.IsPrime(n));
        if (n >= int.MinValue)
            Assert.False(Maths.IsPrime((int)n));
    }

    [Fact]
    public void LargeMersennePrimeHasAnIndependentLucasLehmerCertificate()
    {
        long prime = (1L << 61) - 1;
        BigInteger residue = 4;
        for (int i = 0; i < 59; i++)
            residue = (residue * residue - 2) % prime;
        Assert.Equal(BigInteger.Zero, residue);
        Assert.True(Maths.IsPrime(prime));
        Assert.Equal(prime, Maths.Pollard(prime));
        Assert.Equal(new[] { prime }, Maths.Itf(prime));
        Assert.Equal(prime - 1, Maths.Etf(prime));
    }

    [Theory]
    [InlineData(1000000007L, 1000000009L)]
    [InlineData(2147483647L, 2147483647L)]
    [InlineData(2147483647L, 4294967291L)]
    public void LargeSemiprimesSplitIntoTheirPrimeFactors(long p, long q)
    {
        Assert.True(TrialPrime(p));
        Assert.True(TrialPrime(q));
        long n = checked(p * q);
        Assert.Equal(new[] { p, q }, Maths.Itf(n));
        long[] distinct = new[]
        {
            p,
            q
        }.Distinct().ToArray();
        Assert.Equal(distinct, Maths.Itf(n, true));
        Assert.Equal(distinct.Aggregate(1L, (a, b) => a * b), Maths.Radical(n));
        Assert.Equal(p == q ? p * (p - 1) : (p - 1) * (q - 1), Maths.Etf(n));
    }

    [Theory]
    [InlineData(0L)]
    [InlineData(1L)]
    [InlineData(9L)]
    [InlineData(10L)]
    [InlineData(99L)]
    [InlineData(100L)]
    [InlineData(123456789L)]
    [InlineData(987654321012345L)]
    [InlineData(long.MaxValue)]
    public void DigitArraysUseExactIntegerArithmeticAndDocumentedOrdering(long n)
    {
        int[] expectedDecimal = n.ToString().Select(c => c - '0').ToArray();
        Assert.Equal(expectedDecimal, Maths.Numeral2Vector(n));
        Assert.Equal(n, Maths.Vector2Numeral(expectedDecimal));
        foreach (int radix in new[]
        {
            2,
            3,
            8,
            10,
            16,
            36,
            int.MaxValue
        }

        )
        {
            var digits = Maths.Decimal2Base(n, radix);
            BigInteger expected = 0, weight = 1;
            foreach (int digit in digits)
            {
                expected += digit * weight;
                weight *= radix;
            }

            Assert.Equal((BigInteger)n, expected);
            Assert.Equal(n, Maths.Base2Decimal(digits, radix));
            Assert.Equal(digits.Length, Maths.NumLength(n, radix));
            Assert.Equal(digits.Length, Maths.NumLength(-n, radix));
        }
    }

    [Fact]
    public void PowersOfEachBaseHaveExactDigitLengths()
    {
        foreach (int radix in new[]
        {
            2,
            3,
            10,
            16,
            36
        }

        )
        {
            for (long power = radix, length = 2; power <= long.MaxValue / radix; power *= radix, length++)
            {
                Assert.Equal(length - 1, Maths.NumLength(power - 1, radix));
                Assert.Equal(length, Maths.NumLength(power, radix));
                Assert.Equal(length, Maths.NumLength(power + 1, radix));
            }
        }

        Assert.Equal(64, Maths.NumLength(long.MinValue, 2));
        Assert.Equal(19, Maths.NumLength(long.MinValue, 10));
        Assert.Equal(new[] { 0, 1, 0, 1 }, Maths.Decimal2Base(10, 2));
    }

    [Fact]
    public void CoprimeSearchStartsAtTheInclusiveBound()
    {
        for (int n = -20; n <= 20; n++)
            for (int start = -8; start <= 8; start++)
            {
                if (n == 0 && start > 1)
                    continue;
                int expected = Enumerable.Range(start, 50).First(p => BigInteger.GreatestCommonDivisor(n, p) == 1);
                Assert.Equal(expected, Maths.Coprime(n, start));
                Assert.Equal(expected, Maths.Coprime((long)n, start));
            }

        Assert.Equal(long.MinValue + 1, Maths.Coprime(long.MinValue, long.MinValue));
        Assert.Throws<OverflowException>(() => Maths.Coprime(0, 2));
        Assert.Throws<OverflowException>(() => Maths.Coprime(0L, 2));
        Assert.Throws<OverflowException>(() => Maths.Coprime(int.MaxValue, int.MaxValue));
        Assert.Throws<OverflowException>(() => Maths.Coprime(long.MaxValue, long.MaxValue));
    }

    [Fact]
    public void UndefinedInputsAndUnrepresentableResultsFailExplicitly()
    {
        Assert.Empty(Maths.Itf(1));
        Assert.Empty(Maths.Itf(1L, true));
        Assert.Equal(1, Maths.Pollard(1));
        Assert.Equal(1, Maths.Pollard(1L));
        foreach (int n in new[]
        {
            0,
            -1,
            int.MinValue
        }

        )
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Itf(n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Itf((long)n, true));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Pollard(n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Pollard((long)n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Etf(n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Etf((long)n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Radical(n));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Radical((long)n));
        }

        Assert.Throws<DivideByZeroException>(() => Maths.Mod(1, 0));
        Assert.Throws<DivideByZeroException>(() => Maths.Mod(1L, 0));
        Assert.Throws<DivideByZeroException>(() => Maths.ModInv(1, 0));
        Assert.Throws<DivideByZeroException>(() => Maths.ModInv(1L, 0));
        foreach (bool modified in new[]
        {
            false,
            true
        }

        )
        {
            Assert.Throws<DivideByZeroException>(() => Maths.ModPow(2, 3, 0, modified));
            Assert.Throws<DivideByZeroException>(() => Maths.ModPow(2L, 3, 0, modified));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.ModPow(2, -1, 3, modified));
            Assert.Throws<ArgumentOutOfRangeException>(() => Maths.ModPow(2L, -1, 3, modified));
        }

        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Decimal2Base(-1, 10));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Numeral2Vector(long.MinValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.NumLength(1, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Decimal2Base(1, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Base2Decimal(new[] { 1 }, -2));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Base2Decimal(new[] { 2 }, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Vector2Numeral(new[] { -1 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Maths.Vector2Numeral(new[] { 10 }));
        Assert.Throws<ArgumentNullException>(() => Maths.Base2Decimal(null!, 10));
        Assert.Throws<ArgumentNullException>(() => Maths.Vector2Numeral(null!));
        Assert.Throws<OverflowException>(() => Maths.Base2Decimal(Enumerable.Repeat(1, 64).ToArray(), 2));
        Assert.Throws<OverflowException>(() => Maths.Vector2Numeral("9223372036854775808".Select(c => c - '0').ToArray()));
        Assert.Equal(0, Maths.Base2Decimal(Array.Empty<int>(), 10));
        Assert.Equal(0, Maths.Vector2Numeral(Array.Empty<int>()));
    }

    private static bool TrialPrime(long n)
    {
        if (n < 2)
            return false;
        for (long divisor = 2; divisor <= n / divisor; divisor++)
            if (n % divisor == 0)
                return false;
        return true;
    }
}
