using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Core")]
public class NumberTheoryAuditTests
{
    [Theory] [InlineData(12,18)] [InlineData(-12,18)] [InlineData(12,-18)] [InlineData(18,-6)] [InlineData(12345,54321)] [InlineData(50000,50000)]
    public void IntegerGcdLcmAndBezoutAgreeWithExactArithmetic(int a,int b)
    {
        long gcd=(long)BigInteger.GreatestCommonDivisor(a,b);
        long lcm=Math.Abs((long)a/gcd*b);
        Assert.Equal(gcd,Maths.Gcd(a,b)); Assert.Equal(gcd,Maths.Gcd((long)a,b));
        Assert.Equal(lcm,Maths.Lcm(a,b)); Assert.Equal(lcm,Maths.Lcm((long)a,b));
        var v=Maths.Euclidean(a,b); Assert.Equal(gcd,v[0]); Assert.Equal(gcd,(long)a*v[1]+(long)b*v[2]);
        var w=Maths.Euclidean((long)a,b); Assert.Equal(gcd,w[0]); Assert.Equal(gcd,a*w[1]+b*w[2]);
    }

    [Theory] [InlineData(2,0,7)] [InlineData(2,10,17)] [InlineData(17,13,41)] [InlineData(123456789,13,1000000007)]
    public void ModularExponentiationMatchesBigInteger(int a,int exponent,int modulus)
    {
        long expected=(long)BigInteger.ModPow(a,exponent,modulus);
        foreach(bool modified in new[]{false,true})
        {
            Assert.Equal(expected,Maths.ModPow(a,exponent,modulus,modified));
            Assert.Equal(expected,Maths.ModPow((long)a,exponent,modulus,modified));
        }
    }

    [Fact]
    public void ModularArithmeticAndCoprimeSearchSatisfyTheirDefinitions()
    {
        foreach(int a in new[]{-17,-3,1,5,17})foreach(int n in new[]{3,7,11})
        {
            int expected=(a%n+n)%n;
            Assert.Equal(expected,Maths.Mod(a,n)); Assert.Equal(expected,Maths.Mod((long)a,n));
            NumericAssert.Close(expected,Maths.Mod((float)a,n));
            if(BigInteger.GreatestCommonDivisor(a,n)==1)
            {
                Assert.Equal(1,Maths.Mod(a*Maths.ModInv(a,n),n));
                Assert.Equal(1,Maths.Mod((long)a*Maths.ModInv((long)a,n),n));
            }
        }
    }

    [Theory] [InlineData(6)] [InlineData(12)] [InlineData(35)]
    public void CoprimeSearchActuallyReturnsACoprime(int n)
    {
        Assert.Equal(BigInteger.One,BigInteger.GreatestCommonDivisor(n,Maths.Coprime(n,2)));
        Assert.Equal(BigInteger.One,BigInteger.GreatestCommonDivisor(n,Maths.Coprime((long)n,2)));
    }

    [Fact]
    public void PrimeSieveFactorizationTotientAndRadicalAgreeWithIntegerDefinitions()
    {
        bool Prime(int n)=>n>=2&&!Enumerable.Range(2,Math.Max(0,(int)Math.Sqrt(n)-1)).Any(d=>n%d==0);
        var primes=Enumerable.Range(0,101).Where(Prime).ToArray();
        Assert.Equal(primes,Maths.Sieve(100));
        foreach(int n in Enumerable.Range(2,99))
        {
            Assert.Equal(Prime(n),Maths.IsPrime(n)); Assert.Equal(Prime(n),Maths.IsPrime((long)n));
            int totient=Enumerable.Range(1,n).Count(k=>BigInteger.GreatestCommonDivisor(k,n)==1);
            Assert.Equal(totient,Maths.Etf(n)); Assert.Equal(totient,Maths.Etf((long)n));
            int radical=primes.Where(p=>n%p==0).Aggregate(1,(a,b)=>a*b);
            Assert.Equal(radical,Maths.Radical(n)); Assert.Equal(radical,Maths.Radical((long)n));
        }
        foreach(int n in new[]{12,30,72,121})
        {
            var factors=Maths.Itf(n); Assert.Equal(n,factors.Aggregate(1,(a,b)=>a*b));
            var longs=Maths.Itf((long)n); Assert.Equal(n,longs.Aggregate(1L,(a,b)=>a*b));
            Assert.Equal(primes.Where(p=>n%p==0),Maths.Itf(n,true).OrderBy(p=>p));
            Assert.Equal(primes.Where(p=>n%p==0).Select(p=>(long)p),Maths.Itf((long)n,true).OrderBy(p=>p));
        }
    }

    [Theory] [InlineData("IsPrimeInt")] [InlineData("IsPrimeLong")]
    public async Task OneIsNotPrimeAndTheCheckTerminates(string operation)
    {
        Assert.Equal("False",await AuditProcess.RunAsync(operation,"1"));
    }

    [Theory] [InlineData(0L)] [InlineData(1L)] [InlineData(123456789L)] [InlineData(987654321012345L)]
    public void BaseConversionsPreserveExactIntegerValues(long value)
    {
        foreach(int radix in new[]{2,3,8,10,16,36})
        {
            var digits=Maths.Decimal2Base(value,radix);
            Assert.All(digits,d=>Assert.InRange(d,0,radix-1));
            Assert.Equal(value,Maths.Base2Decimal(digits,radix));
        }
        Assert.Equal(value,Maths.Vector2Numeral(Maths.Numeral2Vector(value)));
    }

    // Independent cases ensure one defect cannot conceal another operation or overload.
    public static IEnumerable<object[]> IndependentCases()
    {
        foreach(string operation in new[]{"Gcd","Lcm","Euclidean","ModPow","Coprime","IsPrime","Etf","Radical","Itf","Digits"})
        foreach(int n in new[]{2,3,6,12,35,72,97,100,121,50000})
        foreach(bool wide in new[]{false,true})yield return new object[]{operation,n,wide};
    }
    [Theory] [MemberData(nameof(IndependentCases))]
    public void IntegerOperationsIndependentlyMatchExactArithmetic(string operation,int n,bool wide)
    {
        bool Prime(int p)=>p>=2&&!Enumerable.Range(2,Math.Max(0,(int)Math.Sqrt(p)-1)).Any(d=>p%d==0);
        var primes=Enumerable.Range(2,n-1).Where(p=>n%p==0&&Prime(p)).Select(p=>(long)p).ToArray();
        switch(operation)
        {
            case "Gcd":Assert.Equal(n,wide?Maths.Gcd(-(long)n,n):Maths.Gcd(-n,n));break;
            case "Lcm":Assert.Equal(n,wide?Maths.Lcm((long)n,n):Maths.Lcm(n,n));break;
            case "Euclidean":
                var v=wide?Maths.Euclidean(-(long)n,n+1):Maths.Euclidean(-n,n+1).Select(x=>(long)x).ToArray();
                Assert.Equal(1,v[0]);Assert.Equal(1,-n*v[1]+(n+1)*v[2]);break;
            case "ModPow":Assert.Equal(1,wide?Maths.ModPow((long)n,0,101):Maths.ModPow(n,0,101));break;
            case "Coprime":Assert.Equal(BigInteger.One,BigInteger.GreatestCommonDivisor(n,wide?Maths.Coprime((long)n,2):Maths.Coprime(n,2)));break;
            case "IsPrime":Assert.Equal(Prime(n),wide?Maths.IsPrime((long)n):Maths.IsPrime(n));break;
            case "Etf":Assert.Equal(Enumerable.Range(1,n).Count(k=>BigInteger.GreatestCommonDivisor(k,n)==1),wide?Maths.Etf((long)n):Maths.Etf(n));break;
            case "Radical":Assert.Equal(primes.Aggregate(1L,(a,b)=>a*b),wide?Maths.Radical((long)n):Maths.Radical(n));break;
            case "Itf":
                var f=wide?Maths.Itf((long)n):Maths.Itf(n).Select(p=>(long)p).ToArray();
                Assert.Equal(n,f.Aggregate(1L,(a,b)=>a*b));
                Assert.Equal(primes,(wide?Maths.Itf((long)n,true):Maths.Itf(n,true).Select(p=>(long)p)).OrderBy(p=>p));break;
            case "Digits":Assert.Equal(n,wide?Maths.Vector2Numeral(Maths.Numeral2Vector((long)n)):Maths.Vector2Numeral(Maths.Numeral2Vector(n)));break;
        }
    }

    [Theory] [InlineData(0L)] [InlineData(1L)] [InlineData(123456789L)] [InlineData(987654321012345L)]
    public void DecimalDigitVectorsRoundTripIndependentlyOfBaseConversion(long value)=>Assert.Equal(value,Maths.Vector2Numeral(Maths.Numeral2Vector(value)));

    [Theory] [InlineData(25)] [InlineData(49)] [InlineData(169)] [InlineData(341)] [InlineData(561)]
    public void CompositeNumbersAreNotDeclaredPrimeWhenPollardFailsToSplitThem(int n)
    {
        Assert.False(Maths.IsPrime(n));Assert.False(Maths.IsPrime((long)n));
    }
}
