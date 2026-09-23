using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Regression")]
public class SpecialFunctionBoundaryTests
{
    [Theory]
    [InlineData(.1f)] [InlineData(1f)] [InlineData(30f)] [InlineData(500f)]
    public void IncompleteGammaHasCorrectEndpointsAndComplement(float shape)
    {
        Assert.Equal(0, Special.GammaP(shape, 0));
        Assert.Equal(1, Special.GammaQ(shape, 0));
        Assert.Equal(1, Special.GammaP(shape, float.PositiveInfinity));
        Assert.Equal(0, Special.GammaQ(shape, float.PositiveInfinity));
        foreach (float x in new[] { .01f, shape, shape + 1, 2 * shape })
            Close(1, (double)Special.GammaP(shape, x) + Special.GammaQ(shape, x), 1e-7, 1e-7);
    }

    [Fact]
    public void ErrorFunctionEndpointsAndInvalidInverseInputs()
    {
        Assert.Equal(1, Special.Erf(float.PositiveInfinity));
        Assert.Equal(-1, Special.Erf(float.NegativeInfinity));
        Assert.Equal(0, Special.Erfc(float.PositiveInfinity));
        Assert.Equal(2, Special.Erfc(float.NegativeInfinity));
        Assert.Equal(float.PositiveInfinity, Special.Erf(1, true));
        Assert.Equal(float.NegativeInfinity, Special.Erf(-1, true));
        Assert.True(float.IsNaN(Special.Erf(1.01f, true)));
        Assert.True(float.IsNaN(Special.Erf(-1.01f, true)));
        Assert.True(float.IsNaN(Special.Q(-.01f, true)));
        Assert.True(float.IsNaN(Special.Q(1.01f, true)));
        Assert.Equal(float.PositiveInfinity, Special.Q(0, true));
        Assert.Equal(float.NegativeInfinity, Special.Q(1, true));
    }

    [Fact]
    public void ScaledDawsonStaysFiniteAfterErfiOverflows()
    {
        Assert.Equal(float.PositiveInfinity, Special.Erfi(30));
        Assert.Equal(float.PositiveInfinity, Special.Dawson(30, false));
        Close(.016675941401059176, Special.Dawson(30, true), 1e-9, 1e-6);
        Close(0, Special.Faddeeva(30).Real, 1e-44, 0);
        Close(2 / Math.Sqrt(Math.PI) * Special.Dawson(30, true), Special.Faddeeva(30).Imag);
    }

    [Fact]
    public void IntegerIdentitiesDoNotEvaluateGammaPoles()
    {
        Assert.Equal(1, Special.FactorialUp(0, 0));
        Assert.Equal(0, Special.FactorialUp(0, 10));
        Assert.Equal(0, Special.FactorialUp(-3, 5));
        Assert.Equal(0, Special.FactorialDown(3, 5));
        Assert.Equal(0, Special.Binomial(3, 5));
        Assert.Equal(-10, Special.Binomial(-3, 3));
        Close(1, Special.Ssqrt(1));
        Assert.True(float.IsNaN(Special.Ssqrt(2, 1)));
        Assert.True(float.IsNaN(Special.Ssqrt(.1f)));
        Assert.Equal(0, Special.Harm(0));
        Assert.Equal(1836311903, Special.Fibonacci(46));
        Assert.Equal(1568397607, Special.Lucas(44));
        Assert.Throws<ArgumentOutOfRangeException>(() => Special.Fibonacci(47));
        Assert.Throws<ArgumentOutOfRangeException>(() => Special.Lucas(45));
    }

    [Theory]
    [InlineData(-1f)] [InlineData(0f)] [InlineData(.25f)] [InlineData(.5f)] [InlineData(1f)]
    public void RademacherVanishesAtItsDyadicZeros(float x)
    {
        Assert.Equal(0, Special.Rademacher(x, 2));
        Assert.Equal(Complex32.Zero, Special.Rademacher(new Complex32(x, 0), 2));
        Assert.Equal(0, Special.Rademacher(x, 200));
    }

    [Fact]
    public void HypergeometricSentinelOverloadsRetainTheirMeaning()
    {
        Close(Math.Exp(.3), Special.Hypergeom(float.NaN, float.NaN, .3f));
        Close(Math.Pow(.7, -2), Special.Hypergeom(2, float.NaN, .3f));
        Close(Math.Sinh(2) / 2, Special.Hypergeom(float.NaN, 1.5f, 1));
        Assert.True(float.IsNaN(Special.Hypergeom(1, 1, 2, 2)));
        Close(new Complex(0, -Math.PI / 2), Special.Hypergeom(new Complex32(1, 0), new Complex32(1, 0), new Complex32(2, 0), new Complex32(2, 0)));
    }

    [Fact]
    public void GammaAndZetaLimitsRemainConsistent()
    {
        Assert.Equal(float.PositiveInfinity, Special.Gamma(float.PositiveInfinity));
        Assert.Equal(float.PositiveInfinity, Special.LogGamma(float.PositiveInfinity));
        Assert.Equal(float.PositiveInfinity, Special.DiGamma(float.PositiveInfinity));
        Assert.Equal(0, Special.TriGamma(float.PositiveInfinity));
        Assert.Equal(float.PositiveInfinity, Special.Zeta(1));
        Assert.Equal(1, Special.Zeta(float.PositiveInfinity));
        Assert.Equal(-.5f, Special.Zeta(0));
        Assert.Equal(0, Special.Zeta(-2));
    }

    [Theory]
    [InlineData(1e-20f)] [InlineData(1e-38f)]
    public void GammaUpperTailSurvivesVanishingShape(float shape)
    {
        // Q(a,x)/a tends to E1(x) as a tends to zero.
        Close(.5597735947761608, Special.GammaQ(shape, .5f) / (double)shape, 1e-6, 1e-6);
        Close(.5597735947761608, Special.GammaQ(new Complex32(shape, 0), new Complex32(.5f, 0)).Real / (double)shape, 1e-6, 1e-6);
    }

    [Theory]
    [InlineData(0)] [InlineData(1)] [InlineData(5)] [InlineData(20)]
    public void PolynomialEndpointsAndNegativeOrders(int n)
    {
        Close(1, Special.ChebyshevT(1f, n));
        Close(n + 1, Special.ChebyshevU(1f, n));
        Close(Special.ChebyshevT(.3f, n), Special.ChebyshevT(.3f, -n));
        Close(-Special.ChebyshevU(.3f, n), Special.ChebyshevU(.3f, -n - 2));
        Close(Special.Euler(n) / Math.Pow(2, n), Special.Euler(n, .5f));
    }

    [Fact]
    public void BesselWronskianHoldsAtComplexArgument()
    {
        Complex32 z = new(10, 10);
        int n = 20;
        Complex left = (Complex)Special.J(z, n + 1) * (Complex)Special.Y(z, n)
            - (Complex)Special.J(z, n) * (Complex)Special.Y(z, n + 1);
        Close(2 / (Math.PI * (Complex)z), left, 1e-7, 1e-5);
    }

    [Fact]
    public void MinkowskiFixesIntegerArguments()
    {
        foreach (long n in new[] { long.MinValue, -100, -1, 0, 1, 100, long.MaxValue }) Assert.Equal((float)n, Special.Minkowski(n));
    }
}
