using System.Numerics;
using UMapx.Analysis;
using UMapx.Colorspace;
using UMapx.Core;
using UMapx.Distribution;
using UMapx.Response;
using UMapx.Wavelet;
using UMapx.Window;
using Xunit;
using static UMapx.Tests.NumericAssert;
using Poisson = UMapx.Distribution.Poisson;

namespace UMapx.Tests;

// These tests assert the correct mathematics. Failures are intentional evidence
// of existing library defects, not assertions of the current incorrect behavior.
[Trait("Category", "Regression")]
public class MathematicalRegressionTests
{
    [Fact] public void GammaSeriesMustIncludeGammaNormalization() =>
        Close(0.5595067149347877, Special.GammaP(5, 5));

    [Fact] public void GammaContinuedFractionMustUseNextDenominator() =>
        Close(4 * Math.Exp(-3), Special.GammaQ(2, 3));

    [Theory]
    [InlineData(-1f)] [InlineData(0f)] [InlineData(1f)]
    public void QuadraticHandlesZeroLinearCoefficient(float constant)
    {
        var roots = Maths.Quadratic(1, 0, constant);
        Assert.Equal(2, roots.Length);
        foreach (var root in roots) Close(0, (Complex)root * (Complex)root + constant);
        Close(0, (Complex)roots[0] + (Complex)roots[1]);
        Close(constant, (Complex)roots[0] * (Complex)roots[1]);
    }

    [Theory]
    [InlineData(0f, 0f, -1f)] [InlineData(0f, -3f, -2f)]
    public void CubicUsesRealCubeRootsForNegativeRadicands(float a, float b, float c)
    {
        var roots = Maths.Cubic(a, b, c);
        Assert.Equal(3, roots.Length);
        foreach (Complex root in roots) Close(0, ((root + a) * root + b) * root + c);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void BiorthogonalWaveletReconstructsImpulse(bool normalized)
    {
        var x = new float[8]; x[0] = 1;
        var wavelet = new WaveletDecomposition(WaveletPacket.Bior13, 1, normalized);
        Close(x, wavelet.Backward(wavelet.Forward(x)));
    }

    [Fact] public void ComplexArccotangentAgreesWithPositiveRealBranch() =>
        Close(new Complex(Math.PI / 4, 0), Maths.Actan(new Complex32(1, 0)));

    [Fact] public void ComplexAcoshUsesPrincipalBranchOnNegativeRealAxis() =>
        Close(new Complex(Math.Acosh(2), Math.PI), Maths.Acosh(new Complex32(-2, 0)));

    [Fact] public void ComplexErfAgreesWithEntireFunctionReference() =>
        Close(new Complex(-1.0035022433130363, .004740903031294336), Special.Erf(new Complex32(-2, .5f)));

    [Fact] public void LambertWReturnsRequestedPrincipalBranch() =>
        Close(new Complex(.21453028412564626, .3510955496497862), Special.LambertW(new Complex32(.1f, .5f), 0));

    [Fact] public void BesselYIncludesBothExponentialTerms() => Close(.08825696421567696, Special.Y(1f, 0));

    [Fact] public void BesselJAsymptoticAccuracyDependsOnOrder() => Close(.14853180559607407, Special.J(21f, 10));

    [Fact] public void BesselJAtImaginaryArgumentSatisfiesConnectionToI() =>
        Close(new Complex(-145831809975.96713, 0), Special.J(new Complex32(0, 30), 10));

    [Fact] public void ComplexBesselKResolvesOscillatoryIntegral() =>
        Close(new Complex(3.015242162785077, -31.062215312887222), Special.K(new Complex32(1, 5), 10));

    [Fact] public void BetaAvoidsIntermediateGammaOverflow() => Close(1.7681885473062026e-15, Special.Beta(20f, 30f), 1e-20);

    [Fact] public void HypergeometricSeriesAvoidsIntermediatePochhammerOverflow() =>
        Close(21.789416887313024, Special.Hypergeom(2f, 3f, 4f, .9f));

    [Fact]
    public void PadeSupportsDenominatorDegreeLargerThanNumeratorDegree()
    {
        var (p, q) = new Pade(1, 3).Compute(new[] { 1f, 1f, .5f, 1f / 6, 1f / 24 });
        Close(new[] { 1f, .25f }, p);
        Close(new[] { 1f, -.75f, .25f, -1f / 24 }, q);
    }

    [Theory] [InlineData(0f)] [InlineData(1f)]
    public void BilinearInterpolationIsLinearAlongGridEdges(float x) =>
        Close(x + .5, new Interpolation().Compute(new[] { 0f, 1f }, new[] { 0f, 1f }, new float[,] { { 0, 1 }, { 1, 2 } }, x, .5f));

    [Fact] public void ComplexVarianceUsesSquaredMagnitudes() => Close(new Complex(2, 0), new[] { new Complex32(0, 1), new Complex32(0, -1) }.Var());

    [Fact] public void ComplexVectorModulusCannotCancelNonzeroComponents() => Close(new Complex(Math.Sqrt(2), 0), new[] { new Complex32(1, 0), new Complex32(0, 1) }.Abs());

    [Theory] [InlineData(0)] [InlineData(5)]
    public void BinomialCertainSuccessHasUnitMass(int n) => Close(1, new Binomial(n, 1).Function(n));

    [Fact] public void PoissonMassIsFiniteAtItsModeForLambda100() => Close(.0398609968091471, new Poisson(100).Function(100));

    [Fact] public void BinomialMedianSatisfiesBothHalfProbabilityInequalities() => Close(1, new Binomial(2, .7f).Median);

    [Fact] public void PoissonMedianCanBeOneBelowLambdaOne() => Close(1, new Poisson(.9f).Median);

    [Fact] public void IirStabilityUsesTheSamePolynomialAsReaction() => Assert.False(new IIR(new[] { 1f }, new[] { .5f, -.75f }).Stability);

    [Fact] public void EvenNormalWindowIsSymmetric() { var w = new Normal(8).GetWindow(); Close(w, w.Reverse().ToArray()); }

    [Fact] public void EvenConfinedWindowIsSymmetric() { var w = new Confined(8).GetWindow(); Close(w, w.Reverse().ToArray()); }

    [Fact] public void XyzCanRepresentItsD65WhitePoint() => Close(1.089, XYZ.White.Z);

    [Fact] public void RybConversionPreservesWhite() { var color = RYB.FromRGB(255, 255, 255).ToRGB; Assert.Equal((byte)255, color.Red); Assert.Equal((byte)255, color.Green); Assert.Equal((byte)255, color.Blue); }

    [Fact] public void TanhAvoidsInfinityDividedByInfinity() => Close(1, Maths.Tanh(100f));

    [Fact] public void AsinhAvoidsCancellationOnNegativeArguments() => Close(Math.Asinh(-10000), Maths.Asinh(-10000f));

    [Theory] [InlineData(1e-20f)] [InlineData(1e20f)]
    public void ComplexDivisionIsInvariantToScaling(float scale) { var z = new Complex32(scale, scale); Close(Complex.One, z / z); }

    [Fact] public void ChebyshevUHasFiniteEndpointValue() => Close(3, Special.ChebyshevU(1f, 2));

    [Fact] public void ChebyshevPolynomialIsDefinedOutsideUnitInterval() => Close(7, Special.ChebyshevT(2f, 2));

    [Fact] public void PositiveOrderPochhammerOfZeroIsZero() => Close(0, Special.FactorialUp(0, 2));

    [Fact] public void NegativeRealBaseSupportsComplexExponent() => Close(Complex.ImaginaryOne, Maths.Pow(-1f, new Complex32(.5f, 0)));

    [Fact] public void InverseChiSquareEntropyAgreesWithInverseGamma() => Close(1 - Math.Log(2) + 2 * .57721566490153286, new InverseChiSquare(2).Entropy);
}
