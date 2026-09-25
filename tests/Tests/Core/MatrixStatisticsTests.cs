using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixStatisticsTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(5)]
    [InlineData(17)]
    public void RealStatisticsMatchSampleDefinitions(int n)
    {
        var x = Enumerable.Range(0, n).Select(i => (float)(Math.Sin(i * .7) + i * .2)).ToArray();
        double mean = x.Average(v => (double)v), variance = x.Sum(v => Math.Pow(v - mean, 2)) / (n - 1);
        Close(x.Sum(v => (double)v), Matrice.Sum(x));
        Close(mean, Matrice.Mean(x));
        Close(variance, Matrice.Var(x));
        Close(Math.Sqrt(variance), Matrice.StnDev(x));
        Close(variance, Matrice.Cov(x));
        var sorted = x.OrderBy(v => v).ToArray();
        Close(sorted[0], Matrice.Min(x));
        Close(sorted[^1], Matrice.Max(x));
        var a = new float[n, 3];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < 3; j++)
                a[i, j] = x[i] * (j + 1) + j;
        var means = a.Mean();
        var vars = a.Var();
        var deviations = a.StnDev();
        var sums = a.Sum();
        var cov = a.Cov();
        for (int j = 0; j < 3; j++)
        {
            Close(mean * (j + 1) + j, means[j]);
            Close(variance * (j + 1) * (j + 1), vars[j]);
            Close(Math.Sqrt(variance) * (j + 1), deviations[j]);
            Close((mean * (j + 1) + j) * n, sums[j]);
            for (int k = 0; k < 3; k++)
                Close(variance * (j + 1) * (k + 1), cov[j, k]);
        }

        var p = Enumerable.Range(1, n).Select(i => (float)(2.0 * i / (n * (n + 1)))).ToArray();
        Close(-p.Sum(v => v * Math.Log2(v)), Matrice.Entropy(p));
        var normalized = x.Normalized();
        for (int i = 0; i < n; i++)
            Close((x[i] - sorted[0]) / (sorted[^1] - sorted[0]), normalized[i]);
    }

    [Theory]
    [InlineData("Variance")]
    [InlineData("CovarianceVector")]
    [InlineData("CovarianceMatrix")]
    [InlineData("Norm")]
    public void ComplexStatisticsRespectHermitianInnerProducts(string operation)
    {
        var x = new[]
        {
            new Complex32(1, 2),
            new Complex32(-2, 1),
            new Complex32(3, -4),
            new Complex32(2, 3)
        };
        Complex mean = x.Aggregate(Complex.Zero, (s, v) => s + (Complex)v) / x.Length;
        double variance = x.Sum(v => Complex.Abs((Complex)v - mean) * Complex.Abs((Complex)v - mean)) / (x.Length - 1);
        if (operation == "Variance")
            Check(variance, x.Var());
        else if (operation == "CovarianceVector")
            Check(variance, x.Cov());
        else if (operation == "Norm")
            Check(Math.Sqrt(x.Sum(v => Complex.Abs((Complex)v) * Complex.Abs((Complex)v))), x.Abs());
        else
        {
            var a = new Complex32[x.Length, 2];
            for (int i = 0; i < x.Length; i++)
            {
                a[i, 0] = x[i];
                a[i, 1] = new Complex32(0, 1) * x[i];
            }

            var c = a.Cov();
            Check(variance, c[0, 0]);
            Check(variance, c[1, 1]);
            Check(Complex.ImaginaryOne * variance, c[0, 1]);
            Check(-Complex.ImaginaryOne * variance, c[1, 0]);
        }
    }

    public static IEnumerable<object[]> StatisticsCases()
    {
        foreach (int n in new[]
        {
            2,
            3,
            17
        }

        )
            foreach (float scale in new[]
            {
                1e-20f,
                1f,
                1e10f,
                1e20f
            }

            )
                foreach (float offset in new[]
                {
                    0f,
                    10000f
                }

                )
                    yield return new object[]
                    {
                        n,
                        scale,
                        offset
                    };
    }

    [Theory]
    [MemberData(nameof(StatisticsCases))]
    public void HermitianMomentsAccumulateInDoubleAndReturnRealValues(int n, float scale, float offset)
    {
        var v = Enumerable.Range(0, n).Select(i => new Complex32(scale * (offset + i * .7f), scale * (-offset + (i % 3) * 1.3f))).ToArray();
        var other = v.Select(z => new Complex32(-z.Imag, z.Real)).ToArray();
        Complex mean = v.Aggregate(Complex.Zero, (s, z) => s + (Complex)z) / n;
        double variance = v.Sum(z => Complex.Abs((Complex)z - mean) * Complex.Abs((Complex)z - mean)) / (n - 1);
        double norm2 = v.Sum(z => (double)z.Real * z.Real + (double)z.Imag * z.Imag);
        double error = v.Select((z, i) => Complex.Abs((Complex)z - (Complex)other[i])).Sum(x => x * x) / (n - 1);
        RealComplex(variance, v.Var());
        RealComplex(variance, v.Cov());
        RealComplex(Math.Sqrt(variance), v.StnDev());
        RealComplex(norm2, v.Abs(true));
        RealComplex(Math.Sqrt(norm2), v.Abs());
        RealComplex(error, v.Var(other));
        RealComplex(Math.Sqrt(error), v.StnDev(other));
        var matrix = new Complex32[n, 2];
        var rotated = new Complex32[n, 2];
        for (int i = 0; i < n; i++)
        {
            matrix[i, 0] = v[i];
            matrix[i, 1] = other[i];
            rotated[i, 0] = other[i];
            rotated[i, 1] = -v[i];
        }

        var vars = matrix.Var();
        var std = matrix.StnDev();
        var errors = matrix.Var(rotated);
        var rms = matrix.StnDev(rotated);
        for (int j = 0; j < 2; j++)
        {
            RealComplex(variance, vars[j]);
            RealComplex(Math.Sqrt(variance), std[j]);
            RealComplex(error, errors[j]);
            RealComplex(Math.Sqrt(error), rms[j]);
        }

        var covariance = matrix.Cov();
        RealComplex(variance, covariance[0, 0]);
        RealComplex(variance, covariance[1, 1]);
        RealComplex(variance, new Complex32(covariance[0, 1].Imag, covariance[0, 1].Real));
        RealComplex(variance, new Complex32(-covariance[1, 0].Imag, covariance[1, 0].Real));
        var rowNorms = matrix.Abs();
        for (int i = 0; i < n; i++)
            RealComplex(Math.Sqrt(2) * Complex.Abs((Complex)v[i]), rowNorms[i]);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(1)]
    public void SampleVarianceIsUndefinedWithFewerThanTwoObservations(int n)
    {
        var v = new Complex32[n];
        Assert.True(float.IsNaN(v.Var().Real));
        Assert.True(float.IsNaN(v.Cov().Real));
        Assert.True(float.IsNaN(v.StnDev().Real));
        Assert.Equal(0, v.Var().Imag);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(17)]
    public void ComplexCovarianceMatchesCenteredOuterProductsAndIsPositiveSemidefinite(int n)
    {
        var a = new Complex32[n, 3];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < 3; j++)
                a[i, j] = new Complex32((float)(Math.Sin(i * (j + 1)) + i), (float)Math.Cos(i + j * .7));
        var means = new Complex[3];
        for (int j = 0; j < 3; j++)
            for (int i = 0; i < n; i++)
                means[j] += (Complex)a[i, j] / n;
        var covariance = a.Cov();
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
            {
                Complex expected = 0;
                for (int i = 0; i < n; i++)
                    expected += Complex.Conjugate((Complex)a[i, j] - means[j]) * ((Complex)a[i, k] - means[k]) / (n - 1);
                Close(expected, covariance[j, k]);
                Assert.Equal(covariance[j, k].Real, covariance[k, j].Real);
                Assert.Equal(covariance[j, k].Imag, -covariance[k, j].Imag);
            }

        Complex[] weights =
        {
            new(1, 2),
            new(-2, .5),
            new(0, -1)
        };
        Complex quadratic = 0;
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                quadratic += Complex.Conjugate(weights[j]) * (Complex)covariance[j, k] * weights[k];
        Assert.True(quadratic.Real >= 0);
        Close(0, quadratic.Imaginary);
    }

    [Fact]
    public void ComplexVarianceUsesSquaredMagnitudes() => Close(new Complex(2, 0), new[] { new Complex32(0, 1), new Complex32(0, -1) }.Var());
    [Fact]
    public void ComplexVectorModulusCannotCancelNonzeroComponents() => Close(new Complex(Math.Sqrt(2), 0), new[] { new Complex32(1, 0), new Complex32(0, 1) }.Abs());
}
