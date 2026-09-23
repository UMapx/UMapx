using System.Numerics;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

internal static class NumericAssert
{
    internal static void Close(double expected, double actual, double absolute = 2e-6, double relative = 2e-5)
    {
        double tolerance = absolute + relative * Math.Abs(expected);
        Assert.True(double.IsFinite(actual) && Math.Abs(actual - expected) <= tolerance,
            FormattableString.Invariant($"Expected {expected:G17}; actual {actual:G17}; tolerance {tolerance:G6}."));
    }

    internal static void Close(Complex expected, Complex32 actual, double absolute = 2e-6, double relative = 2e-5)
    {
        double error = Complex.Abs((Complex)actual - expected);
        double tolerance = absolute + relative * expected.Magnitude;
        Assert.True(double.IsFinite(error) && error <= tolerance,
            FormattableString.Invariant($"Expected {expected}; actual ({actual.Real:G9}, {actual.Imag:G9}); error {error:G6}."));
    }

    internal static void Close(float[] expected, float[] actual, double tolerance = 1e-4)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (int i = 0; i < expected.Length; i++) Close(expected[i], actual[i], tolerance, 0);
    }

    internal static void Close(float[,] expected, float[,] actual, double tolerance = 1e-4)
    {
        Assert.Equal(expected.GetLength(0), actual.GetLength(0));
        Assert.Equal(expected.GetLength(1), actual.GetLength(1));
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++) Close(expected[i, j], actual[i, j], tolerance, 0);
    }

    // Independent double accumulation: do not verify Matrice.Dot with Matrice.Dot.
    internal static float[,] Product(float[,] a, float[,] b)
    {
        Assert.Equal(a.GetLength(1), b.GetLength(0));
        var result = new float[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++)
            {
                double sum = 0;
                for (int k = 0; k < a.GetLength(1); k++) sum += (double)a[i, k] * b[k, j];
                result[i, j] = (float)sum;
            }
        return result;
    }

    internal static float[,] Transpose(float[,] a)
    {
        var result = new float[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) result[j, i] = a[i, j];
        return result;
    }

    internal static float[,] Diagonal(float[] values)
    {
        var result = new float[values.Length, values.Length];
        for (int i = 0; i < values.Length; i++) result[i, i] = values[i];
        return result;
    }

    internal static float[,] Matrix(int m, int n)
    {
        var random = new Random(731911 + 17 * m + n);
        var a = new float[m, n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) a[i, j] = (float)(2 * random.NextDouble() - 1);
        return a;
    }
}
