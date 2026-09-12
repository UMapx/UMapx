using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using UMapx.Distribution;
using UMapx.Response;
using UMapx.Transform;
using Xunit;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Transform")]
public class TransformFilterRepairTests
{
    public static IEnumerable<object[]> HankelCases()
    {
        using var stream = typeof(TransformFilterRepairTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.meyer-hankel.json");
        using var data = JsonDocument.Parse(stream!);
        foreach (var row in data.RootElement.GetProperty("hankel").EnumerateArray())
            yield return new object[] { row.GetProperty("order").GetInt32(), row.GetProperty("size").GetInt32(), row.GetProperty("values").GetRawText() };
    }

    [Theory, MemberData(nameof(HankelCases))]
    public void HankelUsesConsecutivePositiveRootsAtSmallAndLargeOrders(int order, int size, string reference)
    {
        using var data = JsonDocument.Parse(reference);
        var matrix = HankelTransform.Matrix(size, order);
        for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        {
            Close(data.RootElement[i][j].GetDouble(), matrix[i, j], 3e-5, 3e-5);
            Assert.Equal(matrix[i, j], matrix[j, i]);
        }
    }

    public static IEnumerable<object[]> PoleCases()
    {
        foreach (float a0 in new[] { -.9f, 0f, .9f })
        foreach (float pole in new[] { -1.05f, -.75f, 0f, .5f, .98f, 1f, 1.05f })
        foreach (bool paired in new[] { false, true })
            yield return new object[] { a0, pole, paired };
    }

    [Theory, MemberData(nameof(PoleCases))]
    public void IirStabilityMatchesKnownRealAndComplexPoles(float a0, float radius, bool paired)
    {
        double[] polynomial = paired ? new[] { 1.0, -(double)radius, (double)radius * radius } : new[] { 1.0, -(double)radius };
        // The quadratic has conjugate poles radius * exp(+/- i*pi/3).
        var feedback = new float[polynomial.Length];
        feedback[0] = a0;
        for (int i = 1; i < feedback.Length; i++) feedback[i] = (float)(-(1.0 - a0) * polynomial[i]);
        var filter = new IIR { A = feedback, B = new[] { 1f } };
        Assert.Equal(Math.Abs(radius) < 1, filter.Stability);
    }

    [Fact]
    public void IirHandlesRepeatedPolesZeroOrderAndUndefinedFeedback()
    {
        // (z - 1/2)^4, all roots strictly inside the unit disk.
        Assert.True(new IIR { A = new[] { 0f, 2f, -1.5f, .5f, -.0625f } }.Stability);
        Assert.True(new IIR { A = new[] { 2f } }.Stability);
        Assert.False(new IIR { A = new[] { 1f } }.Stability);
        Assert.False(new IIR { A = new[] { 0f, float.NaN } }.Stability);
        Assert.False(new IIR { A = new[] { float.PositiveInfinity } }.Stability);
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(3)] [InlineData(7)]
    [InlineData(8)] [InlineData(17)] [InlineData(32)]
    public void LaplacianPyramidsRecoverComplexVectorsAndThinMatrices(int size)
    {
        foreach (int levels in new[] { 1, 3, 10 })
        foreach (int radius in new[] { 0, 1, 2 })
        {
            var d = new LaplacianPyramidTransform(levels, radius);
            var x = Enumerable.Range(0, size).Select(i => new Complex32((float)Math.Sin(i), (float)Math.Cos(.3 * i))).ToArray();
            var result = d.Backward(d.Forward(x));
            Close(x.Select(v => v.Real).ToArray(), d.Backward(d.Forward(x.Select(v => v.Real).ToArray())), 2e-6);
            for (int i = 0; i < size; i++) Close((Complex)x[i], result[i], 2e-6);
            var thin = new Complex32[1, size];
            for (int i = 0; i < size; i++) thin[0, i] = x[i];
            var restored = d.Backward(d.Forward(thin));
            for (int i = 0; i < size; i++) Close((Complex)x[i], restored[0, i], 2e-6);
        }
    }

    [Theory]
    [InlineData(-1f)] [InlineData(0f)] [InlineData(1f)]
    public void ComplexGridFilteringCommutesWithPhaseRotationAndConjugation(float factor)
    {
        var filter = new BilateralGridFilter(2, .125f, factor);
        var x = Enumerable.Range(0, 24).Select(i => new Complex32((float)(.4 + .2 * Math.Sin(i)), (float)(.15 * Math.Cos(i)))).ToArray();
        var rotation = new Complex32(0, 1);
        var rotated = x.Select(v => v * rotation).ToArray();
        var conjugate = x.Select(v => new Complex32(v.Real, -v.Imag)).ToArray();
        filter.Apply(x); filter.Apply(rotated); filter.Apply(conjugate);
        for (int i = 0; i < x.Length; i++)
        {
            Close((Complex)(x[i] * rotation), rotated[i], 2e-6);
            Close(Complex.Conjugate((Complex)x[i]), conjugate[i], 2e-6);
        }
        var matrix = new Complex32[4, 6]; var shifted = new Complex32[4, 6];
        for (int i = 0; i < 24; i++) { matrix[i / 6, i % 6] = x[i]; shifted[i / 6, i % 6] = x[i] * rotation; }
        filter.Apply(matrix); filter.Apply(shifted);
        for (int i = 0; i < 24; i++) Close((Complex)(matrix[i / 6, i % 6] * rotation), shifted[i / 6, i % 6], 2e-6);
    }

    [Theory]
    [InlineData(ThresholdMode.Abs)] [InlineData(ThresholdMode.Under)] [InlineData(ThresholdMode.Over)]
    public void ComplexThresholdingUsesTheSameComponentsForEveryShape(ThresholdMode mode)
    {
        var input = new[] { new Complex32(-.5f, 2), new Complex32(2, -.5f), new Complex32(.5f, -.5f), new Complex32(.2f, .2f) };
        var vector = (Complex32[])input.Clone(); var matrix = new Complex32[2, 2];
        for (int i = 0; i < 4; i++) matrix[i / 2, i % 2] = input[i];
        var filter = new ThresholdFilter(.5f, mode);
        filter.Apply(vector); filter.Apply(matrix);
        double F(double x) => mode == ThresholdMode.Under ? (x < .5 ? 0 : x) : (x > .5 ? 0 : x);
        for (int i = 0; i < 4; i++)
        {
            Complex expected = mode == ThresholdMode.Abs ? (input[i].Abs < .5 ? Complex.Zero : (Complex)input[i]) : new Complex(F(input[i].Real), F(input[i].Imag));
            Close(expected, vector[i], 0, 0); Close(expected, matrix[i / 2, i % 2], 0, 0);
        }
    }

    [Theory]
    [InlineData(.125f)] [InlineData(.5f)] [InlineData(2f)]
    public void ConeTimeKernelHasTheMassOfItsFourierValueAtZero(float tau)
    {
        var kernel = new ConeShape(.2f);
        foreach (float signed in new[] { -tau, tau })
        {
            double sum = 0;
            for (int i = 0; i < 1000; i++) sum += kernel.Distribution((float)(-tau / 2.0 + (i + .5) * tau / 1000), signed) * tau / 1000;
            Close(kernel.Function(0, signed), sum, 2e-6, 2e-6);
            Close(0, kernel.Distribution(tau, signed), 0, 0);
        }
    }

    [Theory]
    [InlineData(0f)] [InlineData(.375f)] [InlineData(1f)]
    public void LocalLaplacianPreservesConstantsAtAllLevelsAndDegenerateParameters(float value)
    {
        foreach (int size in new[] { 1, 2, 8, 17 })
        foreach (int steps in new[] { 0, 1, 3, 7 })
        foreach (float sigma in new[] { 0f, .1f })
        {
            var filter = new LocalLaplacianFilter(2, sigma, steps, 10, -1);
            var x = Enumerable.Repeat(value, size).ToArray(); var matrix = new float[size, size];
            for (int y = 0; y < size; y++) for (int z = 0; z < size; z++) matrix[y, z] = value;
            filter.Apply(x); filter.Apply(matrix);
            Assert.All(x, v => Close(value, v, 2e-6, 0));
            foreach (float v in matrix) Close(value, v, 2e-6, 0);
        }
    }

    public static IEnumerable<object[]> LocalCases()
    {
        foreach (int steps in new[] { 3, 7 })
        foreach (float sigma in new[] { .05f, .2f })
        foreach (float factor in new[] { -1f, 1f })
        foreach (int radius in new[] { 2, 3 }) yield return new object[] { steps, sigma, factor, radius };
    }

    [Theory, MemberData(nameof(LocalCases))]
    public void LocalLaplacianMatchesDirectRemappingWithoutLookupTables(int steps, float sigma, float factor, int radius)
    {
        var input = Enumerable.Range(0, 17).Select(i => (float)(.5 + .45 * Math.Sin(.41 * i))).ToArray();
        var expected = DirectLocalLaplacian(input, steps, sigma, factor, radius);
        var actual = (float[])input.Clone();
        var image = new float[32, input.Length];
        for (int y = 0; y < 32; y++) for (int x = 0; x < input.Length; x++) image[y, x] = input[x];
        var filter = new LocalLaplacianFilter(radius, sigma, steps, 3, factor);
        filter.Apply(actual); filter.Apply(image);
        // Budget the interpolation of the 256-entry LUT against direct double exponentials.
        for (int x = 0; x < input.Length; x++)
        {
            Close(expected[x], actual[x], 3e-4, 0);
            for (int y = 0; y < 32; y++) Close(expected[x], image[y, x], 3e-4, 0);
        }
        Assert.Contains(Enumerable.Range(0, input.Length), i => Math.Abs(actual[i] - input[i]) > 1e-4);
    }

    /// <summary>Evaluates the sampled local-Laplacian definition with independent double pyramids and no lookup tables</summary>
    /// <param name="input">Unit-interval signal with at least eight samples.</param>
    /// <param name="steps">Number of intensity intervals.</param>
    /// <param name="sigma">Positive remapping width.</param>
    /// <param name="factor">Signed detail strength.</param>
    /// <param name="radius">Clipped averaging window length.</param>
    /// <returns>The reconstruction of three modified detail/base levels.</returns>
    private static double[] DirectLocalLaplacian(float[] input, int steps, double sigma, double factor, int radius)
    {
        double[] Mean(double[] x) => Enumerable.Range(0, x.Length).Select(i =>
        {
            int first = Math.Max(0, i - radius / 2), last = Math.Min(x.Length - 1, i + (radius - 1) / 2);
            return Enumerable.Range(first, last - first + 1).Average(j => x[j]);
        }).ToArray();
        double[] Up(double[] x) => Mean(x.SelectMany(v => new[] { v, v }).ToArray());
        double[][] Gaussian(double[] x)
        {
            var p = new[] { x, Array.Empty<double>(), Array.Empty<double>() };
            for (int l = 1; l < 3; l++) p[l] = Mean(p[l - 1].Where((_, i) => i % 2 == 0).ToArray());
            return p;
        }
        double[][] Laplacian(double[][] g) => new[]
        {
            g[0].Zip(Up(g[1]), (a, b) => a - b).ToArray(),
            g[1].Zip(Up(g[2]), (a, b) => a - b).ToArray(), (double[])g[2].Clone()
        };
        var gaussian = Gaussian(input.Select(v => (double)v).ToArray());
        var output = Laplacian(gaussian);
        for (int sample = 0; sample <= steps; sample++)
        {
            double center = sample / (double)steps;
            var delta = Laplacian(Gaussian(input.Select(v => factor * (v - center) * Math.Exp(-Math.Pow(v - center, 2) / (2 * sigma * sigma))).ToArray()));
            for (int level = 0; level < 2; level++)
            for (int i = 0; i < output[level].Length; i++)
                output[level][i] += Math.Max(0, 1 - Math.Abs(gaussian[level][i] - center) * steps) * delta[level][i];
        }
        var result = output[2];
        for (int level = 1; level >= 0; level--) result = output[level].Zip(Up(result), (a, b) => a + b).ToArray();
        return result;
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void ComplexGridMatchesAnIndependentSumOfSampleInfluences(bool matrix)
    {
        int height = matrix ? 4 : 1, width = 6;
        var input = new Complex32[height, width];
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
            input[y, x] = new Complex32((float)(.4 * Math.Sin(x + 2 * y)), (float)(.5 * Math.Cos(2 * x - y)));
        var actual = (Complex32[,])input.Clone();
        var vector = Enumerable.Range(0, width).Select(x => input[0, x]).ToArray();
        var filter = new BilateralGridFilter(2, .125f, -1);
        if (matrix) filter.Apply(actual); else filter.Apply(vector);
        int gridZ = 9, gridX = width / 2 + 2, gridY = matrix ? height / 2 + 2 : 1;
        for (int y = 0; y < height; y++) for (int x = 0; x < width; x++)
        {
            double z = input[y, x].Abs / .125, gx = x / 2.0, gy = matrix ? y / 2.0 : 0;
            Complex numerator = 0; double denominator = 0;
            for (int sy = 0; sy < height; sy++) for (int sx = 0; sx < width; sx++)
            {
                double influence = 0;
                for (int cz = 0; cz <= 1; cz++) for (int cx = 0; cx <= 1; cx++) for (int cy = 0; cy <= (matrix ? 1 : 0); cy++)
                {
                    int zi = Math.Clamp((int)z + cz, 0, gridZ - 1), xi = Math.Clamp((int)gx + cx, 0, gridX - 1), yi = Math.Clamp((int)gy + cy, 0, gridY - 1);
                    double weight = (1 - Math.Abs(z - (int)z - cz)) * (1 - Math.Abs(gx - (int)gx - cx)) * (matrix ? 1 - Math.Abs(gy - (int)gy - cy) : 1);
                    for (int dz = -1; dz <= 1; dz++) for (int dx = -1; dx <= 1; dx++) for (int dy = matrix ? -1 : 0; dy <= (matrix ? 1 : 0); dy++)
                        if (Math.Clamp(zi + dz, 0, gridZ - 1) == (int)(input[sy, sx].Abs / .125f) &&
                            Math.Clamp(xi + dx, 0, gridX - 1) == sx / 2 && Math.Clamp(yi + dy, 0, gridY - 1) == sy / 2) influence += weight;
                }
                numerator += influence * (Complex)input[sy, sx]; denominator += influence;
            }
            Close(numerator / denominator, matrix ? actual[y, x] : vector[x], 3e-6, 0);
        }
    }
}
