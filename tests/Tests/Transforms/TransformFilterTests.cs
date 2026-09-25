using System.Numerics;
using UMapx.Core;
using UMapx.Transform;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Transform")]
public class TransformFilterTests
{
    [Theory]
    [InlineData(ThresholdMode.Abs)]
    [InlineData(ThresholdMode.Under)]
    [InlineData(ThresholdMode.Over)]
    public void ThresholdFiltersRespectEqualityAndSignedComponents(ThresholdMode mode)
    {
        var input = new[]
        {
            -2f,
            -1,
            -.5f,
            0,
            .5f,
            1,
            2
        };
        float threshold = 1;
        float F(float v) => mode switch
        {
            ThresholdMode.Abs => Math.Abs(v) < threshold ? 0 : v,
            ThresholdMode.Under => v < threshold ? 0 : v,
            _ => v > threshold ? 0 : v
        };
        var d = new ThresholdFilter(threshold, mode);
        var x = (float[])input.Clone();
        d.Apply(x);
        Close(input.Select(F).ToArray(), x);
        var a = new float[2, input.Length];
        var c = new Complex32[2, input.Length];
        var z = new Complex32[input.Length];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < input.Length; j++)
            {
                a[i, j] = input[j];
                c[i, j] = new(input[j], -.5f * input[j]);
                z[j] = c[i, j];
            }

        var expected = z.Select(v => mode == ThresholdMode.Abs ? (v.Abs < threshold ? Complex32.Zero : v) : new Complex32(F(v.Real), F(v.Imag))).ToArray();
        d.Apply(a);
        d.Apply(c);
        d.Apply(z);
        for (int j = 0; j < input.Length; j++)
        {
            Close((Complex)expected[j], z[j]);
            for (int i = 0; i < 2; i++)
            {
                Close(F(input[j]), a[i, j]);
                Close((Complex)expected[j], c[i, j]);
            }
        }
    }

    public static IEnumerable<object[]> FilterCases()
    {
        foreach (string name in new[]
        {
            "Guided",
            "Bilateral",
            "BilateralGrid",
            "Domain",
            "LocalLaplacian",
            "Laplacian"
        }

        )
            foreach (bool matrix in new[]
            {
                false,
                true
            }

            )
                foreach (bool complex in new[]
                {
                    false,
                    true
                }

                )
                    yield return new object[]
                    {
                        name,
                        matrix,
                        complex
                    };
    }

    [Theory, MemberData(nameof(FilterCases))]
    public void SmoothingAndDetailFiltersPreserveConstants(string name, bool matrix, bool complex)
    {
        IFilter d = name switch
        {
            "Guided" => new GuidedFilter(2),
            "Bilateral" => new BilateralFilter(2, .1f, 8),
            "BilateralGrid" => new BilateralGridFilter(2, .1f),
            "Domain" => new DomainTransformFilter(2, .1f),
            "LocalLaplacian" => new LocalLaplacianFilter(2, .1f, 5, 3),
            _ => new LaplacianPyramidFilter(new LaplacianPyramidTransform(3, 2))
        };
        var x = Enumerable.Repeat(.375f, 16).ToArray();
        var z = x.Select(v => new Complex32(v, .125f)).ToArray();
        var a = new float[8, 12];
        var c = new Complex32[8, 12];
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 12; j++)
            {
                a[i, j] = .375f;
                c[i, j] = new(.375f, .125f);
            }

        if (name == "LocalLaplacian" && complex)
        {
            if (matrix)
                Assert.Throws<NotSupportedException>(() => d.Apply(c));
            else
                Assert.Throws<NotSupportedException>(() => d.Apply(z));
            return;
        }

        if (matrix && complex)
        {
            d.Apply(c);
            Assert.All(c.Cast<Complex32>(), v => Close(new Complex(.375, .125), v, .001));
        }
        else if (matrix)
        {
            d.Apply(a);
            Assert.All(a.Cast<float>(), v => Close(.375, v, .001));
        }
        else if (complex)
        {
            d.Apply(z);
            Assert.All(z, v => Close(new Complex(.375, .125), v, .001));
        }
        else
        {
            d.Apply(x);
            Assert.All(x, v => Close(.375, v, .001));
        }

        if (name is "Guided" or "Bilateral" or "Domain")
        {
            // A small interior impulse must spread and lose height under smoothing.
            // Constant preservation alone would also accept an identity implementation.
            int center = matrix ? 4 * 12 + 6 : 8;
            x[8] += .05f;
            z[8] += new Complex32(.05f, .025f);
            a[4, 6] += .05f;
            c[4, 6] += new Complex32(.05f, .025f);
            Complex32[] response;
            if (matrix && complex)
            {
                d.Apply(c);
                response = c.Cast<Complex32>().ToArray();
            }
            else if (matrix)
            {
                d.Apply(a);
                response = a.Cast<float>().Select(v => new Complex32(v, 0)).ToArray();
            }
            else if (complex)
            {
                d.Apply(z);
                response = z;
            }
            else
            {
                d.Apply(x);
                response = x.Select(v => new Complex32(v, 0)).ToArray();
            }

            Assert.All(response, v =>
            {
                Assert.True(float.IsFinite(v.Real) && float.IsFinite(v.Imag));
                Assert.InRange(v.Real, .375f - 2e-6f, .425f + 2e-6f);
                if (complex)
                    Assert.InRange(v.Imag, .125f - 2e-6f, .15f + 2e-6f);
            });
            Assert.InRange(response[center].Real, .375f + 1e-5f, .425f - 1e-5f);
            Assert.Contains(Enumerable.Range(0, response.Length), i => i != center && response[i].Real > .375f + 1e-5f);
            if (complex)
            {
                Assert.InRange(response[center].Imag, .125f + 1e-5f, .15f - 1e-5f);
                Assert.Contains(Enumerable.Range(0, response.Length), i => i != center && response[i].Imag > .125f + 1e-5f);
            }
        }
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void WaveletFiltersWithZeroFactorPreserveTheirInput(bool complex, bool matrix)
    {
        var type = matrix ? (complex ? typeof(Complex32[,]) : typeof(float[,])) : (complex ? typeof(Complex32[]) : typeof(float[]));
        var input = (Array)MatrixTestSupport.Operand(type, 1, 16, 16);
        // Bilateral grids quantize intensities on [0,1]. Keep magnitudes inside that domain.
        for (int y = 0; y < (matrix ? 16 : 1); y++)
            for (int x = 0; x < 16; x++)
            {
                var v = MatrixTestSupport.Value(input, y, x) / 8;
                object sample = complex ? (object)(Complex32)v : (float)v.Real;
                if (matrix)
                    input.SetValue(sample, y, x);
                else
                    input.SetValue(sample, x);
            }

        foreach (IFilter filter in new IFilter[]
        {
            new WaveletFilter(new WaveletDecomposition(WaveletPacket.D4, 2), 0),
            new EdgeAvoidingWaveletFilter(new EdgeAvoidingWaveletDecomposition(levels: 2), 0)
        }

        )
        {
            var actual = (Array)input.Clone();
            typeof(IFilter).GetMethod("Apply", new[] { type })!.Invoke(filter, new object[] { actual });
            var expected = input.Cast<object>().ToArray();
            var result = actual.Cast<object>().ToArray();
            for (int i = 0; i < expected.Length; i++)
                MatrixTestSupport.Check(MatrixTestSupport.Value(expected[i]), result[i], 2e-4);
        }
    }

    [Theory]
    [InlineData(-1f)]
    [InlineData(0f)]
    [InlineData(1f)]
    public void ComplexGridFilteringCommutesWithPhaseRotationAndConjugation(float factor)
    {
        var filter = new BilateralGridFilter(2, .125f, factor);
        var x = Enumerable.Range(0, 24).Select(i => new Complex32((float)(.4 + .2 * Math.Sin(i)), (float)(.15 * Math.Cos(i)))).ToArray();
        var rotation = new Complex32(0, 1);
        var rotated = x.Select(v => v * rotation).ToArray();
        var conjugate = x.Select(v => new Complex32(v.Real, -v.Imag)).ToArray();
        filter.Apply(x);
        filter.Apply(rotated);
        filter.Apply(conjugate);
        for (int i = 0; i < x.Length; i++)
        {
            Close((Complex)(x[i] * rotation), rotated[i], 2e-6);
            Close(Complex.Conjugate((Complex)x[i]), conjugate[i], 2e-6);
        }

        var matrix = new Complex32[4, 6];
        var shifted = new Complex32[4, 6];
        for (int i = 0; i < 24; i++)
        {
            matrix[i / 6, i % 6] = x[i];
            shifted[i / 6, i % 6] = x[i] * rotation;
        }

        filter.Apply(matrix);
        filter.Apply(shifted);
        for (int i = 0; i < 24; i++)
            Close((Complex)(matrix[i / 6, i % 6] * rotation), shifted[i / 6, i % 6], 2e-6);
    }

    [Theory]
    [InlineData(ThresholdMode.Abs)]
    [InlineData(ThresholdMode.Under)]
    [InlineData(ThresholdMode.Over)]
    public void ComplexThresholdingUsesTheSameComponentsForEveryShape(ThresholdMode mode)
    {
        var input = new[]
        {
            new Complex32(-.5f, 2),
            new Complex32(2, -.5f),
            new Complex32(.5f, -.5f),
            new Complex32(.2f, .2f)
        };
        var vector = (Complex32[])input.Clone();
        var matrix = new Complex32[2, 2];
        for (int i = 0; i < 4; i++)
            matrix[i / 2, i % 2] = input[i];
        var filter = new ThresholdFilter(.5f, mode);
        filter.Apply(vector);
        filter.Apply(matrix);
        double F(double x) => mode == ThresholdMode.Under ? (x < .5 ? 0 : x) : (x > .5 ? 0 : x);
        for (int i = 0; i < 4; i++)
        {
            Complex expected = mode == ThresholdMode.Abs ? (input[i].Abs < .5 ? Complex.Zero : (Complex)input[i]) : new Complex(F(input[i].Real), F(input[i].Imag));
            Close(expected, vector[i], 0, 0);
            Close(expected, matrix[i / 2, i % 2], 0, 0);
        }
    }

    [Theory]
    [InlineData(0f)]
    [InlineData(.375f)]
    [InlineData(1f)]
    public void LocalLaplacianPreservesConstantsAtAllLevelsAndDegenerateParameters(float value)
    {
        foreach (int size in new[]
        {
            1,
            2,
            8,
            17
        }

        )
            foreach (int steps in new[]
            {
                0,
                1,
                3,
                7
            }

            )
                foreach (float sigma in new[]
                {
                    0f,
                    .1f
                }

                )
                {
                    var filter = new LocalLaplacianFilter(2, sigma, steps, 10, -1);
                    var x = Enumerable.Repeat(value, size).ToArray();
                    var matrix = new float[size, size];
                    for (int y = 0; y < size; y++)
                        for (int z = 0; z < size; z++)
                            matrix[y, z] = value;
                    filter.Apply(x);
                    filter.Apply(matrix);
                    Assert.All(x, v => Close(value, v, 2e-6, 0));
                    foreach (float v in matrix)
                        Close(value, v, 2e-6, 0);
                }
    }

    public static IEnumerable<object[]> LocalCases()
    {
        foreach (int steps in new[]
        {
            3,
            7
        }

        )
            foreach (float sigma in new[]
            {
                .05f,
                .2f
            }

            )
                foreach (float factor in new[]
                {
                    -1f,
                    1f
                }

                )
                    foreach (int radius in new[]
                    {
                        2,
                        3
                    }

                    )
                        yield return new object[]
                        {
                            steps,
                            sigma,
                            factor,
                            radius
                        };
    }

    [Theory, MemberData(nameof(LocalCases))]
    public void LocalLaplacianMatchesDirectRemappingWithoutLookupTables(int steps, float sigma, float factor, int radius)
    {
        var input = Enumerable.Range(0, 17).Select(i => (float)(.5 + .45 * Math.Sin(.41 * i))).ToArray();
        var expected = DirectLocalLaplacian(input, steps, sigma, factor, radius);
        var actual = (float[])input.Clone();
        var image = new float[32, input.Length];
        for (int y = 0; y < 32; y++)
            for (int x = 0; x < input.Length; x++)
                image[y, x] = input[x];
        var filter = new LocalLaplacianFilter(radius, sigma, steps, 3, factor);
        filter.Apply(actual);
        filter.Apply(image);
        // Budget the interpolation of the 256-entry LUT against direct double exponentials.
        for (int x = 0; x < input.Length; x++)
        {
            Close(expected[x], actual[x], 3e-4, 0);
            for (int y = 0; y < 32; y++)
                Close(expected[x], image[y, x], 3e-4, 0);
        }

        Assert.Contains(Enumerable.Range(0, input.Length), i => Math.Abs(actual[i] - input[i]) > 1e-4);
    }

    /// <summary>Evaluates the sampled local-Laplacian definition with independent double pyramids and no lookup tables.</summary>
    /// <param name = "input">Unit-interval signal with at least eight samples.</param>
    /// <param name = "steps">Number of intensity intervals.</param>
    /// <param name = "sigma">Positive remapping width.</param>
    /// <param name = "factor">Signed detail strength.</param>
    /// <param name = "radius">Clipped averaging window length.</param>
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
            var p = new[]
            {
                x,
                Array.Empty<double>(),
                Array.Empty<double>()
            };
            for (int l = 1; l < 3; l++)
                p[l] = Mean(p[l - 1].Where((_, i) => i % 2 == 0).ToArray());
            return p;
        }

        double[][] Laplacian(double[][] g) => new[]
        {
            g[0].Zip(Up(g[1]), (a, b) => a - b).ToArray(),
            g[1].Zip(Up(g[2]), (a, b) => a - b).ToArray(),
            (double[])g[2].Clone()
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
        for (int level = 1; level >= 0; level--)
            result = output[level].Zip(Up(result), (a, b) => a + b).ToArray();
        return result;
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void ComplexGridMatchesAnIndependentSumOfSampleInfluences(bool matrix)
    {
        int height = matrix ? 4 : 1, width = 6;
        var input = new Complex32[height, width];
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
                input[y, x] = new Complex32((float)(.4 * Math.Sin(x + 2 * y)), (float)(.5 * Math.Cos(2 * x - y)));
        var actual = (Complex32[,])input.Clone();
        var vector = Enumerable.Range(0, width).Select(x => input[0, x]).ToArray();
        var filter = new BilateralGridFilter(2, .125f, -1);
        if (matrix)
            filter.Apply(actual);
        else
            filter.Apply(vector);
        int gridZ = 9, gridX = width / 2 + 2, gridY = matrix ? height / 2 + 2 : 1;
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                double z = input[y, x].Abs / .125, gx = x / 2.0, gy = matrix ? y / 2.0 : 0;
                Complex numerator = 0;
                double denominator = 0;
                for (int sy = 0; sy < height; sy++)
                    for (int sx = 0; sx < width; sx++)
                    {
                        double influence = 0;
                        for (int cz = 0; cz <= 1; cz++)
                            for (int cx = 0; cx <= 1; cx++)
                                for (int cy = 0; cy <= (matrix ? 1 : 0); cy++)
                                {
                                    int zi = Math.Clamp((int)z + cz, 0, gridZ - 1), xi = Math.Clamp((int)gx + cx, 0, gridX - 1), yi = Math.Clamp((int)gy + cy, 0, gridY - 1);
                                    double weight = (1 - Math.Abs(z - (int)z - cz)) * (1 - Math.Abs(gx - (int)gx - cx)) * (matrix ? 1 - Math.Abs(gy - (int)gy - cy) : 1);
                                    for (int dz = -1; dz <= 1; dz++)
                                        for (int dx = -1; dx <= 1; dx++)
                                            for (int dy = matrix ? -1 : 0; dy <= (matrix ? 1 : 0); dy++)
                                                if (Math.Clamp(zi + dz, 0, gridZ - 1) == (int)(input[sy, sx].Abs / .125f) && Math.Clamp(xi + dx, 0, gridX - 1) == sx / 2 && Math.Clamp(yi + dy, 0, gridY - 1) == sy / 2)
                                                    influence += weight;
                                }

                        numerator += influence * (Complex)input[sy, sx];
                        denominator += influence;
                    }

                Close(numerator / denominator, matrix ? actual[y, x] : vector[x], 3e-6, 0);
            }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void FrequencyMasksPreserveTheRequestedBins(bool complex)
    {
        var matrix = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, 7, 9);
        var expected = (Array)matrix.Clone();
        var vector = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 2, 1, 9);
        var vcopy = (Array)vector.Clone();
        var filter = new FrequencyFilter(-2, 2);
        if (complex)
            filter.Apply((Complex32[])vector);
        else
            filter.Apply((float[])vector);
        for (int i = 0; i < 9; i++)
            Check(i >= 2 && i <= 6 ? Value(vcopy, 0, i) : Complex.Zero, vector.GetValue(i)!);
        filter.FrequencyRange = new RangeInt(1, 2);
        if (complex)
            filter.Apply((Complex32[,])matrix);
        else
            filter.Apply((float[,])matrix);
        // Matrix radii are quantized to integer bins by this API.
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 9; j++)
            {
                int squared = (i - 3) * (i - 3) + (j - 4) * (j - 4);
                Check(squared >= 1 && squared < 9 ? Value(expected, i, j) : Complex.Zero, matrix.GetValue(i, j)!);
            }
    }
}
