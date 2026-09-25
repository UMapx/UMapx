using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using UMapx.Transform;
using UMapx.Wavelet;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Transform")]
public class TransformTests
{
    internal static readonly string[] Names =
    {
        "Cosine",
        "FastCosine",
        "Sine",
        "FastSine",
        "Chebyshev",
        "FastChebyshev",
        "Hartley",
        "FastHartley",
        "HartleyDirect",
        "FastHartleyDirect",
        "WalshHadamard",
        "FastWalshHadamard",
        "Fourier",
        "FastFourier",
        "Laplace",
        "FastLaplace",
        "Delta"
    };
    internal static ITransform Create(string name, Direction direction = Direction.Vertical, bool normalized = true) => name switch
    {
        "Cosine" => new CosineTransform(direction),
        "FastCosine" => new FastCosineTransform(direction),
        "Sine" => new SineTransform(direction),
        "FastSine" => new FastSineTransform(direction),
        "Chebyshev" => new ChebyshevTransform(direction),
        "FastChebyshev" => new FastChebyshevTransform(direction),
        "Hartley" => new HartleyTransform(normalized, SpectrumType.Fourier, direction),
        "FastHartley" => new FastHartleyTransform(normalized, SpectrumType.Fourier, direction),
        "HartleyDirect" => new HartleyTransform(normalized, SpectrumType.Hartley, direction),
        "FastHartleyDirect" => new FastHartleyTransform(normalized, SpectrumType.Hartley, direction),
        "WalshHadamard" => new WalshHadamardTransform(normalized, direction),
        "FastWalshHadamard" => new FastWalshHadamardTransform(normalized, direction),
        "Fourier" => new FourierTransform(normalized, direction),
        "FastFourier" => new FastFourierTransform(normalized, direction),
        "Laplace" => new LaplaceTransform(.07f, normalized, direction),
        "FastLaplace" => new FastLaplaceTransform(.07f, normalized, direction),
        _ => new DeltaTransform(direction)
    };
    private static bool ComplexOnly(string name) => name.Contains("Fourier") || name.Contains("Laplace");
    private static bool Scaled(string name) => name.Contains("Hartley") || name.Contains("Walsh") || ComplexOnly(name);
    private static Complex Weight(string name, int k, int j, int n, bool normalized)
    {
        if (name.Contains("Cosine"))
            return Math.Cos(Math.PI * (j + .5) * k / n) * Math.Sqrt((k == 0 ? 1 : 2.0) / n);
        if (name.Contains("Sine"))
            return Math.Sin(Math.PI * (j + 1) * (k + 1) / (n + 1)) * Math.Sqrt(2.0 / (n + 1));
        if (name.Contains("Chebyshev"))
            return n == 1 ? 1 : Math.Cos(Math.PI * k * j / (n - 1)) * Math.Sqrt(2.0 / (n - 1)) * (k == 0 || k == n - 1 ? 1 / Math.Sqrt(2) : 1) * (j == 0 || j == n - 1 ? 1 / Math.Sqrt(2) : 1);
        double scale = normalized ? Math.Sqrt(n) : 1;
        if (name.Contains("Hartley"))
            return (Math.Cos(2 * Math.PI * j * k / n) + Math.Sin(2 * Math.PI * j * k / n)) / scale;
        if (name.Contains("Walsh"))
            return (BitOperations.PopCount((uint)(j & k)) % 2 == 0 ? 1 : -1) / scale;
        if (ComplexOnly(name))
            return Complex.Exp(-Complex.ImaginaryOne * 2 * Math.PI * j * k / n) * (name.Contains("Laplace") ? Math.Exp(-.07f * j) : 1) / scale;
        return j == k ? 1 : j == k - 1 ? -1 : 0;
    }

    public static IEnumerable<object[]> VectorCases()
    {
        foreach (string name in Names)
            foreach (int n in new[]
            {
                1,
                2,
                3,
                4,
                5,
                8,
                17
            }

            )
            {
                if (name.Contains("Walsh") && (n & (n - 1)) != 0)
                    continue;
                foreach (bool normalized in new[]
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
                        if (complex || !ComplexOnly(name))
                            yield return new object[]
                            {
                                name,
                                n,
                                normalized,
                                complex
                            };
            }
    }

    [Theory, MemberData(nameof(VectorCases))]
    public void OneDimensionalTransformsAgreeWithIndependentBasisFunctions(string name, int n, bool normalized, bool complex)
    {
        var x = Enumerable.Range(0, n).Select(i => new Complex32((float)Math.Sin(i * .71) + .3f, complex ? (float)Math.Cos(i * .37) : 0)).ToArray();
        var d = Create(name, normalized: normalized);
        var actual = complex ? d.Forward(x) : d.Forward(x.Select(z => z.Real).ToArray()).Select(v => new Complex32(v, 0)).ToArray();
        Assert.Equal(n, actual.Length);
        for (int k = 0; k < n; k++)
        {
            Complex expected = 0;
            for (int j = 0; j < n; j++)
                expected += Weight(name, k, j, n, normalized) * (Complex)x[j];
            Close(expected, actual[k], 3e-4);
        }

        var restored = complex ? d.Backward(actual) : d.Backward(actual.Select(z => z.Real).ToArray()).Select(v => new Complex32(v, 0)).ToArray();
        double gain = !normalized && Scaled(name) ? n : 1;
        for (int i = 0; i < n; i++)
            Close((Complex)x[i] * gain, restored[i], .001);
    }

    public static IEnumerable<object[]> MatrixCases()
    {
        foreach (string name in Names)
            foreach (var direction in new[]
            {
                Direction.Horizontal,
                Direction.Vertical,
                Direction.Both
            }

            )
                foreach (bool complex in new[]
                {
                    false,
                    true
                }

                )
                    if (complex || !ComplexOnly(name))
                        yield return new object[]
                        {
                            name,
                            direction,
                            complex
                        };
    }

    [Theory, MemberData(nameof(MatrixCases))]
    public void MatrixTransformsUseTheRequestedAxesAndPreserveComplexComponents(string name, Direction direction, bool complex)
    {
        int m = name.Contains("Walsh") ? 4 : 3, n = name.Contains("Walsh") ? 8 : 5;
        var a = new Complex32[m, n];
        var real = new float[m, n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                real[i, j] = (float)Math.Sin(i + j * .3);
                a[i, j] = new(real[i, j], complex ? (float)Math.Cos(i * .7 - j * .4) : 0);
            }

        var d = Create(name, direction);
        var actual = complex ? d.Forward(a) : ToComplex(d.Forward(real));
        Assert.Equal(m, actual.GetLength(0));
        Assert.Equal(n, actual.GetLength(1));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                Complex expected = 0;
                // Complex matrix transforms use the two-sided basis convention U*A*V^H.
                for (int r = 0; r < m; r++)
                    for (int c = 0; c < n; c++)
                        expected += (Complex)a[r, c] * (direction == Direction.Horizontal ? (r == i ? 1 : 0) : Weight(name, i, r, m, true)) * (direction == Direction.Vertical ? (c == j ? 1 : 0) : Complex.Conjugate(Weight(name, j, c, n, true)));
                Close(expected, actual[i, j], .001);
            }

        var restored = complex ? d.Backward(actual) : ToComplex(d.Backward(ToReal(actual)));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Close((Complex)a[i, j], restored[i, j], .002);
    }

    private static float[,] ToReal(Complex32[,] a)
    {
        var r = new float[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < r.GetLength(0); i++)
            for (int j = 0; j < r.GetLength(1); j++)
                r[i, j] = a[i, j].Real;
        return r;
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(8)]
    [InlineData(17)]
    public void HilbertTransformHasTheCorrectPhaseAndRemovesDcAndNyquist(int n)
    {
        var x = Enumerable.Range(0, n).Select(j => new Complex32(2 + (float)Math.Cos(2 * Math.PI * j / n), .25f * (float)Math.Sin(2 * Math.PI * j / n))).ToArray();
        foreach (ITransform d in new ITransform[]
        {
            new HilbertTransform(),
            new FastHilbertTransform()
        }

        )
        {
            var actual = d.Forward(x);
            var restored = d.Backward(actual);
            for (int j = 0; j < n; j++)
            {
                Close(new Complex(Math.Sin(2 * Math.PI * j / n), -.25 * Math.Cos(2 * Math.PI * j / n)), actual[j], 1e-4);
                Close((Complex)x[j] - 2, restored[j], 1e-4);
            }
        }
    }

    [Theory]
    [InlineData(4, 4)]
    [InlineData(8, 12)]
    [InlineData(7, 9)]
    [InlineData(16, 16)]
    public void LaplacianPyramidsReconstructMatricesIncludingOddDimensions(int m, int n)
    {
        var x = NumericAssert.Matrix(m, n);
        var d = new LaplacianPyramidTransform(3, 2);
        Close(x, d.Backward(d.Forward(x)), .001f);
        var z = ToComplex(x);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                z[i, j].Imag = (i - j) * .1f;
        var restored = d.Backward(d.Forward(z));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Close((Complex)z[i, j], restored[i, j], .001);
    }

    [Theory]
    [InlineData(4, false)]
    [InlineData(8, false)]
    [InlineData(9, false)]
    [InlineData(4, true)]
    [InlineData(8, true)]
    [InlineData(9, true)]
    public void LaplacianPyramidsReconstructVectors(int n, bool complex)
    {
        var x = Enumerable.Range(0, n).Select(i => (float)Math.Sin(i * .7)).ToArray();
        var d = new LaplacianPyramidTransform(3, 2);
        if (!complex)
            Close(x, d.Backward(d.Forward(x)), .001f);
        else
        {
            var z = x.Select((v, i) => new Complex32(v, .1f * i)).ToArray();
            var r = d.Backward(d.Forward(z));
            for (int i = 0; i < n; i++)
                Close((Complex)z[i], r[i], .001);
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void GaussianPyramidPreservesConstantSignalsAtEveryScale(int levels)
    {
        var d = new GaussianPyramidTransform(levels, 2);
        var a = new float[16, 24];
        for (int i = 0; i < 16; i++)
            for (int j = 0; j < 24; j++)
                a[i, j] = .375f;
        var x = Enumerable.Repeat(.375f, 32).ToArray();
        var z = x.Select(v => new Complex32(v, .125f)).ToArray();
        var c = ToComplex(a);
        var p = d.Forward(a);
        Assert.Equal(levels, p.Length);
        foreach (var v in p)
            Assert.All(v.Cast<float>(), t => Close(.375, t, 1e-5));
        var q = d.Forward(x);
        foreach (var v in q)
            Assert.All(v, t => Close(.375, t, 1e-5));
        var r = d.Forward(z);
        foreach (var v in r)
            Assert.All(v, t => Close(new Complex(.375, .125), t, 1e-5));
        var s = d.Forward(c);
        foreach (var v in s)
            Assert.All(v.Cast<Complex32>(), t => Close(new Complex(.375, 0), t, 1e-5));
        Assert.Throws<NotSupportedException>(() => d.Backward(p));
        Assert.Throws<NotSupportedException>(() => d.Backward(q));
        Assert.Throws<NotSupportedException>(() => d.Backward(r));
        Assert.Throws<NotSupportedException>(() => d.Backward(s));
    }

    public static IEnumerable<object[]> HankelCases()
    {
        using var stream = typeof(TransformTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.hankel.json");
        using var data = JsonDocument.Parse(stream!);
        foreach (var item in data.RootElement.EnumerateArray())
            yield return new object[]
            {
                item.GetProperty("order").GetInt32(),
                item.GetProperty("size").GetInt32(),
                item.GetProperty("values").GetRawText()
            };
    }

    [Theory]
    [MemberData(nameof(HankelCases))]
    public void HankelMatrixMatchesHighPrecisionBesselZeros(int order, int size, string values)
    {
        using var data = JsonDocument.Parse(values);
        var actual = HankelTransform.Matrix(size, order);
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
                Close(data.RootElement[i][j].GetDouble(), actual[i, j], 5e-4, 5e-4);
        var x = Enumerable.Range(0, size).Select(i => (float)Math.Sin(.7 * i)).ToArray();
        var d = new HankelTransform(order);
        var result = d.Forward(x);
        for (int i = 0; i < size; i++)
        {
            double sum = 0;
            for (int j = 0; j < size; j++)
                sum += data.RootElement[i][j].GetDouble() * x[j];
            Close(sum, result[i], .001, .001);
        }
    }

    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void MultichannelWrappersPreserveChannelSeparationAndInvertTransforms(bool complex, bool matrix)
    {
        var element = matrix ? (complex ? typeof(Complex32[,]) : typeof(float[,])) : (complex ? typeof(Complex32[]) : typeof(float[]));
        var input = Array.CreateInstance(element, 3);
        for (int c = 0; c < 3; c++)
            input.SetValue(MatrixTestSupport.Operand(element, c, 8, 16), c);
        var d = new MultidimensionalTransform(new CosineTransform());
        var p = new MultidimensionalPyramidTransform(new WaveletDecomposition(WaveletPacket.D4, 2));
        foreach (object transform in new object[]
        {
            d,
            p
        }

        )
        {
            var forward = (Array)transform.GetType().GetMethod("Forward", new[] { input.GetType() })!.Invoke(transform, new object[] { input })!;
            var restored = (Array)transform.GetType().GetMethod("Backward", new[] { forward.GetType() })!.Invoke(transform, new object[] { forward })!;
            Assert.Equal(3, restored.Length);
            for (int c = 0; c < 3; c++)
            {
                var a = ((Array)input.GetValue(c)!).Cast<object>().ToArray();
                var b = ((Array)restored.GetValue(c)!).Cast<object>().ToArray();
                Assert.Equal(a.Length, b.Length);
                for (int i = 0; i < a.Length; i++)
                    MatrixTestSupport.Check(MatrixTestSupport.Value(a[i]), b[i], 2e-4);
            }
        }

        var filter = new MultidimensionalFilter(new ThresholdFilter(100));
        filter.GetType().GetMethod("Apply", new[] { input.GetType() })!.Invoke(filter, new object[] { input });
        for (int c = 0; c < 3; c++)
            Assert.All(((Array)input.GetValue(c)!).Cast<object>(), v => MatrixTestSupport.Check(0, v));
    }

    public static IEnumerable<object[]> HankelBoundaryCases()
    {
        using var stream = typeof(TransformTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.meyer-hankel.json");
        using var data = JsonDocument.Parse(stream!);
        foreach (var row in data.RootElement.GetProperty("hankel").EnumerateArray())
            yield return new object[]
            {
                row.GetProperty("order").GetInt32(),
                row.GetProperty("size").GetInt32(),
                row.GetProperty("values").GetRawText()
            };
    }

    [Theory, MemberData(nameof(HankelBoundaryCases))]
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

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(17)]
    [InlineData(32)]
    public void LaplacianPyramidsRecoverComplexVectorsAndThinMatrices(int size)
    {
        foreach (int levels in new[]
        {
            1,
            3,
            10
        }

        )
            foreach (int radius in new[]
            {
                0,
                1,
                2
            }

            )
            {
                var d = new LaplacianPyramidTransform(levels, radius);
                var x = Enumerable.Range(0, size).Select(i => new Complex32((float)Math.Sin(i), (float)Math.Cos(.3 * i))).ToArray();
                var result = d.Backward(d.Forward(x));
                Close(x.Select(v => v.Real).ToArray(), d.Backward(d.Forward(x.Select(v => v.Real).ToArray())), 2e-6);
                for (int i = 0; i < size; i++)
                    Close((Complex)x[i], result[i], 2e-6);
                var thin = new Complex32[1, size];
                for (int i = 0; i < size; i++)
                    thin[0, i] = x[i];
                var restored = d.Backward(d.Forward(thin));
                for (int i = 0; i < size; i++)
                    Close((Complex)x[i], restored[0, i], 2e-6);
            }
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(15)]
    [InlineData(16)]
    [InlineData(17)]
    [InlineData(31)]
    public void FftAgreesWithIndependentDftAndInverts(int n)
    {
        var random = new Random(731911 + n);
        var x = Enumerable.Range(0, n).Select(_ => new Complex32((float)(random.NextDouble() * 2 - 1), (float)(random.NextDouble() * 2 - 1))).ToArray();
        var fft = new FastFourierTransform();
        var spectrum = fft.Forward(x);
        var restored = fft.Backward(spectrum);
        for (int k = 0; k < n; k++)
        {
            Complex sum = 0;
            for (int j = 0; j < n; j++)
                sum += (Complex)x[j] * Complex.Exp(-Complex.ImaginaryOne * 2 * Math.PI * j * k / n);
            Close(sum / Math.Sqrt(n), spectrum[k], 1e-4);
            Close((Complex)x[k], restored[k], 1e-4);
        }
    }

    [Theory]
    [InlineData(Direction.Horizontal)]
    [InlineData(Direction.Vertical)]
    [InlineData(Direction.Both)]
    public void RectangularFftInvertsInEveryDirection(Direction direction)
    {
        var x = new Complex32[3, 5];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 5; j++)
                x[i, j] = new Complex32(i + j * .3f, i * j - 1);
        var fft = new FastFourierTransform(true, direction);
        var actual = fft.Backward(fft.Forward(x));
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 5; j++)
                Close((Complex)x[i, j], actual[i, j], 1e-4);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    [InlineData(8)]
    [InlineData(16)]
    public void OrthogonalRealTransformsInvert(int n)
    {
        var x = Enumerable.Range(0, n).Select(i => (float)Math.Sin(i * .37)).ToArray();
        ITransform[] transforms =
        {
            new CosineTransform(),
            new FastCosineTransform(),
            new SineTransform(),
            new FastSineTransform(),
            new HartleyTransform(),
            new FastHartleyTransform(),
            new WalshHadamardTransform(),
            new FastWalshHadamardTransform(),
            new ChebyshevTransform(),
            new FastChebyshevTransform()
        };
        foreach (var transform in transforms)
            Close(x, transform.Backward(transform.Forward(x)));
    }
}
