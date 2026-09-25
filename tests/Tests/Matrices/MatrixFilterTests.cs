using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixFilterTests
{
    [Theory]
    [InlineData(false, false)]
    [InlineData(false, true)]
    [InlineData(true, false)]
    [InlineData(true, true)]
    public void FilteringMatchesCenteredCorrelationWithTruncatedSupport(bool complex, bool normalize)
    {
        // Conv uses correlation orientation; this audit records that API convention.
        var v = Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, 9);
        var k = new float[]
        {
            1,
            2,
            4
        };
        var result = (Array)Invoke(typeof(Matrice).GetMethod("Conv", new[] { v.GetType(), typeof(float[]), typeof(bool) })!, v, k, normalize);
        for (int i = 0; i < 9; i++)
        {
            Complex sum = 0;
            double weight = 0;
            for (int j = 0; j < 3; j++)
                if (i + j - 1 >= 0 && i + j - 1 < 9)
                {
                    sum += Value(v, 0, i + j - 1) * k[j];
                    weight += k[j];
                }

            Check(normalize ? sum / weight : sum, result.GetValue(i)!);
        }
    }

    [Theory]
    [InlineData(MorphologyMode.Median)]
    [InlineData(MorphologyMode.Erosion)]
    [InlineData(MorphologyMode.Dilatation)]
    public void MorphologyMatchesSortingOfReplicatedEdgeNeighborhoods(MorphologyMode mode)
    {
        // Public Matrice arguments are window sizes; the internal filter takes half sizes.
        var a = Matrix(5, 7);
        var actual = a.Morph(3, 5, mode);
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 7; j++)
            {
                var v = new List<float>();
                for (int di = -1; di <= 1; di++)
                    for (int dj = -2; dj <= 2; dj++)
                        v.Add(a[Math.Clamp(i + di, 0, 4), Math.Clamp(j + dj, 0, 6)]);
                v.Sort();
                Close(mode == MorphologyMode.Erosion ? v[0] : mode == MorphologyMode.Dilatation ? v[^1] : v[v.Count / 2], actual[i, j]);
            }

        var row = Enumerable.Range(0, 7).Select(j => a[2, j]).ToArray();
        var filtered = row.Morph(5, mode);
        for (int j = 0; j < 7; j++)
        {
            var v = Enumerable.Range(-2, 5).Select(d => row[Math.Clamp(j + d, 0, 6)]).OrderBy(x => x).ToArray();
            Close(mode == MorphologyMode.Erosion ? v[0] : mode == MorphologyMode.Dilatation ? v[^1] : v[2], filtered[j]);
        }
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    public void LocalMeanHasTheCorrectInteriorImpulseResponse(bool complex)
    {
        var x = new float[]
        {
            0,
            0,
            0,
            3,
            0,
            0,
            0,
            0,
            0
        };
        var v = complex ? (object)x.ToComplex() : x;
        var result = (Array)Invoke(typeof(Matrice).GetMethod("Mean", new[] { v.GetType(), typeof(int) })!, v, 3);
        for (int i = 1; i < x.Length - 1; i++)
            Check((x[i - 1] + x[i] + x[i + 1]) / 3.0, result.GetValue(i)!);
    }

    [Theory]
    [InlineData(MorphologyMode.Median)]
    [InlineData(MorphologyMode.Dilatation)]
    public void MorphologyUsesZeroBasedRanksInAThreeSampleWindow(MorphologyMode mode)
    {
        float[] input =
        {
            1,
            2,
            3,
            4,
            5
        };
        var actual = input.Morph(3, mode);
        Close(mode == MorphologyMode.Median ? 3 : 4, actual[2]);
    }

    public static IEnumerable<object[]> ConvolutionCases()
    {
        foreach (bool a in new[]
        {
            false,
            true
        }

        )
            foreach (bool b in new[]
            {
                false,
                true
            }

            )
                foreach (bool normalized in new[]
                {
                    false,
                    true
                }

                )
                    foreach (string direction in new[]
                    {
                        "Matrix",
                        "Horizontal",
                        "Vertical",
                        "Both"
                    }

                    )
                        yield return new object[]
                        {
                            a,
                            b,
                            normalized,
                            direction
                        };
    }

    static Complex[,] Correlate(Complex[,] a, Complex[,] kernel, bool normalized)
    {
        int h = a.GetLength(0), w = a.GetLength(1), kh = kernel.GetLength(0), kw = kernel.GetLength(1);
        var r = new Complex[h, w];
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                Complex sum = 0, weight = 0;
                for (int dy = 0; dy < kh; dy++)
                    for (int dx = 0; dx < kw; dx++)
                    {
                        int yy = y + dy - kh / 2, xx = x + dx - kw / 2;
                        if (yy < 0 || yy >= h || xx < 0 || xx >= w)
                            continue;
                        sum += a[yy, xx] * kernel[dy, dx];
                        weight += kernel[dy, dx];
                    }

                r[y, x] = normalized ? sum / weight : sum;
            }

        return r;
    }

    [Theory]
    [MemberData(nameof(ConvolutionCases))]
    public void MixedMatrixConvolutionOverloadsMatchIndependentNeighborhoodSums(bool complexData, bool complexKernel, bool normalized, string direction)
    {
        var a = (Array)MatrixTestSupport.Operand(complexData ? typeof(Complex32[,]) : typeof(float[,]), 1, 7, 9);
        var k = (Array)MatrixTestSupport.Operand(direction == "Matrix" ? (complexKernel ? typeof(Complex32[,]) : typeof(float[,])) : (complexKernel ? typeof(Complex32[]) : typeof(float[])), 2, 3, 3);
        object[] args = direction == "Matrix" ? new object[]
        {
            a,
            k,
            normalized
        }

        : new object[]
        {
            a,
            k,
            Enum.Parse<Direction>(direction),
            normalized
        };
        var actual = (Array)MatrixTestSupport.Invoke(typeof(Matrice).GetMethod("Conv", args.Select(v => v.GetType()).ToArray())!, args);
        var expected = new Complex[7, 9];
        for (int y = 0; y < 7; y++)
            for (int x = 0; x < 9; x++)
                expected[y, x] = MatrixTestSupport.Value(a, y, x);
        if (direction == "Matrix")
        {
            var kernel = new Complex[3, 3];
            for (int y = 0; y < 3; y++)
                for (int x = 0; x < 3; x++)
                    kernel[y, x] = MatrixTestSupport.Value(k, y, x);
            expected = Correlate(expected, kernel, normalized);
        }
        else
        {
            foreach (bool vertical in direction == "Both" ? new[]
            {
                false,
                true
            }

            : new[]
            {
                direction == "Vertical"
            }

            )
            {
                var kernel = new Complex[vertical ? 3 : 1, vertical ? 1 : 3];
                for (int i = 0; i < 3; i++)
                    kernel[vertical ? i : 0, vertical ? 0 : i] = MatrixTestSupport.Value(k, 0, i);
                expected = Correlate(expected, kernel, normalized);
            }
        }

        for (int y = 0; y < 7; y++)
            for (int x = 0; x < 9; x++)
                MatrixTestSupport.Check(expected[y, x], actual.GetValue(y, x)!, 3e-4);
    }

    static Complex[] ReferenceMean(Complex[] values, Complex[]? weights, int window)
    {
        if (values.Length < 2 || window < 2)
            return values.ToArray();
        var result = new Complex[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            Complex sum = 0, total = 0;
            for (int j = 0; j < values.Length; j++)
                if (j >= (long)i - window / 2 && j <= (long)i + (window - 1) / 2)
                {
                    Complex w = weights?[j] ?? Complex.One;
                    sum += values[j] * w;
                    total += w;
                }

            result[i] = total == 0 ? 0 : sum / total;
        }

        return result;
    }

    public static IEnumerable<object[]> MeanCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (bool weighted in new[]
            {
                false,
                true
            }

            )
                foreach (int n in new[]
                {
                    0,
                    1,
                    2,
                    3,
                    9
                }

                )
                    foreach (int window in new[]
                    {
                        1,
                        2,
                        3,
                        4,
                        9,
                        20,
                        int.MaxValue
                    }

                    )
                        yield return new object[]
                        {
                            complex,
                            weighted,
                            n,
                            window
                        };
    }

    [Theory]
    [MemberData(nameof(MeanCases))]
    public void LocalMeansMatchDirectWindowSumsAtEveryBoundary(bool complex, bool weighted, int n, int window)
    {
        var values = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, n);
        var weights = (Array)Operand(values.GetType(), 3, 1, n);
        // Tiny weights expose the former artificial epsilon in the denominator.
        for (int i = 0; i < n; i++)
            weights.SetValue(complex ? (object)(Complex32)(Value(weights, 0, i) * 1e-20) : (float)(Value(weights, 0, i).Real * 1e-20), i);
        var expected = ReferenceMean(values.Cast<object>().Select(v => Value(v)).ToArray(), weighted ? weights.Cast<object>().Select(v => Value(v)).ToArray() : null, window);
        var actual = (Array)(weighted ? Call("Mean", values, weights, window) : Call("Mean", values, window));
        Assert.Equal(n, actual.Length);
        for (int i = 0; i < n; i++)
            Check(expected[i], actual.GetValue(i)!);
    }

    public static IEnumerable<object[]> MatrixMeanCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (bool weighted in new[]
            {
                false,
                true
            }

            )
                foreach (var shape in new[]
                {
                    (0, 3),
                    (3, 0),
                    (1, 5),
                    (5, 1),
                    (3, 7),
                    (7, 3)
                }

                )
                    foreach (var window in new[]
                    {
                        (1, 1),
                        (3, 5),
                        (4, 2),
                        (20, 20)
                    }

                    )
                        yield return new object[]
                        {
                            complex,
                            weighted,
                            shape.Item1,
                            shape.Item2,
                            window.Item1,
                            window.Item2
                        };
    }

    [Theory]
    [MemberData(nameof(MatrixMeanCases))]
    public void MatrixMeansPreserveTheSeparableHorizontalThenVerticalContract(bool complex, bool weighted, int rows, int cols, int heightWindow, int widthWindow)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var weights = (Array)Operand(a.GetType(), 3, rows, cols);
        var expected = new Complex[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            var line = ReferenceMean(Enumerable.Range(0, cols).Select(j => Value(a, i, j)).ToArray(), weighted ? Enumerable.Range(0, cols).Select(j => Value(weights, i, j)).ToArray() : null, widthWindow);
            for (int j = 0; j < cols; j++)
                expected[i, j] = line[j];
        }

        for (int j = 0; j < cols; j++)
        {
            var line = ReferenceMean(Enumerable.Range(0, rows).Select(i => expected[i, j]).ToArray(), weighted ? Enumerable.Range(0, rows).Select(i => Value(weights, i, j)).ToArray() : null, heightWindow);
            for (int i = 0; i < rows; i++)
                expected[i, j] = line[i];
        }

        var result = (Array)(weighted ? Call("Mean", a, weights, heightWindow, widthWindow) : Call("Mean", a, heightWindow, widthWindow));
        Assert.Equal(rows, result.GetLength(0));
        Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                Check(expected[i, j], result.GetValue(i, j)!);
    }

    public static IEnumerable<object[]> MorphCases()
    {
        foreach (var mode in Enum.GetValues<MorphologyMode>())
            foreach (var shape in new[]
            {
                (1, 1),
                (1, 7),
                (7, 1),
                (3, 5),
                (5, 3)
            }

            )
                foreach (var window in new[]
                {
                    (1, 1),
                    (3, 3),
                    (5, 7)
                }

                )
                    yield return new object[]
                    {
                        mode,
                        shape.Item1,
                        shape.Item2,
                        window.Item1,
                        window.Item2
                    };
    }

    [Theory]
    [MemberData(nameof(MorphCases))]
    public void MorphologyMatchesSortedWindowsIncludingDuplicateValues(MorphologyMode mode, int rows, int cols, int heightWindow, int widthWindow)
    {
        var a = new float[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                a[i, j] = (i * 17 + j * 7) % 5 - 2;
        var result = a.Morph(heightWindow, widthWindow, mode);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                var samples = new List<float>();
                for (int di = -heightWindow / 2; di <= heightWindow / 2; di++)
                    for (int dj = -widthWindow / 2; dj <= widthWindow / 2; dj++)
                        samples.Add(a[Math.Clamp(i + di, 0, rows - 1), Math.Clamp(j + dj, 0, cols - 1)]);
                samples.Sort();
                Assert.Equal(samples[mode == MorphologyMode.Erosion ? 0 : mode == MorphologyMode.Dilatation ? samples.Count - 1 : samples.Count / 2], result[i, j]);
            }

        var vector = Enumerable.Range(0, cols).Select(j => a[0, j]).ToArray();
        var filtered = vector.Morph(widthWindow, mode);
        for (int j = 0; j < cols; j++)
        {
            var samples = Enumerable.Range(-widthWindow / 2, widthWindow).Select(k => vector[Math.Clamp(j + k, 0, cols - 1)]).OrderBy(x => x).ToArray();
            Assert.Equal(samples[mode == MorphologyMode.Erosion ? 0 : mode == MorphologyMode.Dilatation ? samples.Length - 1 : samples.Length / 2], filtered[j]);
        }
    }
}
