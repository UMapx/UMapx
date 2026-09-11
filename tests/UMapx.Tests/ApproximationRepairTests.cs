using System.Drawing;
using System.Drawing.Imaging;
using System.Numerics;
using System.Runtime.Versioning;
using UMapx.Analysis;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.NumericAssert;
using static UMapx.Tests.MatrixAuditTests;

namespace UMapx.Tests;

[Trait("Category", "Analysis")]
public class ApproximationRepairTests
{
    static object Call(string name, params object[] args) => Invoke(typeof(Matrice).GetMethod(name, args.Select(x => x.GetType()).ToArray())!, args);

    public static IEnumerable<object[]> PadeCases()
    {
        foreach (bool complex in new[] { false, true })
            for (int m = 1; m <= 4; m++) for (int n = 1; n <= 5; n++) yield return new object[] { complex, m, n };
    }

    [Theory]
    [MemberData(nameof(PadeCases))]
    public void PadeReconstructsIndependentTaylorProductEquations(bool complex, int m, int n)
    {
        // Taylor coefficients of exp(s*x), with s chosen independently of the solver.
        Complex s = complex ? new Complex(.6, .4) : .8;
        var coefficients = new Complex32[m + n + 1]; coefficients[0] = 1;
        Complex term = 1;
        for (int k = 1; k < coefficients.Length; k++) { term *= s / k; coefficients[k] = (Complex32)term; }
        var pade = new Pade(m, n);
        Complex32[] p, q;
        if (complex) (p, q) = pade.Compute(coefficients);
        else
        {
            var result = pade.Compute(coefficients.Select(z => z.Real).ToArray());
            p = result.NumeratorCoeffs.Select(x => new Complex32(x, 0)).ToArray();
            q = result.DenominatorCoeffs.Select(x => new Complex32(x, 0)).ToArray();
        }
        Assert.Equal(m + 1, p.Length); Assert.Equal(n + 1, q.Length); Assert.Equal(1, q[0].Real);
        for (int k = 0; k <= m + n; k++)
        {
            Complex product = 0;
            for (int j = 0; j <= Math.Min(k, n); j++) product += (Complex)q[j] * (Complex)coefficients[k - j];
            Close(k <= m ? (Complex)p[k] : Complex.Zero, (Complex32)product, 2e-6, 2e-5);
        }
        // Closed-form coefficient ratios follow from matching the Taylor product, not Matrice.Solve.
        double Factorial(int k) { double value = 1; for (int j = 2; j <= k; j++) value *= j; return value; }
        for (int j = 0; j <= n; j++)
        {
            Complex expected = Complex.Pow(-s, j) * Factorial(m + n - j) * Factorial(n) / (Factorial(m + n) * Factorial(n - j) * Factorial(j));
            Close(expected, q[j], 3e-5, 3e-4);
        }
    }

    public static IEnumerable<object[]> GridCases()
    {
        foreach (float x in new[] { -20f, -3f, -1f, 0f, 2f, 4f, 7f, 9f, 20f })
            foreach (float y in new[] { -20f, -4f, -2f, 1f, 2f, 3f, 8f, 20f }) yield return new object[] { x, y };
    }

    [Theory]
    [MemberData(nameof(GridCases))]
    public void BilinearGridReproducesMixedLinearSurfacesWithCoordinateClamping(float x, float y)
    {
        float[] xs = { -3, 0, 4, 9 }, ys = { -4, 1, 3, 8 };
        double Surface(double a, double b) => 2 + 3 * a - .5 * b + .25 * a * b;
        var z = new float[xs.Length, ys.Length];
        for (int i = 0; i < xs.Length; i++) for (int j = 0; j < ys.Length; j++) z[i, j] = (float)Surface(xs[i], ys[j]);
        Close(Surface(Math.Clamp(x, -3, 9), Math.Clamp(y, -4, 8)), new Interpolation().Compute(xs, ys, z, x, y));
    }

    [Fact]
    public void BilinearDifferencesDoNotOverflowBeforeInterpolation()
    {
        float[] x = { -3e38f, 3e38f }, y = { 0, 1 };
        float[,] z = { { -3e38f, -3e38f }, { 3e38f, 3e38f } };
        Assert.Equal(0, new Interpolation().Compute(x, y, z, 0, .5f));
    }

    static Complex[] ReferenceMean(Complex[] values, Complex[]? weights, int window)
    {
        if (values.Length < 2 || window < 2) return values.ToArray();
        var result = new Complex[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            Complex sum = 0, total = 0;
            for (int j = 0; j < values.Length; j++)
                if (j >= (long)i - window / 2 && j <= (long)i + (window - 1) / 2)
                { Complex w = weights?[j] ?? Complex.One; sum += values[j] * w; total += w; }
            result[i] = total == 0 ? 0 : sum / total;
        }
        return result;
    }

    public static IEnumerable<object[]> MeanCases()
    {
        foreach (bool complex in new[] { false, true }) foreach (bool weighted in new[] { false, true })
            foreach (int n in new[] { 0, 1, 2, 3, 9 }) foreach (int window in new[] { 1, 2, 3, 4, 9, 20, int.MaxValue })
                yield return new object[] { complex, weighted, n, window };
    }

    [Theory]
    [MemberData(nameof(MeanCases))]
    public void LocalMeansMatchDirectWindowSumsAtEveryBoundary(bool complex, bool weighted, int n, int window)
    {
        var values = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, n);
        var weights = (Array)Operand(values.GetType(), 3, 1, n);
        // Tiny weights expose the former artificial epsilon in the denominator.
        for (int i = 0; i < n; i++) weights.SetValue(complex ? (object)(Complex32)(Value(weights, 0, i) * 1e-20) : (float)(Value(weights, 0, i).Real * 1e-20), i);
        var expected = ReferenceMean(values.Cast<object>().Select(v => Value(v)).ToArray(), weighted ? weights.Cast<object>().Select(v => Value(v)).ToArray() : null, window);
        var actual = (Array)(weighted ? Call("Mean", values, weights, window) : Call("Mean", values, window));
        Assert.Equal(n, actual.Length);
        for (int i = 0; i < n; i++) Check(expected[i], actual.GetValue(i)!);
    }

    public static IEnumerable<object[]> MatrixMeanCases()
    {
        foreach (bool complex in new[] { false, true }) foreach (bool weighted in new[] { false, true })
            foreach (var shape in new[] { (0, 3), (3, 0), (1, 5), (5, 1), (3, 7), (7, 3) })
                foreach (var window in new[] { (1, 1), (3, 5), (4, 2), (20, 20) })
                    yield return new object[] { complex, weighted, shape.Item1, shape.Item2, window.Item1, window.Item2 };
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
            for (int j = 0; j < cols; j++) expected[i, j] = line[j];
        }
        for (int j = 0; j < cols; j++)
        {
            var line = ReferenceMean(Enumerable.Range(0, rows).Select(i => expected[i, j]).ToArray(), weighted ? Enumerable.Range(0, rows).Select(i => Value(weights, i, j)).ToArray() : null, heightWindow);
            for (int i = 0; i < rows; i++) expected[i, j] = line[i];
        }
        var result = (Array)(weighted ? Call("Mean", a, weights, heightWindow, widthWindow) : Call("Mean", a, heightWindow, widthWindow));
        Assert.Equal(rows, result.GetLength(0)); Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) Check(expected[i, j], result.GetValue(i, j)!);
    }

    public static IEnumerable<object[]> MorphCases()
    {
        foreach (var mode in Enum.GetValues<MorphologyMode>())
            foreach (var shape in new[] { (1, 1), (1, 7), (7, 1), (3, 5), (5, 3) })
                foreach (var window in new[] { (1, 1), (3, 3), (5, 7) })
                    yield return new object[] { mode, shape.Item1, shape.Item2, window.Item1, window.Item2 };
    }

    [Theory]
    [MemberData(nameof(MorphCases))]
    public void MorphologyMatchesSortedWindowsIncludingDuplicateValues(MorphologyMode mode, int rows, int cols, int heightWindow, int widthWindow)
    {
        var a = new float[rows, cols];
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) a[i, j] = (i * 17 + j * 7) % 5 - 2;
        var result = a.Morph(heightWindow, widthWindow, mode);
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++)
        {
            var samples = new List<float>();
            for (int di = -heightWindow / 2; di <= heightWindow / 2; di++) for (int dj = -widthWindow / 2; dj <= widthWindow / 2; dj++) samples.Add(a[Math.Clamp(i + di, 0, rows - 1), Math.Clamp(j + dj, 0, cols - 1)]);
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

    // Independent cubic Hermite polynomial with centered endpoint slopes.
    static Complex Cubic(Complex a, Complex b, Complex c, Complex d, double t) =>
        ((.5 * (-a + 3 * b - 3 * c + d) * t + a - 2.5 * b + 2 * c - .5 * d) * t + .5 * (c - a)) * t + b;

    static Complex Resample(Complex[,] input, double y, double x)
    {
        int iy = (int)Math.Floor(y), ix = (int)Math.Floor(x);
        Complex At(int i, int j) => input[Math.Clamp(i, 0, input.GetLength(0) - 1), Math.Clamp(j, 0, input.GetLength(1) - 1)];
        Complex Row(int i) => Cubic(At(i, ix - 1), At(i, ix), At(i, ix + 1), At(i, ix + 2), x - ix);
        return Cubic(Row(iy - 1), Row(iy), Row(iy + 1), Row(iy + 2), y - iy);
    }

    public static IEnumerable<object[]> ResizeCases()
    {
        foreach (bool complex in new[] { false, true })
            foreach (var shape in new[] { (1, 1), (1, 5), (5, 1), (3, 5), (6, 4) })
                foreach (var target in new[] { (1, 1), (3, 5), (7, 9) }) yield return new object[] { complex, shape.Item1, shape.Item2, target.Item1, target.Item2 };
    }

    [Theory]
    [MemberData(nameof(ResizeCases))]
    public void BicubicResizeMatchesIndependentHermiteInterpolation(bool complex, int rows, int cols, int targetRows, int targetCols)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 3, rows, cols);
        var reference = new Complex[rows, cols];
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) reference[i, j] = Value(a, i, j);
        var actual = (Array)Call("Resize", a, targetRows, targetCols, InterpolationMode.Bicubic);
        for (int i = 0; i < targetRows; i++) for (int j = 0; j < targetCols; j++)
            Check(Resample(reference, (i + .5) * rows / targetRows - .5, (j + .5) * cols / targetCols - .5), actual.GetValue(i, j)!);
        var vector = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 3, 1, cols);
        var line = new Complex[1, cols]; for (int j = 0; j < cols; j++) line[0, j] = Value(vector, 0, j);
        var resized = (Array)Call("Resize", vector, targetCols, InterpolationMode.Bicubic);
        for (int j = 0; j < targetCols; j++) Check(Resample(line, 0, (j + .5) * cols / targetCols - .5), resized.GetValue(j)!);
    }

    public static IEnumerable<object[]> RotationCases()
    {
        foreach (bool complex in new[] { false, true }) foreach (var mode in Enum.GetValues<InterpolationMode>())
            foreach (int n in new[] { 1, 2, 5, 6 }) foreach (float angle in new[] { -540f, -90f, 0f, 90f, 180f, 270f, 360f }) yield return new object[] { complex, mode, n, angle };
    }

    [Theory]
    [MemberData(nameof(RotationCases))]
    public void OrthogonalRotationsMatchExactIndexPermutations(bool complex, InterpolationMode mode, int n, float angle)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, n, n);
        var result = (Array)Call("Rotate", a, angle, mode);
        int turn = ((int)angle % 360 + 360) % 360;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        {
            var p = turn switch { 90 => (j, n - 1 - i), 180 => (n - 1 - i, n - 1 - j), 270 => (n - 1 - j, i), _ => (i, j) };
            Check(Value(a, p.Item1, p.Item2), result.GetValue(i, j)!);
        }
    }

    [Theory]
    [InlineData(1, 1)] [InlineData(3, 5)] [InlineData(7, 9)]
    [SupportedOSPlatform("windows")]
    public void BitmapBicubicResizeMatchesIndependentChannelInterpolation(int height, int width)
    {
        using var input = new Bitmap(5, 3, PixelFormat.Format32bppArgb);
        using var output = new Bitmap(width, height, PixelFormat.Format32bppArgb);
        var planes = Enumerable.Range(0, 4).Select(_ => new Complex[3, 5]).ToArray();
        for (int i = 0; i < 3; i++) for (int j = 0; j < 5; j++)
        {
            var color = Color.FromArgb(70 + i * 20 + j * 7, i * 45 + j * 17, 180 - i * 12 - j * 21, 20 + i * 25 + j * 31);
            input.SetPixel(j, i, color);
            int[] channels = { color.A, color.R, color.G, color.B };
            for (int k = 0; k < 4; k++) planes[k][i, j] = channels[k];
        }
        new UMapx.Imaging.Resize(width, height, InterpolationMode.Bicubic).Apply(output, input);
        for (int i = 0; i < height; i++) for (int j = 0; j < width; j++)
        {
            var c = output.GetPixel(j, i); int[] channels = { c.A, c.R, c.G, c.B };
            for (int k = 0; k < 4; k++)
            {
                double expected = Math.Clamp(Resample(planes[k], (i + .5) * 3 / height - .5, (j + .5) * 5 / width - .5).Real, 0, 255);
                Close(expected, channels[k], 1.001, 0);
            }
        }
    }
}
