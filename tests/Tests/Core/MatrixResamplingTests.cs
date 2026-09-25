using System.Numerics;
using UMapx.Core;
using Xunit;
using static UMapx.Tests.MatrixTestSupport;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixResamplingTests
{
    // Independent cubic Hermite polynomial with centered endpoint slopes.
    static Complex Cubic(Complex a, Complex b, Complex c, Complex d, double t) => ((.5 * (-a + 3 * b - 3 * c + d) * t + a - 2.5 * b + 2 * c - .5 * d) * t + .5 * (c - a)) * t + b;
    static Complex Resample(Complex[,] input, double y, double x)
    {
        int iy = (int)Math.Floor(y), ix = (int)Math.Floor(x);
        Complex At(int i, int j) => input[Math.Clamp(i, 0, input.GetLength(0) - 1), Math.Clamp(j, 0, input.GetLength(1) - 1)];
        Complex Row(int i) => Cubic(At(i, ix - 1), At(i, ix), At(i, ix + 1), At(i, ix + 2), x - ix);
        return Cubic(Row(iy - 1), Row(iy), Row(iy + 1), Row(iy + 2), y - iy);
    }

    public static IEnumerable<object[]> ResizeCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var shape in new[]
            {
                (1, 1),
                (1, 5),
                (5, 1),
                (3, 5),
                (6, 4)
            }

            )
                foreach (var target in new[]
                {
                    (1, 1),
                    (3, 5),
                    (7, 9)
                }

                )
                    yield return new object[]
                    {
                        complex,
                        shape.Item1,
                        shape.Item2,
                        target.Item1,
                        target.Item2
                    };
    }

    [Theory]
    [MemberData(nameof(ResizeCases))]
    public void BicubicResizeMatchesIndependentHermiteInterpolation(bool complex, int rows, int cols, int targetRows, int targetCols)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 3, rows, cols);
        var reference = new Complex[rows, cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                reference[i, j] = Value(a, i, j);
        var actual = (Array)Call("Resize", a, targetRows, targetCols, InterpolationMode.Bicubic);
        for (int i = 0; i < targetRows; i++)
            for (int j = 0; j < targetCols; j++)
                Check(Resample(reference, (i + .5) * rows / targetRows - .5, (j + .5) * cols / targetCols - .5), actual.GetValue(i, j)!);
        var vector = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 3, 1, cols);
        var line = new Complex[1, cols];
        for (int j = 0; j < cols; j++)
            line[0, j] = Value(vector, 0, j);
        var resized = (Array)Call("Resize", vector, targetCols, InterpolationMode.Bicubic);
        for (int j = 0; j < targetCols; j++)
            Check(Resample(line, 0, (j + .5) * cols / targetCols - .5), resized.GetValue(j)!);
    }

    public static IEnumerable<object[]> RotationCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var mode in Enum.GetValues<InterpolationMode>())
                foreach (int n in new[]
                {
                    1,
                    2,
                    5,
                    6
                }

                )
                    foreach (float angle in new[]
                    {
                        -540f,
                        -90f,
                        0f,
                        90f,
                        180f,
                        270f,
                        360f
                    }

                    )
                        yield return new object[]
                        {
                            complex,
                            mode,
                            n,
                            angle
                        };
    }

    [Theory]
    [MemberData(nameof(RotationCases))]
    public void OrthogonalRotationsMatchExactIndexPermutations(bool complex, InterpolationMode mode, int n, float angle)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, n, n);
        var result = (Array)Call("Rotate", a, angle, mode);
        int turn = ((int)angle % 360 + 360) % 360;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                var p = turn switch
                {
                    90 => (j, n - 1 - i),
                    180 => (n - 1 - i, n - 1 - j),
                    270 => (n - 1 - j, i),
                    _ => (i, j)
                };
                Check(Value(a, p.Item1, p.Item2), result.GetValue(i, j)!);
            }
    }

    public static IEnumerable<object[]> HalfTurnCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var mode in Enum.GetValues<InterpolationMode>())
                foreach (float angle in new[]
                {
                    0f,
                    180f
                }

                )
                    yield return new object[]
                    {
                        complex,
                        mode,
                        angle
                    };
    }

    [Theory]
    [MemberData(nameof(HalfTurnCases))]
    public void MatrixRotationsAtExactHalfTurnsMatchIndexReversal(bool complex, InterpolationMode mode, float angle)
    {
        var a = (Array)MatrixTestSupport.Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, 5, 7);
        var result = (Array)MatrixTestSupport.Invoke(typeof(Matrice).GetMethod("Rotate", new[] { a.GetType(), typeof(float), typeof(InterpolationMode) })!, a, angle, mode);
        for (int y = 0; y < 5; y++)
            for (int x = 0; x < 7; x++)
                MatrixTestSupport.Check(MatrixTestSupport.Value(a, angle == 0 ? y : 4 - y, angle == 0 ? x : 6 - x), result.GetValue(y, x)!, 2e-4);
    }

    public static IEnumerable<object[]> IdentityResizeCases()
    {
        foreach (bool complex in new[]
        {
            false,
            true
        }

        )
            foreach (var mode in Enum.GetValues<InterpolationMode>())
                yield return new object[]
                {
                    complex,
                    mode
                };
    }

    [Theory]
    [MemberData(nameof(IdentityResizeCases))]
    public void ResizingAnArrayToItsCurrentSizePreservesSamples(bool complex, InterpolationMode mode)
    {
        var a = (Array)MatrixTestSupport.Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, 3, 5);
        var result = (Array)MatrixTestSupport.Invoke(typeof(Matrice).GetMethod("Resize", new[] { a.GetType(), typeof(int), typeof(int), typeof(InterpolationMode) })!, a, 3, 5, mode);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 5; j++)
                MatrixTestSupport.Check(MatrixTestSupport.Value(a, i, j), result.GetValue(i, j)!);
    }
}
