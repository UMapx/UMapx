using System.Numerics;
using UMapx.Core;
using UMapx.Distance;
using Xunit;
using static UMapx.Tests.MatrixAuditTests;
using static UMapx.Tests.NumericAssert;

namespace UMapx.Tests;

[Trait("Category", "Matrix")]
public class MatrixRepairTests
{
    static object Call(string name, params object[] args) => Invoke(typeof(Matrice).GetMethod(name, args.Select(x => x.GetType()).ToArray())!, args);
    static int Wrap(long index, int length) => (int)((index % length + length) % length);

    public static IEnumerable<object[]> DiagonalCases()
    {
        foreach (var flags in MatrixAuditTests.DiagonalCases())
            foreach (var shape in new[] { (0, 3), (3, 0), (1, 7), (7, 1), (2, 4), (5, 3) })
                yield return flags.Concat(new object[] { shape.Item1, shape.Item2 }).ToArray();
    }

    [Theory]
    [MemberData(nameof(DiagonalCases))]
    public void AllDiagonalProductsRespectRectangularDimensionsAndZeroPseudoinverses(bool left, bool inverse, bool complexMatrix, bool complexVector, int rows, int cols)
    {
        var a = (Array)Operand(complexMatrix ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var v = (Array)Operand(complexVector ? typeof(Complex32[]) : typeof(float[]), 2, 1, left ? rows : cols);
        if (v.Length > 0) v.SetValue(complexVector ? (object)new Complex32(0, 0) : 0f, 0);
        var result = (Array)(left ? Call("Dot", v, a, inverse) : Call("Dot", a, v, inverse));
        Assert.Equal(rows, result.GetLength(0)); Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++)
        {
            Complex diagonal = Value(v, 0, left ? i : j);
            Complex expected = diagonal == 0 ? 0 : inverse ? Value(a, i, j) / diagonal : Value(a, i, j) * diagonal;
            Check(expected, result.GetValue(i, j)!);
        }
    }

    public static IEnumerable<object[]> ShiftCases()
    {
        foreach (bool complex in new[] { false, true })
            foreach (var shape in new[] { (0, 3), (3, 0), (1, 7), (7, 1), (2, 5), (5, 3) })
                foreach (int shift in new[] { int.MinValue, -123456789, -1, 0, 1, int.MaxValue })
                    yield return new object[] { complex, shape.Item1, shape.Item2, shift };
    }

    [Theory]
    [MemberData(nameof(ShiftCases))]
    public void ShiftsWrapBothAxesWithoutIntegerOverflow(bool complex, int rows, int cols, int shift)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, rows, cols);
        var result = (Array)Call("Shift", a, shift, shift);
        Assert.Equal(rows, result.GetLength(0)); Assert.Equal(cols, result.GetLength(1));
        for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++)
            Check(Value(a, Wrap((long)i - shift, rows), Wrap((long)j - shift, cols)), result.GetValue(i, j)!);
        var v = (Array)Operand(complex ? typeof(Complex32[]) : typeof(float[]), 1, 1, cols);
        var shifted = (Array)Call("Shift", v, shift);
        for (int j = 0; j < cols; j++) Check(Value(v, 0, Wrap((long)j - shift, cols)), shifted.GetValue(j)!);
    }

    public static IEnumerable<object[]> MergeCases()
    {
        foreach (int kind in new[] { 0, 1, 2 })
            foreach (bool vector in new[] { false, true })
                foreach (var offset in new[] { (-2, -1), (0, 0), (1, 2), (3, 6), (4, 7), (-10, 2), (int.MinValue, int.MinValue), (int.MaxValue, int.MaxValue) })
                    yield return new object[] { kind, vector, offset.Item1, offset.Item2 };
    }

    [Theory]
    [MemberData(nameof(MergeCases))]
    public void MergeClipsThePatchIntersectionAndPreservesBothInputs(int kind, bool vector, int top, int left)
    {
        Type Element(bool complex) => vector ? complex ? typeof(Complex32[]) : typeof(float[]) : complex ? typeof(Complex32[,]) : typeof(float[,]);
        var a = (Array)Operand(Element(kind != 0), 1, 4, 7);
        var b = (Array)Operand(Element(kind == 2), 8, 3, 4);
        var original = (Array)a.Clone(); var patch = (Array)b.Clone();
        var result = (Array)(vector ? Call("Merge", a, b, left, 4) : Call("Merge", a, b, top, left, 3, 4));
        Assert.NotSame(a, result);
        for (int i = 0; i < (vector ? 1 : 4); i++) for (int j = 0; j < 7; j++)
        {
            long pi = vector ? 0 : (long)i - top, pj = (long)j - left;
            Complex expected = pi >= 0 && pi < (vector ? 1 : 3) && pj >= 0 && pj < 4 ? Value(patch, (int)pi, (int)pj) : Value(original, i, j);
            Check(expected, vector ? result.GetValue(j)! : result.GetValue(i, j)!);
            Check(Value(original, i, j), vector ? a.GetValue(j)! : a.GetValue(i, j)!);
        }
        Assert.Equal(patch.Cast<object>(), b.Cast<object>());
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void EmptyPatchesAreIdentityInsertions(bool complex)
    {
        var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1);
        foreach (var size in new[] { (0, 2), (2, 0), (0, 0) })
        {
            var b = Array.CreateInstance(a.GetType().GetElementType()!, size.Item1, size.Item2);
            var result = (Array)Call("Merge", a, b, 1, 1, size.Item1, size.Item2);
            Assert.Equal(a.Cast<object>(), result.Cast<object>());
        }
    }

    [Theory]
    [InlineData(false, false)] [InlineData(false, true)] [InlineData(true, false)] [InlineData(true, true)]
    public void NonseparableMeshesEvaluateEveryCartesianPair(bool complexX, bool complexY)
    {
        foreach (var size in new[] { (0, 4), (3, 0), (1, 7), (5, 3) })
        {
            var x = (Array)Operand(complexX ? typeof(Complex32[]) : typeof(float[]), 2, 1, size.Item1);
            var y = (Array)Operand(complexY ? typeof(Complex32[]) : typeof(float[]), 6, 1, size.Item2);
            Delegate function = complexX || complexY ? (Delegate)new IMeshComplex32((a, b) => a * b + a) : new IMeshFloat((a, b) => a * b + a);
            var result = (Array)Call("Compute", x, y, function);
            Assert.Equal(size.Item1, result.GetLength(0)); Assert.Equal(size.Item2, result.GetLength(1));
            for (int i = 0; i < size.Item1; i++) for (int j = 0; j < size.Item2; j++)
                Check(Value(x, 0, i) * Value(y, 0, j) + Value(x, 0, i), result.GetValue(i, j)!);
        }
    }

    public static IEnumerable<object[]> StatisticsCases()
    {
        foreach (int n in new[] { 2, 3, 17 })
            foreach (float scale in new[] { 1e-20f, 1f, 1e10f, 1e20f })
                foreach (float offset in new[] { 0f, 10000f }) yield return new object[] { n, scale, offset };
    }

    static void RealComplex(double expected, Complex32 actual)
    {
        Assert.Equal(0, actual.Imag);
        if (float.IsInfinity((float)expected)) Assert.Equal((float)expected, actual.Real);
        else Close(expected, actual.Real, 4 * (double)float.Epsilon, 2e-5);
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
        RealComplex(variance, v.Var()); RealComplex(variance, v.Cov()); RealComplex(Math.Sqrt(variance), v.StnDev());
        RealComplex(norm2, v.Abs(true)); RealComplex(Math.Sqrt(norm2), v.Abs());
        RealComplex(error, v.Var(other)); RealComplex(Math.Sqrt(error), v.StnDev(other));
        var matrix = new Complex32[n, 2]; var rotated = new Complex32[n, 2];
        for (int i = 0; i < n; i++) { matrix[i, 0] = v[i]; matrix[i, 1] = other[i]; rotated[i, 0] = other[i]; rotated[i, 1] = -v[i]; }
        var vars = matrix.Var(); var std = matrix.StnDev(); var errors = matrix.Var(rotated); var rms = matrix.StnDev(rotated);
        for (int j = 0; j < 2; j++) { RealComplex(variance, vars[j]); RealComplex(Math.Sqrt(variance), std[j]); RealComplex(error, errors[j]); RealComplex(Math.Sqrt(error), rms[j]); }
        var covariance = matrix.Cov();
        RealComplex(variance, covariance[0, 0]); RealComplex(variance, covariance[1, 1]);
        RealComplex(variance, new Complex32(covariance[0, 1].Imag, covariance[0, 1].Real));
        RealComplex(variance, new Complex32(-covariance[1, 0].Imag, covariance[1, 0].Real));
        var rowNorms = matrix.Abs();
        for (int i = 0; i < n; i++) RealComplex(Math.Sqrt(2) * Complex.Abs((Complex)v[i]), rowNorms[i]);
    }

    [Theory]
    [InlineData(0)] [InlineData(1)]
    public void SampleVarianceIsUndefinedWithFewerThanTwoObservations(int n)
    {
        var v = new Complex32[n];
        Assert.True(float.IsNaN(v.Var().Real)); Assert.True(float.IsNaN(v.Cov().Real)); Assert.True(float.IsNaN(v.StnDev().Real));
        Assert.Equal(0, v.Var().Imag);
    }

    [Theory]
    [InlineData(3)] [InlineData(17)]
    public void ComplexCovarianceMatchesCenteredOuterProductsAndIsPositiveSemidefinite(int n)
    {
        var a = new Complex32[n, 3];
        for (int i = 0; i < n; i++) for (int j = 0; j < 3; j++)
            a[i, j] = new Complex32((float)(Math.Sin(i * (j + 1)) + i), (float)Math.Cos(i + j * .7));
        var means = new Complex[3];
        for (int j = 0; j < 3; j++) for (int i = 0; i < n; i++) means[j] += (Complex)a[i, j] / n;
        var covariance = a.Cov();
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++)
        {
            Complex expected = 0;
            for (int i = 0; i < n; i++) expected += Complex.Conjugate((Complex)a[i, j] - means[j]) * ((Complex)a[i, k] - means[k]) / (n - 1);
            Close(expected, covariance[j, k]);
            Assert.Equal(covariance[j, k].Real, covariance[k, j].Real);
            Assert.Equal(covariance[j, k].Imag, -covariance[k, j].Imag);
        }
        Complex[] weights = { new(1, 2), new(-2, .5), new(0, -1) };
        Complex quadratic = 0;
        for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) quadratic += Complex.Conjugate(weights[j]) * (Complex)covariance[j, k] * weights[k];
        Assert.True(quadratic.Real >= 0); Close(0, quadratic.Imaginary);
    }

    [Theory]
    [InlineData(false)] [InlineData(true)]
    public void SwappingTheLastAxisIndexIsAnInvolutionForRectangularMatrices(bool complex)
    {
        foreach (var shape in new[] { (1, 7), (7, 1), (3, 5), (5, 3) })
            foreach (Direction direction in Enum.GetValues<Direction>())
            {
                var a = (Array)Operand(complex ? typeof(Complex32[,]) : typeof(float[,]), 1, shape.Item1, shape.Item2);
                var original = (Array)a.Clone();
                int last = (direction == Direction.Horizontal ? shape.Item1 : direction == Direction.Vertical ? shape.Item2 : Math.Min(shape.Item1, shape.Item2)) - 1;
                Call("Swap", a, 0, last, direction);
                for (int i = 0; i < shape.Item1; i++) for (int j = 0; j < shape.Item2; j++)
                {
                    int row = direction != Direction.Vertical && (i == 0 || i == last) ? last - i : i;
                    int col = direction != Direction.Horizontal && (j == 0 || j == last) ? last - j : j;
                    Check(Value(original, row, col), a.GetValue(i, j)!);
                }
                Call("Swap", a, 0, last, direction);
                Assert.Equal(original.Cast<object>(), a.Cast<object>());
            }
    }

    public static IEnumerable<object[]> ContingencyCases()
    {
        for (int tt = 0; tt < 4; tt++) for (int tf = 0; tf < 4; tf++) for (int ft = 0; ft < 4; ft++)
            yield return new object[] { tt, tf, ft };
    }

    [Theory]
    [MemberData(nameof(ContingencyCases))]
    public void SokalSneathMatchesIndependentContingencyCounts(int tt, int tf, int ft)
    {
        float[] p = Enumerable.Repeat(2f, tt + tf).Concat(Enumerable.Repeat(0f, ft + 3)).ToArray();
        float[] q = Enumerable.Repeat(-3f, tt).Concat(Enumerable.Repeat(0f, tf)).Concat(Enumerable.Repeat(-3f, ft)).Concat(new float[3]).ToArray();
        double expected = tt + tf + ft == 0 ? 0 : 2.0 * (tf + ft) / (tt + 2 * (tf + ft));
        var distance = new SokalSneath();
        Close(expected, distance.Compute(p, q)); Close(expected, distance.Compute(q, p));
        var pc = p.Select(x => new Complex32(0, x)).ToArray(); var qc = q.Select(x => new Complex32(0, x)).ToArray();
        RealComplex(expected, distance.Compute(pc, qc)); RealComplex(expected, distance.Compute(qc, pc));
    }
}
