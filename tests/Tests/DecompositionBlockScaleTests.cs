using System.Numerics;
using UMapx.Core;
using UMapx.Decomposition;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class DecompositionBlockScaleTests
{
    [Theory]
    [InlineData(1e-20f, 1f, 0)] [InlineData(1e-20f, 1f, 5)] [InlineData(1e-20f, 1f, 10)]
    [InlineData(1e-30f, 1e30f, 0)] [InlineData(1e-30f, 1e30f, 5)] [InlineData(1e-30f, 1e30f, 10)]
    public void ExceptionalShiftsPreserveAnIndependentRoot(float small, float large, int isolated)
    {
        // Each dense block needs exceptional QR shifts; the tiny root must never be shifted with it.
        float[,] block = {
            { -.8880402f, .3476284f, .26115668f, -.08890075f, .82201105f },
            { .78620726f, -.031750247f, -.99099845f, .11013534f, .83499014f },
            { .10187164f, -.29250166f, -.8633035f, .4379137f, .76222056f },
            { -.7247197f, -.1686885f, .8706444f, .3575709f, .28234488f },
            { .46623564f, .8510192f, -.5098249f, .91684544f, -.34769323f }
        };
        var a = new float[11, 11];
        int[] rest = Enumerable.Range(0, 11).Where(i => i != isolated).ToArray();
        a[isolated, isolated] = small;
        for (int k = 0; k < 2; k++) for (int i = 0; i < 5; i++) for (int j = 0; j < 5; j++)
            a[rest[5 * k + i], rest[5 * k + j]] = block[i, j] * large * (k == 0 ? 1 : .25f);
        var evd = EVD.Decompose(a);
        var schur = Schur.Decompose(a);
        foreach (var values in new[] { evd.D, Schur.Eigenvalues(schur.T) })
            NumericAssert.Close(new Complex(small, 0), values.MinBy(z => Complex.Abs((Complex)z - small)), 0, 2e-5);
        var rowScales = Enumerable.Range(0, 11).Select(i => (double)(i == isolated ? small : large)).ToArray();
        CheckVectors(a, Identity(11), Unpack(evd.V, evd.D), evd.D, Enumerable.Repeat(1f, 11).ToArray(), rowScales);
        CheckReconstruction(a, ComplexMatrix(schur.Q), ComplexMatrix(schur.T), ComplexMatrix(schur.Q), rowScales);
    }

    public static IEnumerable<object[]> RepeatedBlocks()
    {
        foreach (bool generalized in new[] { false, true })
            foreach (bool conjugate in new[] { false, true })
                foreach (bool first in new[] { false, true })
                    foreach (var scales in new[] { (1e-20f, 1f), (1e-30f, 1e30f) })
                        yield return new object[] { generalized, conjugate, first, scales.Item1, scales.Item2 };
    }

    [Theory, MemberData(nameof(RepeatedBlocks))]
    public void RepeatedRootVectorsHaveSmallLocalResiduals(bool generalized, bool conjugate, bool first, float small, float large)
    {
        int count = conjugate ? 4 : 2, n = count + 1, start = first ? 0 : 1, isolated = first ? count : 0;
        var a = new float[n, n];
        a[isolated, isolated] = large;
        if (conjugate)
        {
            for (int k = 0; k < 4; k += 2)
            {
                a[start + k, start + k + 1] = -small;
                a[start + k + 1, start + k] = small;
            }
            a[start, start + 2] = small;
            a[start + 1, start + 3] = small;
        }
        else
        {
            a[start, start] = a[start + 1, start + 1] = small;
            a[start, start + 1] = small;
        }
        var b = Identity(n);
        Complex32[] alpha;
        float[] beta;
        Complex[,] vectors;
        if (generalized)
        {
            var d = GEVD.Decompose(a, b);
            alpha = d.Alpha; beta = d.Beta; vectors = Unpack(d.V, alpha);
        }
        else
        {
            var d = EVD.Decompose(a);
            alpha = d.D; beta = Enumerable.Repeat(1f, n).ToArray(); vectors = Unpack(d.V, alpha);
        }
        var scales = Enumerable.Range(0, n).Select(i => (double)(i == isolated ? large : small)).ToArray();
        CheckVectors(a, b, vectors, alpha, beta, scales);
        var values = GEVD.Eigenvalues(alpha, beta);
        Assert.Equal(count, values.Count(z => ((Complex)z).Magnitude < 2 * small));
        foreach (var value in values.Where(z => ((Complex)z).Magnitude < 2 * small))
            NumericAssert.Close(small, conjugate ? Math.Abs(value.Imag) : value.Real, 0, 2e-5);
    }

    [Theory]
    [InlineData(2, false, false)] [InlineData(2, true, false)]
    [InlineData(3, false, false)] [InlineData(3, true, false)]
    [InlineData(2, false, true)] [InlineData(2, true, true)]
    public void GeneralizedBlocksRetainFiniteAndInfiniteRoots(int count, bool first, bool singular)
    {
        const float small = 1e-20f;
        int n = count + 1, start = first ? 0 : 1, isolated = first ? count : 0;
        var a = new float[n, n]; var b = new float[n, n];
        a[isolated, isolated] = b[isolated, isolated] = 1;
        for (int i = 0; i < count; i++) b[start + i, start + i] = small;
        if (count == 2)
        {
            a[start, start + 1] = -small;
            a[start + 1, start] = small;
            if (singular)
            {
                a[start, start] = a[start + 1, start + 1] = 2 * small;
                b[start + 1, start + 1] = 0;
            }
        }
        else for (int i = 0; i < count; i++) a[start + (i + 1) % count, start + i] = small;
        var rowScales = Enumerable.Range(0, n).Select(i => i == isolated ? 1.0 : small).ToArray();
        var real = GEVD.Decompose(a, b);
        var complex = GEVD.Decompose(ToComplex(a), ToComplex(b));
        CheckVectors(a, b, Unpack(real.V, real.Alpha), real.Alpha, real.Beta, rowScales);
        CheckVectors(a, b, ComplexMatrix(complex.V), complex.Alpha, complex.Beta, rowScales);
        foreach (var spectrum in new[] { (real.Alpha, real.Beta), (complex.Alpha, complex.Beta) })
        {
            Assert.Equal(singular ? 1 : 0, spectrum.Item2.Count(x => x == 0));
            var values = GEVD.Eigenvalues(spectrum.Item1, spectrum.Item2);
            Complex[] expected = singular ? new Complex[] { 1, 2.5 } : count == 2
                ? new Complex[] { 1, Complex.ImaginaryOne, -Complex.ImaginaryOne }
                : new Complex[] { 1, 1, new(-.5, Math.Sqrt(3) / 2), new(-.5, -Math.Sqrt(3) / 2) };
            foreach (var target in expected)
                NumericAssert.Close(target, values.Where(z => float.IsFinite(z.Real)).MinBy(z => Complex.Abs((Complex)z - target)), 2e-5, 2e-5);
        }
        var qr = QZ.Decompose(a, b);
        CheckReconstruction(a, ComplexMatrix(qr.Q), ComplexMatrix(qr.S), ComplexMatrix(qr.Z), rowScales);
        CheckReconstruction(b, ComplexMatrix(qr.Q), ComplexMatrix(qr.T), ComplexMatrix(qr.Z), rowScales);
        var qc = QZ.Decompose(ToComplex(a), ToComplex(b));
        CheckReconstruction(a, ComplexMatrix(qc.Q), ComplexMatrix(qc.S), ComplexMatrix(qc.Z), rowScales);
        CheckReconstruction(b, ComplexMatrix(qc.Q), ComplexMatrix(qc.T), ComplexMatrix(qc.Z), rowScales);
    }

    private static void CheckVectors(float[,] a, float[,] b, Complex[,] v, Complex32[] alpha, float[] beta, double[] rowScales)
    {
        int n = a.GetLength(0);
        for (int j = 0; j < n; j++)
        {
            double norm = Math.Sqrt(Enumerable.Range(0, n).Sum(i => v[i, j].Magnitude * v[i, j].Magnitude));
            Assert.True(double.IsFinite(norm) && norm > 0);
            for (int i = 0; i < n; i++)
            {
                Complex av = 0, bv = 0;
                double bScale = 0;
                for (int k = 0; k < n; k++) { av += a[i, k] * v[k, j]; bv += b[i, k] * v[k, j]; bScale += Math.Abs(b[i, k]); }
                double bound = (Math.Abs(beta[j]) * rowScales[i] + ((Complex)alpha[j]).Magnitude * bScale) * norm;
                double error = Complex.Abs(beta[j] * av - (Complex)alpha[j] * bv);
                Assert.True(double.IsFinite(error) && error <= 2e-5 * bound, $"Local eigenvector residual: row {i}, column {j}, error {error:G6}, bound {bound:G6}");
            }
        }
    }

    private static void CheckReconstruction(float[,] a, Complex[,] q, Complex[,] t, Complex[,] z, double[] rowScales)
    {
        int n = a.GetLength(0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        {
            Complex actual = 0;
            for (int k = 0; k < n; k++) for (int l = 0; l < n; l++) actual += q[i, k] * t[k, l] * Complex.Conjugate(z[j, l]);
            Assert.True(Complex.Abs(actual - a[i, j]) <= 2e-5 * rowScales[i]);
        }
    }

    private static Complex[,] Unpack(float[,] v, Complex32[] alpha)
    {
        int n = v.GetLength(0); var result = new Complex[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
            result[i, j] = alpha[j].Imag > 0 ? new Complex(v[i, j], v[i, j + 1]) : alpha[j].Imag < 0 ? new Complex(v[i, j - 1], -v[i, j]) : v[i, j];
        return result;
    }

    private static Complex[,] ComplexMatrix(Array a)
    {
        int n = a.GetLength(0); var result = new Complex[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) result[i, j] = a.GetValue(i, j) is float f ? f : (Complex)(Complex32)a.GetValue(i, j)!;
        return result;
    }

    private static Complex32[,] ToComplex(float[,] a)
    {
        int n = a.GetLength(0); var result = new Complex32[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) result[i, j] = a[i, j];
        return result;
    }

    private static float[,] Identity(int n)
    {
        var result = new float[n, n];
        for (int i = 0; i < n; i++) result[i, i] = 1;
        return result;
    }
}
