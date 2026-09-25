using System.Numerics;
using System.Text.Json;
using UMapx.Core;
using UMapx.Decomposition;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class DecompositionReferenceTests
{
    public static IEnumerable<object[]> ReferenceCases()
    {
        using var stream = typeof(DecompositionReferenceTests).Assembly.GetManifestResourceStream("UMapx.Tests.Data.decomposition.json")!;
        using var document = JsonDocument.Parse(stream);
        foreach (var item in document.RootElement.GetProperty("cases").EnumerateArray())
            yield return new object[]
            {
                item.GetProperty("name").GetString()!,
                item.GetRawText()
            };
    }

    [Theory, MemberData(nameof(ReferenceCases))]
    public void SpectraAndFactorsAgreeWithIndependentLapackReferences(string name, string json)
    {
        using var document = JsonDocument.Parse(json);
        var item = document.RootElement;
        int m = item.GetProperty("rows").GetInt32(), n = item.GetProperty("columns").GetInt32();
        bool real = item.GetProperty("real").GetBoolean();
        string operation = item.GetProperty("operation").GetString()!;
        var a = ReadMatrix(item.GetProperty("a"), m, n);
        var input = Single(a);
        var expected = item.GetProperty("spectrum").EnumerateArray().Select(ReadValue).ToArray();
        Complex[] values;
        if (operation == "SVD")
        {
            (Complex[,] u, float[] s, Complex[,] v) d;
            if (real)
            {
                var r = SVD.Decompose(Real(a), 100);
                d = (Work(r.U), r.S, Work(r.V));
            }
            else
            {
                var c = SVD.Decompose(input, 100);
                d = (Work(c.U), c.S, Work(c.V));
            }

            Close(a, Product(Product(d.u, Diag(d.s.Select(x => new Complex(x, 0)).ToArray())), Adjoint(d.v)));
            Orthonormal(d.u);
            Orthonormal(d.v);
            values = d.s.Select(x => new Complex(x, 0)).ToArray();
            for (int j = 0; j < values.Length; j++)
                Assert.True(Complex.Abs(values[j] - expected[j]) <= 3e-5 * expected[0].Magnitude, name);
        }
        else if (operation == "GEVD")
        {
            var b = ReadMatrix(item.GetProperty("b"), n, n);
            Complex[,] v;
            Complex32[] alpha;
            float[] beta;
            if (real)
            {
                var d = GEVD.Decompose(Real(a), Real(b));
                v = UnpackEigenvectors(d.V, d.Alpha);
                alpha = d.Alpha;
                beta = d.Beta;
            }
            else
            {
                var d = GEVD.Decompose(input, Single(b));
                v = Work(d.V);
                alpha = d.Alpha;
                beta = d.Beta;
            }

            Homogeneous(a, b, v, alpha, beta);
            values = alpha.Select((x, i) => (Complex)x / beta[i]).ToArray();
        }
        else
        {
            Complex[,] v;
            if (real)
            {
                var d = EVD.Decompose(Real(a));
                v = UnpackEigenvectors(d.V, d.D);
                values = d.D.Select(x => (Complex)x).ToArray();
            }
            else
            {
                var d = EVD.Decompose(input);
                v = Work(d.V);
                values = d.D.Select(x => (Complex)x).ToArray();
            }

            Eigenvectors(a, v, values);
            if (operation == "Hermitian")
            {
                Orthonormal(v);
                Assert.All(values, x => Assert.Equal(0, x.Imaginary));
            }

            // The same reference spectrum also checks the independent Schur entry points.
            if (real)
            {
                var d = Schur.Decompose(Real(a));
                Close(a, Product(Product(Work(d.Q), Work(d.T)), Adjoint(Work(d.Q))));
                Orthonormal(Work(d.Q));
                Match(expected, Schur.Eigenvalues(d.T).Select(x => (Complex)x).ToArray(), name);
            }
            else
            {
                var d = Schur.Decompose(input);
                Close(a, Product(Product(Work(d.Q), Work(d.T)), Adjoint(Work(d.Q))));
                Orthonormal(Work(d.Q));
                Match(expected, Schur.Eigenvalues(d.T).Select(x => (Complex)x).ToArray(), name);
            }
        }

        Match(expected, values, name);
        Assert.Equal(Single(a).Cast<Complex32>(), input.Cast<Complex32>());
    }

    [Theory]
    [InlineData(48, 32, false)]
    [InlineData(32, 48, false)]
    [InlineData(64, 64, true)]
    [InlineData(257, 17, true)]
    [InlineData(17, 257, true)]
    public void ComplexSvdPreservesKnownSpectraAndNullspaces(int m, int n, bool deficient)
    {
        int k = Math.Min(m, n);
        var u = Fourier(m, k);
        var v = Fourier(n, k);
        var spectrum = Enumerable.Range(0, k).Select(i => new Complex(deficient && i >= k / 2 ? 0 : 1 + i / 3, 0)).ToArray();
        var a = Single(Product(Product(u, Diag(spectrum)), Adjoint(v)));
        var d = SVD.Decompose(a, 100);
        Close(Work(a), Product(Product(Work(d.U), Diag(d.S.Select(x => new Complex(x, 0)).ToArray())), Adjoint(Work(d.V))));
        Orthonormal(Work(d.U));
        Orthonormal(Work(d.V));
        Match(spectrum, d.S.Select(x => new Complex(x, 0)).ToArray(), "known SVD");
        var inverse = Work(SVD.PseudoInverse(d.U, d.S, d.V));
        Close(Work(a), Product(Product(Work(a), inverse), Work(a)), 1e-4);
    }

    [Theory]
    [InlineData(32)]
    [InlineData(96)]
    public void HermitianReductionAndEvdResolveRepeatedAndSignedSpectra(int n)
    {
        var q = Fourier(n, n);
        var spectrum = Enumerable.Range(0, n).Select(i => new Complex(i / 3 - n / 6, 0)).ToArray();
        var a = Single(Product(Product(q, Diag(spectrum)), Adjoint(q)));
        MakeHermitian(a);
        var reduction = Householder.Decompose(a);
        Close(Work(a), Product(Product(Work(reduction.H), Work(reduction.T)), Adjoint(Work(reduction.H))));
        Orthonormal(Work(reduction.H));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                Assert.Equal(reduction.T[i, j], reduction.T[j, i].Conjugate);
                if (Math.Abs(i - j) > 1)
                    Assert.Equal(default, reduction.T[i, j]);
            }

        var d = EVD.Decompose(a);
        var values = d.D.Select(x => (Complex)x).ToArray();
        Match(spectrum, values, "Hermitian spectrum");
        Eigenvectors(Work(a), Work(d.V), values);
        Orthonormal(Work(d.V));
        Assert.All(d.D, x => Assert.Equal(0, x.Imag));
    }

    [Theory]
    [InlineData(1e-30f, 1e30f)]
    [InlineData(1e30f, 1e-30f)]
    public void InterleavedHermitianComponentsKeepTheirOwnScales(float small, float large)
    {
        var a = new Complex32[5, 5];
        a[0, 0] = small;
        a[2, 2] = 3 * small;
        a[0, 2] = new(0, small);
        a[2, 0] = new(0, -small);
        a[1, 1] = large;
        a[3, 3] = 3 * large;
        a[1, 3] = new(large, 0);
        a[3, 1] = a[1, 3];
        a[4, 4] = 0;
        var d = EVD.Decompose(a);
        var expected = new[]
        {
            0.0,
            (2 - Math.Sqrt(2)) * small,
            (2 + Math.Sqrt(2)) * small,
            (2 - Math.Sqrt(2)) * large,
            (2 + Math.Sqrt(2)) * large
        }.Order().ToArray();
        for (int i = 0; i < expected.Length; i++)
            Assert.True(Math.Abs(d.D[i].Real - expected[i]) <= Math.Abs(expected[i]) * 2e-5);
        Eigenvectors(Work(a), Work(d.V), d.D.Select(x => (Complex)x).ToArray());
        Orthonormal(Work(d.V));
    }

    [Theory]
    [InlineData(24, false)]
    [InlineData(48, false)]
    [InlineData(24, true)]
    public void GeneralPencilsPreserveFiniteAndInfiniteEigenpairs(int n, bool singular)
    {
        var q = Fourier(n, n);
        var a0 = new Complex[n, n];
        var b0 = new Complex[n, n];
        for (int i = 0; i < n; i++)
        {
            a0[i, i] = new Complex(i + 1, i % 3 - 1);
            b0[i, i] = singular && i % 4 == 0 ? 0 : 1 + i / (double)n;
            if (i + 1 < n)
            {
                a0[i, i + 1] = new Complex(.2, -.1);
                b0[i, i + 1] = .1;
            }
        }

        // A shared unitary similarity keeps the expected generalized spectrum known.
        var a = Single(Product(Product(q, a0), Adjoint(q)));
        var b = Single(Product(Product(q, b0), Adjoint(q)));
        if (singular)
        {
            // Keep exact zero diagonals in B: rounding a dense singular B would perturb infinity.
            a = Single(a0);
            b = Single(b0);
        }

        var d = QZ.Decompose(a, b);
        Close(Work(a), Product(Product(Work(d.Q), Work(d.S)), Adjoint(Work(d.Z))));
        Close(Work(b), Product(Product(Work(d.Q), Work(d.T)), Adjoint(Work(d.Z))));
        Orthonormal(Work(d.Q));
        Orthonormal(Work(d.Z));
        var e = GEVD.Decompose(a, b);
        Homogeneous(Work(a), Work(b), Work(e.V), e.Alpha, e.Beta);
        if (singular)
            Assert.Equal(Enumerable.Range(0, n).Count(i => i % 4 == 0), e.Beta.Count(x => x == 0));
        else
            Match(Enumerable.Range(0, n).Select(i => a0[i, i] / b0[i, i]).ToArray(), e.Alpha.Select((x, i) => (Complex)x / e.Beta[i]).ToArray(), "pencil");
    }

    [Theory]
    [InlineData("GramSchmidt", false)]
    [InlineData("GramSchmidt", true)]
    [InlineData("Arnoldi", false)]
    [InlineData("Arnoldi", true)]
    [InlineData("Lanczos", false)]
    [InlineData("Lanczos", true)]
    public void BasisCompletionPreservesBothDomainsAtBreakdown(string algorithm, bool rankOne)
    {
        const int n = 32;
        var real = new float[n, n];
        if (rankOne)
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    real[i, j] = (i + 1) * (j + 1);
        var complex = Single(Work(real));
        (float[,] q, float[,] h) r = algorithm switch
        {
            "GramSchmidt" => GramSchmidt.Decompose(real),
            "Arnoldi" => Arnoldi.Decompose(real),
            _ => Lanczos.Decompose(real)
        };
        (Complex32[,] q, Complex32[,] h) c = algorithm switch
        {
            "GramSchmidt" => GramSchmidt.Decompose(complex),
            "Arnoldi" => Arnoldi.Decompose(complex),
            _ => Lanczos.Decompose(complex)
        };
        Orthonormal(Work(r.q));
        Orthonormal(Work(c.q));
        var realProduct = Product(Work(r.q), Work(r.h));
        var complexProduct = Product(Work(c.q), Work(c.h));
        if (algorithm != "GramSchmidt")
        {
            realProduct = Product(realProduct, Adjoint(Work(r.q)));
            complexProduct = Product(complexProduct, Adjoint(Work(c.q)));
        }

        Close(Work(real), realProduct);
        Close(Work(complex), complexProduct);
    }

    private static void Homogeneous(Complex[,] a, Complex[,] b, Complex[,] v, Complex32[] alpha, float[] beta)
    {
        var av = Product(a, v);
        var bv = Product(b, v);
        for (int j = 0; j < alpha.Length; j++)
        {
            double error = 0, vectorNorm = 0;
            for (int i = 0; i < alpha.Length; i++)
            {
                error += Math.Pow(Complex.Abs(beta[j] * av[i, j] - (Complex)alpha[j] * bv[i, j]), 2);
                vectorNorm += v[i, j].Magnitude * v[i, j].Magnitude;
            }

            double bound = (Math.Abs(beta[j]) * Norm(a) + ((Complex)alpha[j]).Magnitude * Norm(b)) * Math.Sqrt(vectorNorm);
            Assert.True(vectorNorm > 0 && Math.Sqrt(error) <= 5e-5 * bound, $"Homogeneous residual {Math.Sqrt(error) / bound}");
        }
    }

    private static void Eigenvectors(Complex[,] a, Complex[,] v, Complex[] values)
    {
        var av = Product(a, v);
        for (int j = 0; j < values.Length; j++)
        {
            double error = 0, vn = 0;
            for (int i = 0; i < values.Length; i++)
            {
                error += Math.Pow(Complex.Abs(av[i, j] - v[i, j] * values[j]), 2);
                vn += v[i, j].Magnitude * v[i, j].Magnitude;
            }

            Assert.True(vn > 0 && Math.Sqrt(error) <= 4e-5 * (Norm(a) + values[j].Magnitude) * Math.Sqrt(vn));
        }
    }

    private static void Match(Complex[] expected, Complex[] actual, string name)
    {
        Assert.Equal(expected.Length, actual.Length);
        var remaining = actual.ToList();
        double scale = expected.Max(x => x.Magnitude);
        foreach (var value in expected.OrderBy(x => x.Real).ThenBy(x => x.Imaginary))
        {
            int best = Enumerable.Range(0, remaining.Count).MinBy(i => Complex.Abs(remaining[i] - value));
            Assert.True(Complex.Abs(remaining[best] - value) <= 8e-5 * scale, $"{name}: expected {value}, got {remaining[best]}");
            remaining.RemoveAt(best);
        }
    }

    private static Complex[,] UnpackEigenvectors(float[,] v, Complex32[] eigenvalues)
    {
        var result = Work(v);
        for (int j = 0; j < eigenvalues.Length; j++)
        {
            if (eigenvalues[j].Imag <= 0)
                continue;
            for (int i = 0; i < v.GetLength(0); i++)
            {
                result[i, j] = new(v[i, j], v[i, j + 1]);
                result[i, j + 1] = Complex.Conjugate(result[i, j]);
            }

            j++;
        }

        return result;
    }

    private static Complex ReadValue(JsonElement x) => new(x[0].GetDouble(), x[1].GetDouble());
    private static Complex[,] ReadMatrix(JsonElement values, int m, int n)
    {
        var a = new Complex[m, n];
        int k = 0;
        foreach (var x in values.EnumerateArray())
        {
            a[k / n, k % n] = ReadValue(x);
            k++;
        }

        return a;
    }

    private static Complex[,] Fourier(int m, int n)
    {
        var a = new Complex[m, n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                a[i, j] = Complex.FromPolarCoordinates(1 / Math.Sqrt(m), 2 * Math.PI * i * j / m);
        return a;
    }

    private static void MakeHermitian(Complex32[,] a)
    {
        for (int i = 0; i < a.GetLength(0); i++)
        {
            a[i, i] = a[i, i].Real;
            for (int j = 0; j < i; j++)
                a[j, i] = a[i, j].Conjugate;
        }
    }

    private static Complex[,] Work(Complex32[,] a)
    {
        var b = new Complex[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                b[i, j] = (Complex)a[i, j];
        return b;
    }

    private static Complex[,] Work(float[,] a)
    {
        var b = new Complex[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                b[i, j] = a[i, j];
        return b;
    }

    private static Complex32[,] Single(Complex[,] a)
    {
        var b = new Complex32[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                b[i, j] = new((float)a[i, j].Real, (float)a[i, j].Imaginary);
        return b;
    }

    private static float[,] Real(Complex[,] a)
    {
        var b = new float[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                b[i, j] = (float)a[i, j].Real;
        return b;
    }

    private static Complex[,] Diag(Complex[] d)
    {
        var a = new Complex[d.Length, d.Length];
        for (int i = 0; i < d.Length; i++)
            a[i, i] = d[i];
        return a;
    }

    private static Complex[,] Adjoint(Complex[,] a)
    {
        var b = new Complex[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                b[j, i] = Complex.Conjugate(a[i, j]);
        return b;
    }

    private static Complex[,] Product(Complex[,] a, Complex[,] b)
    {
        var c = new Complex[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < c.GetLength(0); i++)
            for (int j = 0; j < c.GetLength(1); j++)
                for (int k = 0; k < a.GetLength(1); k++)
                    c[i, j] += a[i, k] * b[k, j];
        return c;
    }

    private static double Norm(Complex[,] a) => Math.Sqrt(a.Cast<Complex>().Sum(x => x.Magnitude * x.Magnitude));
    private static void Close(Complex[,] a, Complex[,] b, double tolerance = 2e-5)
    {
        double error = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                error += Math.Pow(Complex.Abs(a[i, j] - b[i, j]), 2);
        Assert.True(Math.Sqrt(error) <= tolerance * Norm(a), $"Relative residual {Math.Sqrt(error) / Norm(a)}");
    }

    private static void Orthonormal(Complex[,] q) => Close(Diag(Enumerable.Repeat(Complex.One, q.GetLength(1)).ToArray()), Product(Adjoint(q), q));
}
