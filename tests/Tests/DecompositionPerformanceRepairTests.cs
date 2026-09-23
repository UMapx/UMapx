using UMapx.Core;
using UMapx.Decomposition;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class DecompositionPerformanceRepairTests
{
    [Theory]
    [InlineData(32, 1e-30f)] [InlineData(32, 1f)] [InlineData(32, 1e30f)]
    [InlineData(96, 1f)]
    public void RealKernelsRetainAccuracyAtLargerOrdersAndExtremeScales(int n, float scale)
    {
        var a = Sample(n, n, 173, scale, symmetric: true);
        var original = (float[,])a.Clone();
        var cholesky = Cholesky.Decompose(a);
        Relative(a, Product(Work(cholesky), Transpose(Work(cholesky))));
        var ldl = LDL.Decompose(a);
        Relative(a, Product(Product(Work(ldl.L), Diagonal(ldl.D)), Transpose(Work(ldl.L))));
        var udl = UDL.Decompose(a);
        Relative(a, Product(Product(Work(udl.U), Diagonal(udl.D)), Transpose(Work(udl.U))));
        var householder = Householder.Decompose(a);
        Similarity(a, householder.H, householder.T);
        var hessenberg = Hessenberg.Decompose(a);
        Similarity(a, hessenberg.P, hessenberg.H);
        var arnoldi = Arnoldi.Decompose(a);
        Similarity(a, arnoldi.Q, arnoldi.H);
        var gram = GramSchmidt.Decompose(a);
        Relative(a, Product(Work(gram.Q), Work(gram.R))); Orthogonal(gram.Q);
        var lu = LU.Decompose(a);
        var permuted = new float[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) permuted[i, j] = a[lu.P[i], j];
        Relative(permuted, Product(Work(lu.L), Work(lu.U)));
        Assert.Equal(original.Cast<float>(), a.Cast<float>());
    }

    [Theory]
    [InlineData(128, 12, 1e-30f)] [InlineData(12, 128, 1e30f)]
    [InlineData(40, 40, 1f)] [InlineData(128, 12, 1f)]
    public void RealRectangularKernelsProduceEconomyAndFullFactors(int m, int n, float scale)
    {
        var a = Sample(m, n, 317, scale);
        var qr = QR.Decompose(a);
        Assert.Equal(Math.Min(m, n), qr.Q.GetLength(1));
        Relative(a, Product(Work(qr.Q), Work(qr.R))); Orthogonal(qr.Q);
        var svd = SVD.Decompose(a, 100);
        Relative(a, Product(Product(Work(svd.U), Diagonal(svd.S)), Transpose(Work(svd.V))));
        Orthogonal(svd.U); Orthogonal(svd.V);
        var bidiagonal = Bidiagonal.Decompose(a);
        Relative(a, Product(Product(Work(bidiagonal.U), Work(bidiagonal.B)), Transpose(Work(bidiagonal.V))));
        Orthogonal(bidiagonal.U); Orthogonal(bidiagonal.V);
    }

    public static IEnumerable<object[]> Pencils()
    {
        foreach (int n in new[] { 1, 2, 3, 12, 40 })
            foreach (string kind in new[] { "Dense", "ZeroA", "ZeroB", "BothZero", "SingularB", "SharedNullspace" })
                foreach (var scale in new[] { (1f, 1f), (1e30f, 1e-30f), (1e-30f, 1e30f) })
                    yield return new object[] { n, kind, scale.Item1, scale.Item2 };
    }

    [Theory, MemberData(nameof(Pencils))]
    public void AccumulatedQzFactorsHandleSingularPencilsAndIndependentUnits(int n, string kind, float scaleA, float scaleB)
    {
        var a = Sample(n, n, 719, scaleA); var b = Sample(n, n, 827, scaleB);
        if (kind is "ZeroA" or "BothZero") Array.Clear(a);
        if (kind is "ZeroB" or "BothZero") Array.Clear(b);
        if (kind == "SingularB") for (int i = 0; i < n; i++) b[i, n - 1] = 0;
        if (kind == "SharedNullspace") for (int j = 0; j < n; j++) a[n - 1, j] = b[n - 1, j] = 0;
        var originalA = (float[,])a.Clone(); var originalB = (float[,])b.Clone();
        var d = QZ.Decompose(a, b);
        Relative(a, Product(Product(Work(d.Q), Work(d.S)), Transpose(Work(d.Z))));
        Relative(b, Product(Product(Work(d.Q), Work(d.T)), Transpose(Work(d.Z))));
        Orthogonal(d.Q); Orthogonal(d.Z);
        for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) Assert.Equal(0, d.T[i, j]);
        Assert.Equal(originalA.Cast<float>(), a.Cast<float>());
        Assert.Equal(originalB.Cast<float>(), b.Cast<float>());
    }

    [Theory]
    [InlineData(32)] [InlineData(128)]
    public void LanczosReorthogonalizesWhenProjectionCancelsMostOfTheVector(int n)
    {
        var a = Sample(n, n, 173 + 13 * n + n, 1, symmetric: true);
        // A cluster of eigenvalues near n+1 caused 25-30% reconstruction error with one pass.
        var d = Lanczos.Decompose(a);
        Similarity(a, d.Q, d.T);
        var complex = new Complex32[n, n];
        // A diagonal unitary similarity produces a genuinely complex Hermitian matrix
        // with the same clustered spectrum, testing the shared adaptive rule.
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        {
            double angle = (i - j) * .37;
            complex[i, j] = new Complex32((float)(a[i, j] * Math.Cos(angle)), (float)(a[i, j] * Math.Sin(angle)));
        }
        var c = Lanczos.Decompose(complex);
        var reconstructed = ProductComplex(ProductComplex(ComplexWork(c.Q), ComplexWork(c.T)), Adjoint(ComplexWork(c.Q)));
        double error = 0, norm = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++)
        {
            double delta = System.Numerics.Complex.Abs((System.Numerics.Complex)complex[i, j] - reconstructed[i, j]);
            error += delta * delta; norm += (double)a[i, j] * a[i, j];
        }
        Assert.True(Math.Sqrt(error / norm) < 5e-6);
    }

    [Theory]
    [InlineData(16, false)] [InlineData(16, true)] [InlineData(40, false)]
    public void RealGsvdPreservesComplementarySubspacesBeforeNarrowing(int n, bool zero)
    {
        var a = Sample(n + 7, n, 1103, 1); var b = Sample(n + 11, n, 1201, 1);
        for (int i = 0; i < a.GetLength(0); i++) a[i, 0] = 0;
        for (int i = 0; i < b.GetLength(0); i++) for (int j = 1; j < n; j++) b[i, j] *= zero ? 0 : 1e-10f;
        var d = GSVD.Decompose(a, b, 100);
        Relative(a, Product(Product(Work(d.U1), Diagonal(d.S1)), Work(d.X)));
        Relative(b, Product(Product(Work(d.U2), Diagonal(d.S2)), Work(d.X)));
        Orthogonal(d.U1); Orthogonal(d.U2);
    }

    [Fact]
    public void NmfReusesItsIterationWorkspaceAndStillImprovesTheApproximation()
    {
        var a = Sample(32, 24, 1297, 1);
        for (int i = 0; i < 32; i++) for (int j = 0; j < 24; j++) a[i, j] = Math.Abs(a[i, j]);
        var initial = NMF.Decompose(a, 4, 1); var refined = NMF.Decompose(a, 4, 100);
        Assert.True(Error(a, Product(Work(refined.W), Work(refined.H))) < Error(a, Product(Work(initial.W), Work(initial.H))));
        long Measure(int iterations)
        {
            long before = GC.GetAllocatedBytesForCurrentThread();
            var d = NMF.Decompose(a, 4, iterations);
            long result = GC.GetAllocatedBytesForCurrentThread() - before;
            GC.KeepAlive(d);
            return result;
        }
        long shortRun = Measure(2), longRun = Measure(100);
        // A workspace should depend on dimensions, not the number of updates.
        Assert.True(longRun <= shortRun + 4096, $"Two iterations allocated {shortRun} bytes; 100 allocated {longRun}.");
    }

    private static float[,] Sample(int m, int n, int seed, float scale, bool symmetric = false)
    {
        var random = new Random(seed); var a = new float[m, n];
        for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) a[i, j] = (float)(2 * random.NextDouble() - 1) * scale;
        if (symmetric) for (int i = 0; i < n; i++) for (int j = 0; j <= i; j++) a[i, j] = a[j, i] = i == j ? (n + 1) * scale : a[i, j];
        return a;
    }

    private static void Similarity(float[,] a, float[,] q, float[,] t)
    {
        Relative(a, Product(Product(Work(q), Work(t)), Transpose(Work(q)))); Orthogonal(q);
    }
    private static void Orthogonal(float[,] q)
    {
        var actual = Product(Transpose(Work(q)), Work(q)); int n = actual.GetLength(0);
        double error = 0;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) error += Math.Pow(actual[i, j] - (i == j ? 1 : 0), 2);
        Assert.True(Math.Sqrt(error / n) < 5e-6, $"Orthogonality error {Math.Sqrt(error / n)}.");
    }
    private static void Relative(float[,] a, double[,] b) => Assert.True(Error(a, b) < 5e-6);
    private static double Error(float[,] a, double[,] b)
    {
        Assert.Equal(a.GetLength(0), b.GetLength(0)); Assert.Equal(a.GetLength(1), b.GetLength(1));
        double error = 0, norm = 0;
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++)
        {
            Assert.True(double.IsFinite(b[i, j]));
            double delta = a[i, j] - b[i, j]; error += delta * delta; norm += (double)a[i, j] * a[i, j];
        }
        return norm == 0 ? Math.Sqrt(error) : Math.Sqrt(error / norm);
    }
    private static double[,] Work(float[,] a)
    {
        var b = new double[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[i, j] = a[i, j]; return b;
    }
    private static double[,] Diagonal(float[] a)
    {
        var b = new double[a.Length, a.Length]; for (int i = 0; i < a.Length; i++) b[i, i] = a[i]; return b;
    }
    private static double[,] Transpose(double[,] a)
    {
        var b = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[j, i] = a[i, j]; return b;
    }
    private static double[,] Product(double[,] a, double[,] b)
    {
        Assert.Equal(a.GetLength(1), b.GetLength(0)); var c = new double[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < c.GetLength(0); i++) for (int k = 0; k < a.GetLength(1); k++)
            for (int j = 0; j < c.GetLength(1); j++) c[i, j] += a[i, k] * b[k, j]; return c;
    }
    private static System.Numerics.Complex[,] ComplexWork(Complex32[,] a)
    {
        var b = new System.Numerics.Complex[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[i, j] = (System.Numerics.Complex)a[i, j]; return b;
    }
    private static System.Numerics.Complex[,] Adjoint(System.Numerics.Complex[,] a)
    {
        var b = new System.Numerics.Complex[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[j, i] = System.Numerics.Complex.Conjugate(a[i, j]); return b;
    }
    private static System.Numerics.Complex[,] ProductComplex(System.Numerics.Complex[,] a, System.Numerics.Complex[,] b)
    {
        var c = new System.Numerics.Complex[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < c.GetLength(0); i++) for (int k = 0; k < a.GetLength(1); k++)
            for (int j = 0; j < c.GetLength(1); j++) c[i, j] += a[i, k] * b[k, j]; return c;
    }
}
