using System.Numerics;
using UMapx.Core;
using UMapx.Decomposition;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class ComplexDecompositionTests
{
    [Theory]
    [InlineData(1e-4f)] [InlineData(1e-10f)] [InlineData(1e-20f)]
    [InlineData(1e-30f)] [InlineData(1e10f)] [InlineData(1e20f)] [InlineData(1e30f)]
    public void GsvdPreservesSmallInputsAndBothOrthonormalBases(float scale)
    {
        foreach (bool complex in new[] { false, true })
        foreach (bool swap in new[] { false, true })
        {
            var a = new Complex32[,] { { 1, 0 }, { 0, 1 } };
            var b = new Complex32[,] { { scale, new(scale, complex ? scale : 0) }, { 0, 2 * scale } };
            if (swap) (a, b) = (b, a);
            var d = GSVD.Decompose(a, b, 100);
            Relative(Work(a), Product(Product(Work(d.U1), Diagonal(d.S1)), Work(d.X)));
            Relative(Work(b), Product(Product(Work(d.U2), Diagonal(d.S2)), Work(d.X)));
            Orthonormal(d.U1); Orthonormal(d.U2);
            Assert.All(d.S2, x => Assert.True(x > 0));
            if (!complex)
            {
                var ra = new float[,] { { 1, 0 }, { 0, 1 } };
                var rb = new float[,] { { scale, scale }, { 0, 2 * scale } };
                if (swap) (ra, rb) = (rb, ra);
                var real = GSVD.Decompose(ra, rb, 100);
                for (int j = 0; j < 2; j++) NumericAssert.Close(d.S2[j], real.S2[j], 0, 2e-6);
            }
        }
    }

    [Theory]
    [InlineData("Rotation,1e30,1e-16")]
    [InlineData("Rotation,1e-30,1e-16")]
    [InlineData("Symmetric,1e30,1e-16")]
    [InlineData("Symmetric,1e-30,1e-16")]
    [InlineData("Triangular,1,0")]
    public async Task RealEigenvaluesPreserveScaleAndTerminateAtZeroTolerance(string argument)
    {
        Assert.Equal("True", await AuditProcess.RunAsync("EigenScale", argument));
    }

    [Theory]
    [InlineData(1e-8f)] [InlineData(1e-30f)]
    public void ComplexEigenvaluesRetainSmallImaginaryComponents(float imaginary)
    {
        var a = new Complex32[,] { { 1, 0 }, { 0, new(0, imaginary) } };
        var d = EVD.Decompose(a);
        Assert.Contains(d.D, z => z.Real == 0 && z.Imag == imaginary);
        a[1, 1] = new(1, imaginary);
        d = EVD.Decompose(a);
        Assert.Contains(d.D, z => z.Real == 1 && z.Imag == imaginary);
    }

    public static IEnumerable<object[]> RealEigenCases()
    {
        foreach (int n in new[] { 1, 2, 5, 12, 24 })
        foreach (bool symmetric in new[] { false, true })
        foreach (float scale in new[] { 1e-30f, 1f, 1e30f })
        foreach (float eps in new[] { 0f, 1e-16f, 1e-7f })
            yield return new object[] { n, symmetric, scale, eps };
    }

    [Theory, MemberData(nameof(RealEigenCases))]
    public void RealEigenvectorsSatisfyScaledEquations(int n, bool symmetric, float scale, float eps)
    {
        var a = NumericAssert.Matrix(n, n);
        if (symmetric)
            for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) a[i, j] = a[j, i];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) a[i, j] *= scale;
        var original = (float[,])a.Clone();
        var d = EVD.Decompose(a, eps);
        var v = Work(d.V);
        var blocks = new Complex[n, n];
        for (int j = 0; j < n; j++)
        {
            blocks[j, j] = d.D[j].Real;
            if (d.D[j].Imag > 0) blocks[j, j + 1] = d.D[j].Imag;
            if (d.D[j].Imag < 0) blocks[j, j - 1] = d.D[j].Imag;
            Assert.True(Enumerable.Range(0, n).Sum(i => v[i, j].Magnitude) > 0);
        }
        Relative(Product(Work(a), v), Product(v, blocks), 2e-5);
        if (symmetric)
            Relative(Diagonal(Enumerable.Repeat(1f, n).ToArray()), Product(Adjoint(v), v));
        Assert.Equal(original.Cast<float>(), a.Cast<float>());
    }

    public static IEnumerable<object[]> RectangularCases()
    {
        foreach (string algorithm in new[] { "QR", "LQ", "QL", "RQ", "SVD", "Bidiagonal", "Polar", "GramSchmidt" })
            foreach (var shape in new[] { (1, 1), (1, 5), (5, 1), (3, 5), (5, 3), (5, 5), (9, 7) })
                foreach (string kind in new[] { "Dense", "Zero", "RankOne", "Diagonal" })
                    foreach (float scale in new[] { 1e-30f, 1f, 1e30f })
                    {
                        if (algorithm == "GramSchmidt" && shape.Item1 < shape.Item2) continue;
                        yield return new object[] { algorithm, shape.Item1, shape.Item2, kind, scale };
                    }
    }

    [Theory, MemberData(nameof(RectangularCases))]
    public void RectangularFactorsReconstructAndHaveUnitaryBases(string algorithm, int m, int n, string kind, float scale)
    {
        var a = Sample(m, n, kind, scale);
        var original = (Complex32[,])a.Clone();
        var expected = Work(a);
        Complex[,] actual;
        switch (algorithm)
        {
            case "QR":
            {
                var (q, r) = QR.Decompose(a);
                actual = Product(Work(q), Work(r)); Orthonormal(q);
                Band(r, 0, n); break;
            }
            case "LQ":
            {
                var (l, q) = LQ.Decompose(a);
                actual = Product(Work(l), Work(q)); Orthonormal(q, false);
                Band(l, m, 0); break;
            }
            case "QL":
            {
                var (q, l) = QL.Decompose(a);
                actual = Product(Work(q), Work(l)); Orthonormal(q);
                Band(l, m, Math.Max(0, n - m)); break;
            }
            case "RQ":
            {
                var (r, q) = RQ.Decompose(a);
                actual = Product(Work(r), Work(q)); Orthonormal(q, false);
                Band(r, Math.Max(0, m - n), n); break;
            }
            case "GramSchmidt":
            {
                var (q, r) = GramSchmidt.Decompose(a);
                actual = Product(Work(q), Work(r)); Orthonormal(q);
                Band(r, 0, n); break;
            }
            case "Bidiagonal":
            {
                var (u, b, v) = Bidiagonal.Decompose(a);
                actual = Product(Product(Work(u), Work(b)), Adjoint(Work(v)));
                Orthonormal(u); Orthonormal(v); Band(b, 0, 1); break;
            }
            case "Polar":
            {
                var (u, p) = Polar.Decompose(a, 100);
                actual = Product(Work(u), Work(p)); Relative(Work(p), Adjoint(Work(p)));
                // A partial isometry satisfies U U^H U = U even for rectangular, rank-deficient inputs.
                var w = Work(u); Relative(w, Product(Product(w, Adjoint(w)), w)); break;
            }
            default:
            {
                var (u, s, v) = SVD.Decompose(a, 100);
                actual = Product(Product(Work(u), Diagonal(s)), Adjoint(Work(v)));
                Orthonormal(u); Orthonormal(v);
                Assert.Equal(Math.Min(m, n), s.Length);
                Assert.All(s, x => Assert.True(float.IsFinite(x) && x >= 0));
                for (int i = 1; i < s.Length; i++) Assert.True(s[i - 1] >= s[i]);
                var inverse = Work(SVD.PseudoInverse(u, s, v));
                var ap = Product(expected, inverse); var pa = Product(inverse, expected);
                Relative(expected, Product(ap, expected), 2e-4);
                Relative(inverse, Product(pa, inverse), 2e-4);
                Relative(ap, Adjoint(ap), 2e-4); Relative(pa, Adjoint(pa), 2e-4);
                break;
            }
        }
        Relative(expected, actual);
        Assert.Equal(original.Cast<Complex32>(), a.Cast<Complex32>());
    }

    public static IEnumerable<object[]> SquareCases()
    {
        foreach (string algorithm in new[] { "LU", "LDU", "Diagonal", "Hessenberg", "Arnoldi", "Schur", "EVD" })
            foreach (int n in new[] { 1, 2, 3, 5, 8 })
                foreach (string kind in new[] { "Dense", "Diagonal", "Zero", "RankOne", "Jordan" })
                    foreach (float scale in new[] { 1e-30f, 1f, 1e30f })
                    {
                        if ((algorithm is "LDU" or "Diagonal") && (kind is "Zero" or "RankOne")) continue;
                        yield return new object[] { algorithm, n, kind, scale };
                    }
    }

    [Theory, MemberData(nameof(SquareCases))]
    public void SquareFactorsRespectTheirIdentities(string algorithm, int n, string kind, float scale)
    {
        var a = Sample(n, n, kind, scale);
        var original = (Complex32[,])a.Clone();
        var expected = Work(a);
        Complex[,] actual;
        switch (algorithm)
        {
            case "LU":
            {
                var (l, u, p) = LU.Decompose(a);
                expected = Rows(expected, p); actual = Product(Work(l), Work(u));
                Band(l, n, 0); Band(u, 0, n); break;
            }
            case "LDU":
            {
                var (l, d, u, p) = LDU.Decompose(a);
                expected = Rows(expected, p); actual = Product(Product(Work(l), Diagonal(d)), Work(u));
                Band(l, n, 0); Band(u, 0, n); break;
            }
            case "Diagonal":
            {
                var (b, d) = UMapx.Decomposition.Diagonal.Decompose(a);
                actual = Product(Work(b), Diagonal(d)); break;
            }
            case "Hessenberg":
            {
                var (p, h) = Hessenberg.Decompose(a); Orthonormal(p); Band(h, 1, n);
                actual = Product(Product(Work(p), Work(h)), Adjoint(Work(p))); break;
            }
            case "Arnoldi":
            {
                var (q, h) = Arnoldi.Decompose(a); Orthonormal(q); Band(h, 1, n);
                actual = Product(Product(Work(q), Work(h)), Adjoint(Work(q))); break;
            }
            case "Schur":
            {
                var (q, t) = Schur.Decompose(a); Orthonormal(q); Band(t, 0, n);
                actual = Product(Product(Work(q), Work(t)), Adjoint(Work(q)));
                Assert.Equal(Enumerable.Range(0, n).Select(i => t[i, i]), Schur.Eigenvalues(t)); break;
            }
            default:
            {
                var (v, d) = EVD.Decompose(a);
                actual = Product(Work(v), Diagonal(d)); expected = Product(expected, Work(v));
                Assert.All(v.Cast<Complex32>(), z => Assert.True(float.IsFinite(z.Real) && float.IsFinite(z.Imag)));
                break;
            }
        }
        Relative(expected, actual, 1e-4);
        Assert.Equal(original.Cast<Complex32>(), a.Cast<Complex32>());
    }

    [Theory]
    [InlineData(1)] [InlineData(2)] [InlineData(5)] [InlineData(8)]
    public void HermitianFactorsAndReflectionsUseConjugation(int n)
    {
        var x = Work(Sample(n, n, "Dense", 1));
        var positive = Product(x, Adjoint(x));
        for (int i = 0; i < n; i++) positive[i, i] += 1;
        var a = Single(positive);
        var eigen = EVD.Decompose(a);
        Orthonormal(eigen.V);
        Assert.All(eigen.D, value => Assert.Equal(0, value.Imag));
        Relative(Work(a), Product(Product(Work(eigen.V), Diagonal(eigen.D)), Adjoint(Work(eigen.V))));
        var l = Cholesky.Decompose(a);
        Relative(Work(a), Product(Work(l), Work(Cholesky.UpperFactor(l))));
        for (int i = 0; i < n; i++) { Assert.Equal(0, l[i, i].Imag); Assert.True(l[i, i].Real > 0); }
        var ldl = LDL.Decompose(a);
        Relative(Work(a), Product(Product(Work(ldl.L), Diagonal(ldl.D)), Work(LDL.UpperFactor(ldl.L))));
        var udl = UDL.Decompose(a);
        Relative(Work(a), Product(Product(Work(udl.U), Diagonal(udl.D)), Work(UDL.LowerFactor(udl.U))));
        var house = Householder.Decompose(a);
        Relative(Work(a), Product(Product(Work(house.H), Work(house.T)), Adjoint(Work(house.H))));
        Orthonormal(house.H); Band(house.T, 1, 1);
        foreach (bool full in new[] { false, true })
        {
            var lanczos = Lanczos.Decompose(a, full);
            Relative(Work(a), Product(Product(Work(lanczos.Q), Work(lanczos.T)), Adjoint(Work(lanczos.Q))));
            Orthonormal(lanczos.Q); Band(lanczos.T, 1, 1);
        }
        var vector = Enumerable.Range(0, n).Select(i => new Complex32(i + 1, 2 - i)).ToArray();
        var reflection = Householder.Reflection(vector); Orthonormal(reflection);
        var column = new Complex[n, 1]; for (int i = 0; i < n; i++) column[i, 0] = (Complex)vector[i];
        var mapped = Product(Work(reflection), column);
        for (int i = 1; i < n; i++) Assert.True(Complex.Abs(mapped[i, 0]) <= 1e-5 * Norm(column));
    }

    public static IEnumerable<object[]> PencilCases()
    {
        foreach (int n in new[] { 1, 2, 3, 5, 8 })
            foreach (string kind in new[] { "Dense", "ZeroA", "ZeroB", "BothZero", "SingularB", "RankOneB", "Triangular" })
                foreach (float scale in new[] { 1e-25f, 1f, 1e25f })
                    yield return new object[] { n, kind, scale };
    }

    [Theory, MemberData(nameof(PencilCases))]
    public void GeneralizedSchurAndEigenvectorsSatisfyHomogeneousEquations(int n, string kind, float scale)
    {
        var a = Sample(n, n, "Dense", scale);
        var b = Sample(n, n, "Dense", scale, 31);
        if (kind is "ZeroA" or "BothZero") a = new Complex32[n, n];
        if (kind is "ZeroB" or "BothZero") b = new Complex32[n, n];
        if (kind == "RankOneB") b = Sample(n, n, "RankOne", scale);
        if (kind == "SingularB")
        {
            b = Sample(n, n, "Diagonal", scale);
            b[n / 2, n / 2] = 0;
        }
        if (kind == "Triangular")
            for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) { a[i, j] = 0; b[i, j] = 0; }
        var originalA = (Complex32[,])a.Clone(); var originalB = (Complex32[,])b.Clone();
        var (q, s, t, z) = QZ.Decompose(a, b);
        Orthonormal(q); Orthonormal(z); Band(s, 0, n); Band(t, 0, n);
        Relative(Work(a), Product(Product(Work(q), Work(s)), Adjoint(Work(z))), 1e-4);
        Relative(Work(b), Product(Product(Work(q), Work(t)), Adjoint(Work(z))), 1e-4);
        for (int i = 0; i < n; i++) { Assert.Equal(0, t[i, i].Imag); Assert.True(t[i, i].Real >= 0); }
        var (v, alpha, beta) = GEVD.Decompose(a, b);
        var av = Product(Work(a), Work(v)); var bv = Product(Work(b), Work(v));
        for (int j = 0; j < n; j++)
        {
            double residual = 0, reference = 0;
            double normA = Norm(Work(a)), normB = Norm(Work(b));
            for (int i = 0; i < n; i++)
            {
                var error = beta[j] * av[i, j] - (Complex)alpha[j] * bv[i, j];
                residual += error.Magnitude * error.Magnitude;
                reference += ((Complex)v[i, j]).Magnitude * ((Complex)v[i, j]).Magnitude;
            }
            double bound = (Math.Abs(beta[j]) * normA + ((Complex)alpha[j]).Magnitude * normB) * Math.Sqrt(reference);
            Assert.True(Math.Sqrt(residual) <= 1e-4 * bound + 1e-290, $"Homogeneous residual {Math.Sqrt(residual)} exceeds {bound} for column {j}.");
        }
        Assert.Equal(originalA.Cast<Complex32>(), a.Cast<Complex32>());
        Assert.Equal(originalB.Cast<Complex32>(), b.Cast<Complex32>());
    }

    [Theory]
    [InlineData(1, 1, 1)] [InlineData(5, 4, 3)] [InlineData(8, 7, 5)]
    public void ComplexGsvdReconstructsBothInputsAndCompletesZeroSineColumns(int m, int p, int n)
    {
        foreach (bool zeroB in new[] { false, true })
        {
            var a = Sample(m, n, "Dense", 1);
            var b = zeroB ? new Complex32[p, n] : Sample(p, n, "Dense", 2);
            var d = GSVD.Decompose(a, b, 100);
            Relative(Work(a), Product(Product(Work(d.U1), Diagonal(d.S1)), Work(d.X)));
            Relative(Work(b), Product(Product(Work(d.U2), Diagonal(d.S2)), Work(d.X)));
            Orthonormal(d.U1); Orthonormal(d.U2);
            Assert.All(GSVD.Identity(d.S1, d.S2), value => NumericAssert.Close(1, value, 1e-6, 0));
            if (zeroB) Assert.All(GSVD.GeneralizedSingularValues(d.S1, d.S2), x => Assert.True(float.IsPositiveInfinity(x)));
        }
    }

    [Fact]
    public void TinyIsolatedValuesAndDifferentPencilUnitsArePreserved()
    {
        var diagonal = new Complex32[,] { { new(1e30f, 1e30f), 0 }, { 0, new(1e-30f, -1e-30f) } };
        var svd = SVD.Decompose(diagonal);
        NumericAssert.Close(Math.Sqrt(2) * 1e30, svd.S[0], 0, 2e-6);
        NumericAssert.Close(Math.Sqrt(2) * 1e-30, svd.S[1], 0, 2e-6);
        var a = new Complex32[,] { { new(1e30f, 1e30f), 0 }, { 0, new(2e30f, -1e30f) } };
        var b = new Complex32[,] { { 1e-30f, 0 }, { 0, 2e-30f } };
        var d = GEVD.Decompose(a, b);
        Assert.False(GEVD.IsSingular(d.Beta));
        for (int i = 0; i < 2; i++)
        {
            var expected = (Complex)a[i, i] / (Complex)b[i, i];
            var actual = (Complex)d.Alpha[i] / d.Beta[i];
            Assert.True(Complex.Abs(expected - actual) <= 2e-6 * Complex.Abs(expected));
        }
    }

    [Fact]
    public void ComplexPowerAndStaticApiAreUsableWithoutResultObjects()
    {
        var a = new Complex32[,] { { new(3, 4), 0 }, { 0, new(0, 1) } };
        var (v, d) = Power.Decompose(a, 100);
        NumericAssert.Close(new Complex(3, 4), d, 1e-5);
        Assert.True(v[1].Abs < 1e-6);
        string[] names = { "Arnoldi", "Bidiagonal", "Cholesky", "Diagonal", "EVD", "GEVD", "GramSchmidt", "GSVD", "Hessenberg", "Householder", "Lanczos", "LDL", "LDU", "LQ", "LU", "NMF", "Polar", "Power", "QL", "QR", "QZ", "RQ", "Schur", "SVD", "UDL" };
        foreach (string name in names)
        {
            var type = typeof(SVD).Assembly.GetType("UMapx.Decomposition." + name)!;
            Assert.True(type.IsAbstract && type.IsSealed);
            Assert.DoesNotContain(type.GetFields(System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.NonPublic), f => !f.IsLiteral && !f.IsInitOnly);
            Assert.Empty(type.GetProperties());
        }
    }

    [Fact]
    public void InvalidInputsAndIterationLimitsAreRejected()
    {
        var bad = new Complex32[,] { { new(float.NaN, 0) } };
        Assert.Throws<ArgumentNullException>(() => QR.Decompose((Complex32[,])null!));
        Assert.Throws<ArgumentException>(() => SVD.Decompose(bad));
        Assert.Throws<ArgumentException>(() => Schur.Decompose(new Complex32[0, 0]));
        Assert.Throws<ArgumentException>(() => QZ.Decompose(new Complex32[2, 2], new Complex32[3, 3]));
        Assert.Throws<ArgumentException>(() => Cholesky.Decompose(new Complex32[,] { { 1, new(0, 1) }, { 0, 1 } }));
        Assert.Throws<ArgumentException>(() => Cholesky.Decompose(new Complex32[,] { { -1 } }));
        Assert.Throws<ArgumentOutOfRangeException>(() => SVD.Decompose(new Complex32[1, 1], 0));
        Assert.Throws<InvalidOperationException>(() => SVD.Decompose(Sample(8, 8, "Dense", 1), 1));
        Assert.Throws<InvalidOperationException>(() => Schur.Decompose(Sample(8, 8, "Dense", 1), iterations: 1));
        Assert.Throws<InvalidOperationException>(() => QZ.Decompose(Sample(8, 8, "Dense", 1), Sample(8, 8, "Dense", 2, 31), iterations: 1));
        Assert.Throws<ArgumentException>(() => NMF.Decompose(new float[,] { { -1 } }, 1));
    }

    [Theory]
    [InlineData(2)] [InlineData(5)] [InlineData(12)] [InlineData(24)]
    public void IndependentDensePencilsAndKnownSpectraSurviveLargerProblems(int n)
    {
        for (int seed = 0; seed < 3; seed++)
        {
            var a = Sample(n, n, "Dense", 1, seed * 137);
            var b = Sample(n, n, "Dense", 1, seed * 137 + 61);
            var qz = QZ.Decompose(a, b);
            Relative(Work(a), Product(Product(Work(qz.Q), Work(qz.S)), Adjoint(Work(qz.Z))));
            Relative(Work(b), Product(Product(Work(qz.Q), Work(qz.T)), Adjoint(Work(qz.Z))));
            Orthonormal(qz.Q); Orthonormal(qz.Z);
            var q = Work(QR.Decompose(a).Q);
            var spectrum = Enumerable.Range(1, n).Select(i => new Complex32(i, i % 3 - 1)).ToArray();
            var known = Single(Product(Product(q, Diagonal(spectrum)), Adjoint(q)));
            var evd = EVD.Decompose(known);
            var values = evd.D.OrderBy(z => z.Real).ToArray();
            for (int i = 0; i < n; i++) NumericAssert.Close((Complex)spectrum[i], values[i], 2e-4);
        }
    }

    [Theory]
    [InlineData(1)] [InlineData(3)] [InlineData(6)]
    public void RealQzRetainsAnOrthogonalLeftFactorWhenBIsSingular(int n)
    {
        foreach (bool zero in new[] { false, true })
        {
            var a = NumericAssert.Matrix(n, n);
            var b = new float[n, n];
            if (!zero) for (int i = 0; i < n; i++) b[i, i] = i;
            var d = QZ.Decompose(a, b);
            NumericAssert.Close(a, NumericAssert.Product(NumericAssert.Product(d.Q, d.S), NumericAssert.Transpose(d.Z)), 1e-4f);
            NumericAssert.Close(b, NumericAssert.Product(NumericAssert.Product(d.Q, d.T), NumericAssert.Transpose(d.Z)), 1e-4f);
            NumericAssert.Close(NumericAssert.Diagonal(Enumerable.Repeat(1f, n).ToArray()), NumericAssert.Product(NumericAssert.Transpose(d.Q), d.Q), 1e-4f);
        }
    }

    [Fact]
    public void IndependentCallsDoNotShareWorkspaceOrResultArrays()
    {
        var a = Sample(7, 5, "Dense", 1);
        var original = Work(a);
        Parallel.For(0, 16, i =>
        {
            var d = SVD.Decompose(a);
            Relative(original, Product(Product(Work(d.U), Diagonal(d.S)), Adjoint(Work(d.V))));
            d.U[0, 0] = 123;
        });
        Relative(original, Work(a));
        var tall = QR.Decompose(Sample(2000, 2, "Dense", 1));
        Assert.Equal(2, tall.Q.GetLength(1));
        Orthonormal(tall.Q);
    }

    private static Complex32[,] Sample(int m, int n, string kind, float scale, int seed = 0)
    {
        var random = new Random(809 + 13 * m + 17 * n + seed);
        var a = new Complex32[m, n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                Complex z = kind switch
                {
                    "Zero" => Complex.Zero,
                    "Diagonal" => i == j ? new Complex(i + 1, 1 - i) : Complex.Zero,
                    "Jordan" => i == j ? new Complex(2, 1) : j == i + 1 ? Complex.One : Complex.Zero,
                    "RankOne" => new Complex(i + 1, 1) * new Complex(j + 1, -1),
                    _ => new Complex(random.NextDouble() * 2 - 1, random.NextDouble() * 2 - 1)
                };
                a[i, j] = new Complex32((float)(z.Real * scale), (float)(z.Imaginary * scale));
            }
        return a;
    }

    private static Complex[,] Work(float[,] a)
    {
        var b = new Complex[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < b.GetLength(0); i++) for (int j = 0; j < b.GetLength(1); j++) b[i, j] = a[i, j];
        return b;
    }

    private static Complex[,] Work(Complex32[,] a)
    {
        var b = new Complex[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < b.GetLength(0); i++) for (int j = 0; j < b.GetLength(1); j++) b[i, j] = (Complex)a[i, j];
        return b;
    }

    private static Complex32[,] Single(Complex[,] a)
    {
        var b = new Complex32[a.GetLength(0), a.GetLength(1)];
        for (int i = 0; i < b.GetLength(0); i++) for (int j = 0; j < b.GetLength(1); j++) b[i, j] = new((float)a[i, j].Real, (float)a[i, j].Imaginary);
        return b;
    }

    private static Complex[,] Product(Complex[,] a, Complex[,] b)
    {
        Assert.Equal(a.GetLength(1), b.GetLength(0));
        var c = new Complex[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < c.GetLength(0); i++) for (int j = 0; j < c.GetLength(1); j++)
            for (int k = 0; k < a.GetLength(1); k++) c[i, j] += a[i, k] * b[k, j];
        return c;
    }

    private static Complex[,] Adjoint(Complex[,] a)
    {
        var b = new Complex[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) b[j, i] = Complex.Conjugate(a[i, j]);
        return b;
    }

    private static Complex[,] Diagonal(float[] d) => Diagonal(d.Select(x => new Complex32(x, 0)).ToArray());
    private static Complex[,] Diagonal(Complex32[] d)
    {
        var a = new Complex[d.Length, d.Length];
        for (int i = 0; i < d.Length; i++) a[i, i] = (Complex)d[i];
        return a;
    }

    private static Complex[,] Rows(Complex[,] a, int[] permutation)
    {
        var b = new Complex[a.GetLength(0), a.GetLength(1)];
        Assert.Equal(Enumerable.Range(0, permutation.Length), permutation.OrderBy(i => i));
        for (int i = 0; i < b.GetLength(0); i++) for (int j = 0; j < b.GetLength(1); j++) b[i, j] = a[permutation[i], j];
        return b;
    }

    private static double Norm(Complex[,] a) => Math.Sqrt(a.Cast<Complex>().Sum(z => z.Magnitude * z.Magnitude));

    private static void Relative(Complex[,] expected, Complex[,] actual, double tolerance = 3e-5)
    {
        Assert.Equal(expected.GetLength(0), actual.GetLength(0)); Assert.Equal(expected.GetLength(1), actual.GetLength(1));
        double error = 0;
        for (int i = 0; i < expected.GetLength(0); i++) for (int j = 0; j < expected.GetLength(1); j++)
        {
            Assert.True(double.IsFinite(actual[i, j].Real) && double.IsFinite(actual[i, j].Imaginary));
            double r = Complex.Abs(expected[i, j] - actual[i, j]); error += r * r;
        }
        double scale = Math.Max(Norm(expected), Norm(actual));
        Assert.True(Math.Sqrt(error) <= tolerance * scale + 1e-290, $"Relative error {Math.Sqrt(error) / scale} exceeds {tolerance}.");
    }

    private static void Orthonormal(Complex32[,] a, bool columns = true)
    {
        var w = Work(a); int n = a.GetLength(columns ? 1 : 0);
        Relative(Diagonal(Enumerable.Repeat(1f, n).ToArray()), columns ? Product(Adjoint(w), w) : Product(w, Adjoint(w)));
    }

    private static void Band(Complex32[,] a, int lower, int upper)
    {
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++)
            if (i - j > lower || j - i > upper) Assert.Equal(new Complex32(0, 0), a[i, j]);
    }
}
