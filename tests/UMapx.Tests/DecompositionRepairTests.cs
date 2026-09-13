using System.Numerics;
using UMapx.Decomposition;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class DecompositionRepairTests
{
    public static IEnumerable<object[]> SingularCases()
    {
        foreach (var shape in new[] { (1, 1), (1, 7), (7, 1), (3, 5), (5, 3), (4, 4), (8, 8), (9, 13), (13, 9) })
            foreach (string kind in new[] { "Zero", "RankOne", "RankTwo", "Diagonal", "Dense" })
                foreach (float scale in new[] { 1e-30f, 1f, 1e30f })
                    yield return new object[] { shape.Item1, shape.Item2, kind, scale };
    }

    [Theory, MemberData(nameof(SingularCases))]
    public void SingularFactorsAndPseudoinverseRespectRankAndScale(int rows, int columns, string kind, float scale)
    {
        var input = Sample(rows, columns, kind, scale);
        var original = (float[,])input.Clone();
        var decomposition = SVD.Decompose(input);
        var u = Double(decomposition.U);
        var v = Double(decomposition.V);
        var singular = decomposition.S;
        Assert.Equal(Math.Min(rows, columns), singular.Length);
        Assert.Equal(rows, u.GetLength(0));
        Assert.Equal(columns, v.GetLength(0));
        var diagonal = new double[singular.Length, singular.Length];
        for (int i = 0; i < singular.Length; i++)
        {
            Assert.True(float.IsFinite(singular[i]) && singular[i] >= 0);
            if (i > 0) Assert.True(singular[i - 1] >= singular[i]);
            diagonal[i, i] = singular[i];
        }
        Orthonormal(u);
        Orthonormal(v);
        Relative(Double(input), Product(Product(u, diagonal), Transpose(v)), 2e-5);
        var inverse = Double(SVD.PseudoInverse(decomposition.U, decomposition.S, decomposition.V));
        var a = Double(input);
        var ap = Product(a, inverse);
        var pa = Product(inverse, a);
        Relative(a, Product(ap, a), 5e-5);
        Relative(inverse, Product(pa, inverse), 5e-5);
        Relative(ap, Transpose(ap), 5e-5);
        Relative(pa, Transpose(pa), 5e-5);

        if (kind == "RankOne")
        {
            // A = scale*x*y^T => A+ = y*x^T / (scale*(x^T*x)*(y^T*y)).
            double left = Enumerable.Range(1, rows).Sum(i => (double)i * i);
            double right = Enumerable.Range(1, columns).Sum(i => (double)i * i);
            var expected = new double[columns, rows];
            for (int i = 0; i < columns; i++)
                for (int j = 0; j < rows; j++) expected[i, j] = (i + 1.0) * (j + 1) / (scale * left * right);
            Relative(expected, inverse, 2e-5);
        }
        Assert.Equal(original.Cast<float>(), input.Cast<float>());
    }

    [Theory]
    [InlineData(1e-30f)] [InlineData(1f)] [InlineData(1e30f)]
    public void PseudoinverseThresholdUsesRelativeSinglePrecisionRank(float scale)
    {
        var input = new float[,] { { scale, 0, 0 }, { 0, scale * 1e-4f, 0 }, { 0, 0, scale * 1e-8f } };
        var factors = SVD.Decompose(input);
        var inverse = SVD.PseudoInverse(factors.U, factors.S, factors.V);
        NumericAssert.Close(1.0 / input[0, 0], inverse[0, 0], 0, 2e-6);
        NumericAssert.Close(1.0 / input[1, 1], inverse[1, 1], 0, 2e-6);
        Assert.Equal(0, inverse[2, 2]);
    }

    [Theory] [InlineData(0)] [InlineData(-1)]
    public void SingularDecompositionRejectsAnInvalidIterationLimit(int iterations)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SVD.Decompose(new float[,] { { 1 } }, iterations));
    }

    [Fact]
    public void SingularDecompositionReportsNonconvergenceInsteadOfReturningPartialFactors()
    {
        Assert.Throws<InvalidOperationException>(() => SVD.Decompose(Sample(8, 8, "Dense", 1), 1));
    }

    public static IEnumerable<object[]> SchurCases()
    {
        foreach (int size in new[] { 1, 2, 3, 5, 8 })
            foreach (string kind in new[] { "Zero", "RankOne", "Diagonal", "Dense", "Jordan", "Rotation" })
                foreach (float scale in new[] { 1e-30f, 1f, 1e30f })
                    foreach (float tolerance in new[] { 0f, 1e-7f, 1e-16f })
                        yield return new object[] { size, kind, scale, tolerance };
    }

    [Theory, MemberData(nameof(SchurCases))]
    public void SchurHandlesDeflationRepeatedRootsAndScaling(int size, string kind, float scale, float tolerance)
    {
        var input = Sample(size, size, kind, scale);
        var original = (float[,])input.Clone();
        var decomposition = Schur.Decompose(input, tolerance);
        var q = Double(decomposition.Q);
        var t = Double(decomposition.T);
        Orthonormal(q);
        Relative(Double(input), Product(Product(q, t), Transpose(q)), 2e-5);
        for (int i = 0; i < size; i++)
            for (int j = 0; j + 1 < i; j++) Assert.Equal(0, t[i, j]);
        // A real Schur form has isolated blocks of size at most two.
        for (int i = 2; i < size; i++) Assert.True(t[i, i - 1] == 0 || t[i - 1, i - 2] == 0);
        Assert.Equal(original.Cast<float>(), input.Cast<float>());
    }

    public static IEnumerable<object[]> GeneralizedCases()
    {
        foreach (int size in new[] { 1, 2, 3, 5, 8 })
            foreach (string kind in new[] { "Zero", "Diagonal", "Repeated", "Rotation", "Dense" })
                foreach (float scale in new[] { 1e-20f, 1f, 1e20f })
                    yield return new object[] { size, kind, scale };
    }

    [Theory, MemberData(nameof(GeneralizedCases))]
    public void GeneralizedEigenvectorsSatisfyThePencilAndNormalization(int size, string kind, float scale)
    {
        var c = Sample(size, size, kind == "Repeated" ? "Diagonal" : kind, 1);
        if (kind == "Repeated") for (int i = 0; i < size; i++) c[i, i] = 2;
        var random = Double(Sample(size, size, "Dense", 1));
        var bBase = Product(Transpose(random), random);
        for (int i = 0; i < size; i++) bBase[i, i] += 1;
        var a = Single(Product(bBase, Double(c)), scale);
        var b = Single(bBase, scale);
        var originalA = (float[,])a.Clone();
        var originalB = (float[,])b.Clone();
        var decomposition = GEVD.Decompose(a, b);
        Assert.False(GEVD.IsSingular(decomposition.Beta));
        var vectors = Double(decomposition.V);
        var eigenvalues = GEVD.Eigenvalues(decomposition.Alpha, decomposition.Beta);
        Relative(Product(Double(a), vectors), Product(Product(Double(b), vectors), Double(GEVD.RealEigenvalueMatrix(decomposition.Alpha, decomposition.Beta))), 3e-5);
        for (int j = 0; j < size; j++)
        {
            bool complex = decomposition.Alpha[j].Imag > 0;
            double max = 0;
            for (int i = 0; i < size; i++)
                max = Math.Max(max, complex ? Complex.Abs(new Complex(vectors[i, j], vectors[i, j + 1])) : Math.Abs(vectors[i, j]));
            NumericAssert.Close(1, max, 2e-5, 0);
            Assert.True(float.IsFinite(eigenvalues[j].Real) && float.IsFinite(eigenvalues[j].Imag));
            if (complex) j++;
        }
        if (kind is "Zero" or "Diagonal" or "Repeated")
        {
            double[] expected = Enumerable.Range(0, size).Select(i => (double)c[i, i]).OrderBy(x => x).ToArray();
            double[] actual = eigenvalues.Select(z => (double)z.Real).OrderBy(x => x).ToArray();
            for (int i = 0; i < size; i++) NumericAssert.Close(expected[i], actual[i], 2e-5, 2e-5);
        }
        Assert.Equal(originalA.Cast<float>(), a.Cast<float>());
        Assert.Equal(originalB.Cast<float>(), b.Cast<float>());
    }

    [Fact]
    public void GeneralizedInfiniteEigenvalueRetainsAHomogeneousEigenvector()
    {
        var a = new float[,] { { 2, 0, 0 }, { 0, 3, 0 }, { 0, 0, 4 } };
        var b = new float[,] { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 2 } };
        var decomposition = GEVD.Decompose(a, b);
        Assert.True(GEVD.IsSingular(decomposition.Beta));
        Assert.Equal(1, decomposition.Beta.Count(x => x == 0));
        var v = decomposition.V;
        for (int j = 0; j < 3; j++)
        {
            Assert.Contains(Enumerable.Range(0, 3), i => v[i, j] != 0);
            for (int i = 0; i < 3; i++)
                NumericAssert.Close(decomposition.Beta[j] * a[i, i] * v[i, j], decomposition.Alpha[j].Real * b[i, i] * v[i, j]);
        }
        Assert.All(v.Cast<float>(), x => Assert.True(float.IsFinite(x)));
    }

    [Fact]
    public void ScalingDoesNotEraseRepresentableIsolatedValues()
    {
        var diagonal = new float[,] { { 1e30f, 0 }, { 0, 1e-30f } };
        var singular = SVD.Decompose(diagonal).S;
        NumericAssert.Close(diagonal[0, 0], singular[0], 0, 2e-6);
        NumericAssert.Close(diagonal[1, 1], singular[1], 0, 2e-6);
        var schur = Schur.Decompose(diagonal).T;
        NumericAssert.Close(diagonal[0, 0], schur[0, 0], 0, 2e-6);
        NumericAssert.Close(diagonal[1, 1], schur[1, 1], 0, 2e-6);
    }

    [Fact]
    public void GeneralizedHomogeneousValuesSurviveDifferentInputUnits()
    {
        var a = new float[,] { { 1e30f, 0 }, { 0, 2e30f } };
        var b = new float[,] { { 1e-30f, 0 }, { 0, 2e-30f } };
        var decomposition = GEVD.Decompose(a, b);
        Assert.False(GEVD.IsSingular(decomposition.Beta));
        for (int j = 0; j < 2; j++)
        {
            Assert.True(decomposition.Beta[j] > 0);
            double quotient = (double)decomposition.Alpha[j].Real / decomposition.Beta[j];
            NumericAssert.Close((double)a[j, j] / b[j, j], quotient, 0, 2e-6);
        }
    }

    [Theory] [InlineData(float.NaN)] [InlineData(float.PositiveInfinity)] [InlineData(float.NegativeInfinity)]
    public void DecompositionsRejectNonfiniteMatrixEntries(float value)
    {
        var invalid = new float[,] { { value } };
        var valid = new float[,] { { 1 } };
        Assert.Throws<ArgumentException>(() => SVD.Decompose(invalid));
        Assert.Throws<ArgumentException>(() => Schur.Decompose(invalid));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(invalid, valid));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(valid, invalid));
    }

    [Fact]
    public void DecompositionsRejectUndefinedDimensionsAndTolerances()
    {
        Assert.Throws<ArgumentNullException>(() => SVD.Decompose((float[,])null!));
        Assert.Throws<ArgumentNullException>(() => Schur.Decompose((float[,])null!));
        Assert.Throws<ArgumentNullException>(() => GEVD.Decompose(null!, new float[1, 1]));
        Assert.Throws<ArgumentNullException>(() => GEVD.Decompose(new float[1, 1], null!));
        Assert.Throws<ArgumentException>(() => SVD.Decompose(new float[0, 2]));
        Assert.Throws<ArgumentException>(() => SVD.Decompose(new float[2, 0]));
        Assert.Throws<ArgumentException>(() => Schur.Decompose(new float[0, 0]));
        Assert.Throws<ArgumentException>(() => Schur.Decompose(new float[2, 3]));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(new float[0, 0], new float[0, 0]));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(new float[2, 2], new float[3, 3]));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(new float[2, 3], new float[2, 3]));
        Assert.Throws<ArgumentException>(() => GEVD.Decompose(new float[2, 2], new float[2, 3]));
        Assert.Throws<ArgumentOutOfRangeException>(() => Schur.Decompose(new float[1, 1], float.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => GEVD.Decompose(new float[1, 1], new float[1, 1], float.NaN));
    }

    [Theory] [InlineData(2)] [InlineData(3)] [InlineData(5)]
    public void QzKeepsItsTriangularFactorFreeOfScratchStorage(int size)
    {
        var a = Sample(size, size, "Dense", 1);
        var b = NumericAssert.Product(NumericAssert.Transpose(a), a);
        for (int i = 0; i < size; i++) b[i, i] += 1;
        var decomposition = QZ.Decompose(a, b, 1e-7f);
        var q = Double(decomposition.Q);
        var z = Double(decomposition.Z);
        Orthonormal(q);
        Orthonormal(z);
        Relative(Double(a), Product(Product(q, Double(decomposition.S)), Transpose(z)), 2e-5);
        Relative(Double(b), Product(Product(q, Double(decomposition.T)), Transpose(z)), 2e-5);
        for (int i = 1; i < size; i++)
            for (int j = 0; j < i; j++) Assert.Equal(0, decomposition.T[i, j]);
    }

    private static float[,] Sample(int rows, int columns, string kind, float scale)
    {
        var random = new Random(719 + 13 * rows + columns);
        var result = new float[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
            {
                double value = kind switch
                {
                    "Zero" => 0,
                    "RankOne" => (i + 1.0) * (j + 1),
                    "RankTwo" => (i + 1.0) * (j + 1) + (i % 2 == 0 ? 1 : -1) * (j % 3 - 1),
                    "Diagonal" => i == j ? 1 + (i % 3) : 0,
                    "Jordan" => j == i + 1 ? 1 : 0,
                    "Rotation" => i == j ? .25 : i % 2 == 0 && j == i + 1 ? -2 : j % 2 == 0 && i == j + 1 ? 2 : 0,
                    _ => 2 * random.NextDouble() - 1
                };
                result[i, j] = (float)(value * scale);
            }
        return result;
    }

    // Independent double arithmetic keeps large/small scale checks meaningful.
    private static double[,] Double(float[,] source)
    {
        var result = new double[source.GetLength(0), source.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++) result[i, j] = source[i, j];
        return result;
    }

    private static float[,] Single(double[,] source, float scale)
    {
        var result = new float[source.GetLength(0), source.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++) result[i, j] = (float)(source[i, j] * scale);
        return result;
    }

    private static double[,] Product(double[,] a, double[,] b)
    {
        Assert.Equal(a.GetLength(1), b.GetLength(0));
        var result = new double[a.GetLength(0), b.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
            for (int j = 0; j < result.GetLength(1); j++)
                for (int k = 0; k < a.GetLength(1); k++) result[i, j] += a[i, k] * b[k, j];
        return result;
    }

    private static double[,] Transpose(double[,] a)
    {
        var result = new double[a.GetLength(1), a.GetLength(0)];
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++) result[j, i] = a[i, j];
        return result;
    }

    private static void Orthonormal(double[,] matrix)
    {
        var identity = new double[matrix.GetLength(1), matrix.GetLength(1)];
        for (int i = 0; i < identity.GetLength(0); i++) identity[i, i] = 1;
        Relative(identity, Product(Transpose(matrix), matrix), 2e-5);
    }

    private static void Relative(double[,] expected, double[,] actual, double tolerance)
    {
        Assert.Equal(expected.GetLength(0), actual.GetLength(0));
        Assert.Equal(expected.GetLength(1), actual.GetLength(1));
        double norm = 0, error = 0;
        for (int i = 0; i < expected.GetLength(0); i++)
            for (int j = 0; j < expected.GetLength(1); j++)
            {
                Assert.True(double.IsFinite(actual[i, j]), $"Nonfinite entry ({i}, {j}).");
                double difference = expected[i, j] - actual[i, j];
                norm += expected[i, j] * expected[i, j];
                error += difference * difference;
            }
        Assert.True(Math.Sqrt(error) <= tolerance * Math.Sqrt(norm),
            $"Relative Frobenius residual {Math.Sqrt(error / norm):G6} exceeds {tolerance:G6}.");
    }
}
