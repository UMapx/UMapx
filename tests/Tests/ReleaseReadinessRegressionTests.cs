using UMapx.Decomposition;
using UMapx.Core;
using Xunit;

namespace UMapx.Tests;

[Trait("Category", "Decomposition")]
public class ReleaseReadinessRegressionTests
{
    public static IEnumerable<object[]> IndependentBlocks()
    {
        foreach (var scales in new[] { (1e20f, 1f), (1f, 1e-20f), (1e30f, 1e-30f) })
        foreach (int order in new[] { 0, 1, 2 })
        foreach (bool singular in new[] { false, true })
            yield return new object[] { scales.Item1, scales.Item2, order, singular };
    }

    [Theory, MemberData(nameof(IndependentBlocks))]
    public void IndependentBlocksRetainTheirSpectraAndEigenvectors(float large, float small, int order, bool singular)
    {
        int[] permutation = order == 0 ? new[] { 0, 1, 2, 3 }
            : order == 1 ? new[] { 2, 3, 0, 1 } : new[] { 0, 2, 1, 3 };
        var a = new float[4, 4];
        var block = singular ? new float[,] { { 1, 2 }, { 2, 4 } } : new float[,] { { 2, 1 }, { 1, 2 } };
        for (int i = 0; i < 2; i++) for (int j = 0; j < 2; j++)
        {
            a[permutation[i], permutation[j]] = large * (i == j ? 2 : 1);
            a[permutation[2 + i], permutation[2 + j]] = small * block[i, j];
        }
        var eigen = EVD.Decompose(a);
        double[] expected = { singular ? 0 : small, (singular ? 5 : 3) * (double)small, large, 3 * (double)large };
        for (int j = 0; j < 4; j++)
        {
            NumericAssert.Close(expected[j], eigen.D[j].Real, expected[j] == 0 ? small * 2e-6 : 0, 2e-6);
            if (j >= 2) continue;
            double norm = 0;
            for (int i = 0; i < 4; i++) norm += (double)eigen.V[i, j] * eigen.V[i, j];
            NumericAssert.Close(1, norm);
            for (int i = 2; i < 4; i++)
            {
                int row = permutation[i];
                double av = 0;
                for (int k = 2; k < 4; k++) av += (double)a[row, permutation[k]] * eigen.V[permutation[k], j];
                NumericAssert.Close(av, (double)eigen.V[row, j] * eigen.D[j].Real, small * 2e-6, 2e-6);
            }
        }
        foreach (int shape in new[] { 0, 1, 2 })
        {
            var input = new float[shape == 1 ? 6 : 4, shape == 2 ? 6 : 4];
            for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) input[i, j] = a[i, j];
            var d = SVD.Decompose(input);
            for (int j = 0; j < 4; j++)
                NumericAssert.Close(expected[3 - j], d.S[j], expected[3 - j] == 0 ? small * 2e-6 : 0, 2e-6);
            foreach (var q in new[] { d.U, d.V })
                for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++)
                {
                    double dot = 0;
                    for (int k = 0; k < q.GetLength(0); k++) dot += (double)q[k, i] * q[k, j];
                    NumericAssert.Close(i == j ? 1 : 0, dot);
                }
            for (int i = 2; i < 4; i++) for (int j = 2; j < 4; j++)
            {
                int row = permutation[i], column = permutation[j];
                double value = 0;
                for (int k = 0; k < 4; k++) value += (double)d.U[row, k] * d.S[k] * d.V[column, k];
                NumericAssert.Close(input[row, column], value, small * 2e-6, 2e-6);
            }
        }
    }

    [Theory]
    [InlineData(1e-10f)]
    [InlineData(1e-20f)]
    [InlineData(1e-30f)]
    public void QrPreservesSmallLeadingComponents(float small)
    {
        var a = new float[,] { { small }, { 1 } };
        var real = QR.Decompose(a);
        NumericAssert.Close(small, (double)real.Q[0, 0] * real.R[0, 0], 0, 2e-6);
        var complex = QR.Decompose(new Complex32[,] { { new(small, small) }, { 1 } });
        NumericAssert.Close(new System.Numerics.Complex(small, small),
            complex.Q[0, 0] * complex.R[0, 0], 0, 2e-6);
    }

    [Theory]
    [InlineData(1e20f, 1f)]
    [InlineData(1f, 1e-20f)]
    [InlineData(1e30f, 1e-30f)]
    public void RealSymmetricEigenvaluesPreserveCoupledSmallBlocks(float large, float small)
    {
        var a = new float[,] { { large, 0, 0 }, { 0, 2 * small, small }, { 0, small, 2 * small } };
        var d = EVD.Decompose(a);
        var values = d.D.OrderBy(value => value.Real).ToArray();
        NumericAssert.Close(small, values[0].Real, 0, 2e-6);
        NumericAssert.Close(3 * (double)small, values[1].Real, 0, 2e-6);
        NumericAssert.Close(large, values[2].Real, 0, 2e-6);
        for (int j = 0; j < d.D.Length; j++)
        {
            Assert.Equal(0, d.D[j].Imag);
            if (d.D[j].Real > 4 * small) continue;
            for (int i = 1; i < 3; i++)
            {
                double product = 0;
                for (int k = 1; k < 3; k++) product += (double)a[i, k] * d.V[k, j];
                NumericAssert.Close(product, (double)d.V[i, j] * d.D[j].Real, small * 2e-6, 2e-6);
            }
        }
    }

    [Theory]
    [InlineData(1e20f, 1f, false)]
    [InlineData(1e20f, 1f, true)]
    [InlineData(1f, 1e-20f, false)]
    [InlineData(1f, 1e-20f, true)]
    [InlineData(1e30f, 1e-30f, false)]
    public void RealSvdPreservesCoupledSmallBlocks(float large, float small, bool transpose)
    {
        // The independent 2x2 block has exact singular values 3*small and small.
        // A global residual dominated by 'large' cannot detect errors in this block.
        var a = new float[,] { { large, 0, 0 }, { 0, 2 * small, small }, { 0, small, 2 * small }, { 0, 0, 0 } };
        if (transpose) a = NumericAssert.Transpose(a);
        var d = SVD.Decompose(a, 100);
        NumericAssert.Close(large, d.S[0], 0, 2e-6);
        NumericAssert.Close(3 * (double)small, d.S[1], 0, 2e-6);
        NumericAssert.Close(small, d.S[2], 0, 2e-6);
        for (int i = 1; i < 3; i++) for (int j = 1; j < 3; j++)
        {
            double value = 0;
            for (int k = 0; k < d.S.Length; k++) value += (double)d.U[i, k] * d.S[k] * d.V[j, k];
            NumericAssert.Close(a[i, j], value, 0, 2e-6);
        }
    }
}
