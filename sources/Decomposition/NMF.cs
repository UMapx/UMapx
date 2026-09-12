using System;

namespace UMapx.Decomposition
{
    /// <summary>Provides nonnegative factorization of real matrices; nonnegativity is not defined for general complex numbers.</summary>
    public static class NMF
    {
        /// <summary>Approximates a nonnegative matrix as A = W H using multiplicative least-squares updates.</summary>
        /// <param name="matrix">Finite nonempty real matrix with nonnegative entries.</param>
        /// <param name="rank">Positive number of factor columns/rows, at most min(m,n).</param>
        /// <param name="iterations">Positive number of alternating updates.</param>
        /// <returns>Nonnegative W of size m by rank and H of size rank by n.</returns>
        public static (float[,] W, float[,] H) Decompose(float[,] matrix, int rank, int iterations = 100)
        {
            InternalMatrixMath.Validate(matrix);
            int m = matrix.GetLength(0), n = matrix.GetLength(1);
            if (rank < 1 || rank > Math.Min(m, n)) throw new ArgumentOutOfRangeException(nameof(rank));
            if (iterations < 1) throw new ArgumentOutOfRangeException(nameof(iterations));
            double scale = 0;
            foreach (float value in matrix)
            {
                if (value < 0) throw new ArgumentException("NMF requires nonnegative values.", nameof(matrix));
                scale = Math.Max(scale, value);
            }
            var resultW = new float[m, rank];
            var resultH = new float[rank, n];
            if (scale == 0) return (resultW, resultH);
            var a = InternalMatrixMath.CopyJagged(matrix);
            InternalMatrixMath.Divide(a, scale);
            var w = InternalMatrixMath.CreateJagged(m, rank);
            var h = InternalMatrixMath.CreateJagged(rank, n);
            var nextW = InternalMatrixMath.CreateJagged(m, rank);
            var nextH = InternalMatrixMath.CreateJagged(rank, n);
            var numeratorH = InternalMatrixMath.CreateJagged(rank, n);
            var gram = InternalMatrixMath.CreateJagged(rank, rank);
            var random = new Random(1729);
            for (int i = 0; i < m; i++) for (int k = 0; k < rank; k++) w[i][k] = 0.5 + random.NextDouble();
            for (int k = 0; k < rank; k++) for (int j = 0; j < n; j++) h[k][j] = 0.5 + random.NextDouble();
            for (int step = 0; step < iterations; step++)
            {
                for (int k = 0; k < rank; k++)
                {
                    Array.Clear(gram[k], 0, rank);
                    Array.Clear(numeratorH[k], 0, n);
                }
                // W^T W and W^T A: keep the wide inner loop contiguous in memory.
                for (int i = 0; i < m; i++)
                    for (int k = 0; k < rank; k++)
                    {
                        double weight = w[i][k];
                        for (int l = 0; l < rank; l++) gram[k][l] += weight * w[i][l];
                        var numerator = numeratorH[k]; var row = a[i];
                        for (int j = 0; j < n; j++) numerator[j] += weight * row[j];
                    }
                for (int k = 0; k < rank; k++)
                    for (int j = 0; j < n; j++)
                    {
                        double denominator = 0;
                        for (int l = 0; l < rank; l++) denominator += gram[k][l] * h[l][j];
                        nextH[k][j] = denominator == 0 ? 0 : h[k][j] * numeratorH[k][j] / denominator;
                    }
                var oldH = h; h = nextH; nextH = oldH;
                for (int k = 0; k < rank; k++)
                    for (int l = 0; l <= k; l++)
                        gram[k][l] = gram[l][k] = InternalMatrixMath.Dot(h[k], h[l], n);
                for (int i = 0; i < m; i++)
                    for (int k = 0; k < rank; k++)
                    {
                        double numerator = InternalMatrixMath.Dot(a[i], h[k], n);
                        double denominator = InternalMatrixMath.Dot(w[i], gram[k], rank);
                        nextW[i][k] = denominator == 0 ? 0 : w[i][k] * numerator / denominator;
                    }
                var oldW = w; w = nextW; nextW = oldW;
            }
            double factorScale = Math.Sqrt(scale);
            for (int i = 0; i < m; i++) for (int k = 0; k < rank; k++) resultW[i, k] = (float)(w[i][k] * factorScale);
            for (int k = 0; k < rank; k++) for (int j = 0; j < n; j++) resultH[k, j] = (float)(h[k][j] * factorScale);
            return (resultW, resultH);
        }
    }
}
