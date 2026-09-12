using System;
using UMapx.Core;
using C = System.Numerics.Complex;

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
            MatrixMath.Copy(matrix);
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
            var w = new double[m, rank];
            var h = new double[rank, n];
            var random = new Random(1729);
            for (int i = 0; i < m; i++) for (int k = 0; k < rank; k++) w[i, k] = 0.5 + random.NextDouble();
            for (int k = 0; k < rank; k++) for (int j = 0; j < n; j++) h[k, j] = 0.5 + random.NextDouble();
            for (int step = 0; step < iterations; step++)
            {
                var gram = new double[rank, rank];
                for (int k = 0; k < rank; k++)
                    for (int l = 0; l < rank; l++)
                        for (int i = 0; i < m; i++) gram[k, l] += w[i, k] * w[i, l];
                var nextH = new double[rank, n];
                for (int k = 0; k < rank; k++)
                    for (int j = 0; j < n; j++)
                    {
                        double numerator = 0, denominator = 0;
                        for (int i = 0; i < m; i++) numerator += w[i, k] * (matrix[i, j] / scale);
                        for (int l = 0; l < rank; l++) denominator += gram[k, l] * h[l, j];
                        nextH[k, j] = denominator == 0 ? 0 : h[k, j] * numerator / denominator;
                    }
                h = nextH;
                gram = new double[rank, rank];
                for (int k = 0; k < rank; k++)
                    for (int l = 0; l < rank; l++)
                        for (int j = 0; j < n; j++) gram[k, l] += h[k, j] * h[l, j];
                var nextW = new double[m, rank];
                for (int i = 0; i < m; i++)
                    for (int k = 0; k < rank; k++)
                    {
                        double numerator = 0, denominator = 0;
                        for (int j = 0; j < n; j++) numerator += (matrix[i, j] / scale) * h[k, j];
                        for (int l = 0; l < rank; l++) denominator += w[i, l] * gram[l, k];
                        nextW[i, k] = denominator == 0 ? 0 : w[i, k] * numerator / denominator;
                    }
                w = nextW;
            }
            double factorScale = Math.Sqrt(scale);
            for (int i = 0; i < m; i++) for (int k = 0; k < rank; k++) resultW[i, k] = (float)(w[i, k] * factorScale);
            for (int k = 0; k < rank; k++) for (int j = 0; j < n; j++) resultH[k, j] = (float)(h[k, j] * factorScale);
            return (resultW, resultH);
        }
    }
}
