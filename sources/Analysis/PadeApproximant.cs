using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a Pade approximant.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Pad%C3%A9_approximant
    /// </remarks>
    /// Example: exp(x) = 1 + x + x^2/2! + x^3/3! + ...
    /// <code>
    /// float[] taylorExp = new float[] { 1.0f, 1.0f, 0.5f, 1.0f/6.0f, 1.0f/24.0f, 1.0f/120.0f, 1.0f/720.0f };
    /// </code>
    /// </summary>
    [Serializable]
    public class PadeApproximant
    {
        #region Private data
        private readonly float[] numeratorCoeffs;
        private readonly float[] denominatorCoeffs;
        #endregion

        /// <summary>
        /// Initializes the Pade approximant.
        /// </summary>
        /// <param name="taylorCoeffs">Taylor series coefficients</param>
        /// <param name="m">The degree of the numerator of a rational function</param>
        /// <param name="n">The degree of the denominator of a rational function</param>
        /// <exception cref="ArgumentException"></exception>
        public PadeApproximant(float[] taylorCoeffs, int m, int n)
        {
            if (taylorCoeffs.Length < m + n + 1)
                throw new ArgumentException("Not enough Taylor series coefficients for the specified orders.");

            float[,] A = new float[n, n];
            float[] b = new float[n];

            for (int i = 0; i < n; i++)
            {
                b[i] = -taylorCoeffs[m + i + 1];

                for (int j = 0; j < n; j++)
                {
                    A[i, j] = taylorCoeffs[m + i - j];
                }
            }

            float[] q = Matrice.Solve(A, b);

            denominatorCoeffs = new float[n + 1];
            denominatorCoeffs[0] = 1.0f;

            for (int i = 0; i < n; i++)
            {
                denominatorCoeffs[i + 1] = q[i];
            }

            numeratorCoeffs = new float[m + 1];

            for (int k = 0; k <= m; k++)
            {
                float sum = 0.0f;

                for (int j = 0; j <= Math.Min(k, n); j++)
                {
                    sum += denominatorCoeffs[j] * taylorCoeffs[k - j];
                }
                numeratorCoeffs[k] = sum;
            }
        }

        /// <summary>
        /// Evaluates a function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Value</returns>
        public float Evaluate(float x)
        {
            float num = 0;
            float den = 0;

            for (int i = numeratorCoeffs.Length - 1; i >= 0; i--)
            {
                num = num * x + numeratorCoeffs[i];
            }

            for (int i = denominatorCoeffs.Length - 1; i >= 0; i--)
            {
                den = den * x + denominatorCoeffs[i];
            }

            return num / den;
        }
    }

}
