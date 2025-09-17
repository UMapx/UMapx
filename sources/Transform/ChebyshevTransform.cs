using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Discrete Chebyshev transform of the first kind (Type-I), orthonormal.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Chebyshev_transform
    /// </remarks>
    [Serializable]
    public class ChebyshevTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the Chebyshev transform (Type-I, orthonormal).
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public ChebyshevTransform(Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
        }
        #endregion

        #region Chebyshev static components
        /// <summary>
        /// Constructs the orthonormal Chebyshev (DCT-I) matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            if (n == 1)
                return new float[1, 1] { { 1f } };

            int M = n - 1;
            float[,] H = new float[n, n];
            float cm = Maths.Pi / M;                // π / M
            float norm = Maths.Sqrt(2.0f / M);      // √(2/M)
            float invSqrt2 = 1.0f / Maths.Sqrt2;    // 1/√2

            for (int k = 0; k < n; k++)
            {
                float ak = (k == 0 || k == M) ? invSqrt2 : 1.0f;
                for (int j = 0; j < n; j++)
                {
                    float aj = (j == 0 || j == M) ? invSqrt2 : 1.0f;
                    H[k, j] = norm * ak * aj * Maths.Cos(cm * k * j);
                }
            }
            return H;
        }
        #endregion

        #region Chebyshev Transform
        /// <inheritdoc/>
        protected override float[,] TransformationMatrix(int n)
        {
            return Matrix(n);
        }
        #endregion
    }
}
