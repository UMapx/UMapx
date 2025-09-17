using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the cosine transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_cosine_transform
    /// </remarks>
    [Serializable]
    public class CosineTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the cosine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public CosineTransform(Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
        }
        #endregion

        #region Cosine static components
        /// <summary>
        /// Implements the construction of the cosine transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            int j, i;
            float[,] H = new float[n, n];
            float c = Maths.Pi / (2.0f * n);
            float g1 = Maths.Sqrt(1.0f / n);
            float g2 = Maths.Sqrt(2.0f / n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = (i == 0) ? g1 : Maths.Cos((2 * j + 1) * i * c) * g2;
                }
            }

            return H;
        }
        #endregion

        #region Cosine Transform
        /// <inheritdoc/>
        public override float[,] Matrix(int n, bool backward)
        {
            var U = Matrix(n);

            if (backward)
            {
                U = U.Transpose();
            }

            return U;
        }
        #endregion
    }
}
