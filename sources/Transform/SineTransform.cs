using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the sine transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    [Serializable]
    public class SineTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the sine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public SineTransform(Direction direction = Direction.Vertical)
        {
            this.Direction = direction;
        }
        #endregion

        #region Sine static components
        /// <summary>
        /// Implements the construction of the sine transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int n)
        {
            int j, i;
            float[,] H = new float[n, n];
            float n1 = n + 1;
            float scale = Maths.Sqrt(2.0f / n1);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Sin(Maths.Pi * (j + 1) * (i + 1) / n1) * scale;
                }
            }

            return H;
        }
        #endregion

        #region Sine Transform
        /// <inheritdoc/>
        protected override float[,] Matrix(int n, bool backward)
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
