using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Fourier transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_Fourier_transform
    /// </remarks>
    [Serializable]
    public class FourierTransform : TransformBaseMatrixComplex32, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the Fourier transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FourierTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.Normalized = normalized;
            this.Direction = direction;
        }
        #endregion

        #region Fourier static components
        /// <summary>
        /// Implements the construction of the Fourier matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Matrix(int n)
        {
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * i * j / n);
                }
            }
            return H;
        }
        #endregion

        #region Fourier Transform
        /// <inheritdoc/>
        protected override float ForwardNormalizationFactor(int n)
        {
            return Maths.Sqrt(n);
        }
        /// <inheritdoc/>
        protected override float BackwardNormalizationFactor(int n)
        {
            return Maths.Sqrt(n);
        }
        /// <inheritdoc/>
        protected override Complex32[,] TransformationMatrix(int n)
        {
            return Matrix(n);
        }
        #endregion
    }
}
