using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Laplace transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    [Serializable]
    public class LaplaceTransform : TransformBaseMatrixComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Damping factor.
        /// </summary>
        private float sigma;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Laplace transform.
        /// </summary>
        /// <param name="sigma">Damping factor (0, 1)</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public LaplaceTransform(float sigma = 0.0005f, bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.Sigma = sigma; 
            this.Normalized = normalized; 
            this.Direction = direction;
        }
        /// <summary>
        /// Gets or sets the damping factor (0, 1).
        /// </summary>
        /// <remarks>
        /// If σ = 0, then the Laplace transform takes the form of a Fourier transform.
        /// </remarks>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Invalid argument value");

                this.sigma = value;
            }
        }
        #endregion

        #region Laplace static components
        /// <summary>
        /// Implements the construction of the Laplace matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="sigma">Damping factor (0, 1)</param>
        /// <param name="backward">Return backward transformation matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Matrix(int n, float sigma, bool backward = false)
        {
            Complex32[,] H = new Complex32[n, n];
            float factor;
            int i, j;

            // inverse matrix or not?
            if (backward)
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i) / factor;
                    }
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    factor = Maths.Exp(-sigma * i);

                    for (j = 0; j < n; j++)
                    {
                        H[i, j] = factor * Maths.Exp(-2 * Maths.Pi * Maths.I * j / n * i);
                    }
                }
            }
            return H;
        }
        #endregion

        #region Laplace Transform
        /// <inheritdoc/>
        protected override float NormalizationFactor(int n, bool backward = false)
        {
            return Maths.Sqrt(n);
        }
        /// <inheritdoc/>
        protected override Complex32[,] Matrix(int n, bool backward)
        {
            var U = Matrix(n, sigma, backward);

            if (backward)
            {
                U = U.Hermitian();
            }

            return U;
        }
        #endregion
    }
}
