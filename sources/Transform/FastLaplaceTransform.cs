using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Laplace transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Laplace_transform
    /// </remarks>
    [Serializable]
    public class FastLaplaceTransform : TransformBaseComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        /// <summary>
        /// Damping factor.
        /// </summary>
        private float sigma;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Laplace transform.
        /// </summary>
        /// <param name="sigma">Damping factor (0, 1)</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastLaplaceTransform(float sigma = 0.0005f, bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, Direction.Vertical);
            this.Sigma = sigma; 
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
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        #endregion

        #region Laplace Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            // Fourier transform:
            Complex32[] B = FFT.Forward(A);

            // Fourier to Laplace transform:
            ApplyLaplaceOperator(B, sigma);
            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            // Laplace to Fourier transform:
            Complex32[] A = (Complex32[])B.Clone();
            ApplyLaplaceInverseOperator(A, sigma);

            // Fourier transform:
            return FFT.Backward(A);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Applies Laplace transform operator.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="sigma">Sigma</param>
        private static void ApplyLaplaceOperator(Complex32[] v, float sigma = 0.003f)
        {
            // forward laplacian transform
            // for signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = Math.Exp(-sigma * i) * v[i];
            }
        }
        /// <summary>
        /// Applies inverse Laplace transform operator.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="sigma">Sigma</param>
        private static void ApplyLaplaceInverseOperator(Complex32[] v, float sigma = 0.003f)
        {
            // inverse laplacian transform
            // for signal:
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                v[i] = v[i] / Math.Exp(-sigma * i);
            }
        }
        #endregion
    }
}
