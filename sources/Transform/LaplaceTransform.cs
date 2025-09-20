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
    public class LaplaceTransform : TransformBaseComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FourierTransform DFT;
        /// <summary>
        /// Damping factor.
        /// </summary>
        private float sigma;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Laplace transform.
        /// </summary>
        /// <param name="sigma">Non-negative damping factor</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public LaplaceTransform(float sigma = 0.0005f, bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.DFT = new FourierTransform(normalized, Direction.Vertical);
            this.Sigma = sigma;
            this.Direction = direction;
        }
        /// <summary>
        /// Gets or sets the non-negative damping factor.
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
        public override bool Normalized
        {
            get
            {
                return this.DFT.Normalized;
            }
            set
            {
                this.DFT.Normalized = value;
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
            // Apply Laplace damping before Fourier transform:
            Complex32[] weighted = (Complex32[])A.Clone();
            LaplaceTransform.ApplyLaplaceOperator(weighted, sigma);

            // Fourier transform of damped signal:
            return DFT.Forward(weighted);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            // Fourier transform:
            Complex32[] A = DFT.Backward(B);

            // Compensate Laplace damping:
            LaplaceTransform.ApplyLaplaceInverseOperator(A, sigma);
            return A;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Applies Laplace transform operator.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="sigma">Sigma</param>
        internal static void ApplyLaplaceOperator(Complex32[] v, float sigma = 0.003f)
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
        internal static void ApplyLaplaceInverseOperator(Complex32[] v, float sigma = 0.003f)
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
