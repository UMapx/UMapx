using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hilbert transform.
    /// </summary>
    /// <remarks>
    /// NOT RECOMMENDED.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    [Serializable]
    public class HilbertTransform : TransformBaseComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FourierTransform DFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hilbert transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public HilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.DFT = new FourierTransform(normalized, Direction.Both);
            this.Direction = direction;
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

        #region Hilbert Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            var F = DFT.Forward(A);
            ApplyHilbertOperatorMaskInplace(F);
            return DFT.Backward(F);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            var Hx = Forward(B);
            for (int i = 0; i < Hx.Length; i++) Hx[i] = -Hx[i];
            return Hx;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Applies Hilber spectrum operator.
        /// </summary>
        /// <param name="F">Array</param>
        internal static void ApplyHilbertOperatorMaskInplace(Complex32[] F)
        {
            int N = F.Length;
            if (N == 0) return;

            // DC
            F[0] = Complex32.Zero;

            if ((N & 1) == 0)
            {
                int half = N / 2;

                // Nyquist
                F[half] = Complex32.Zero;

                // positive freqs: 1..half-1  -> -i
                for (int k = 1; k < half; k++) F[k] *= new Complex32(0f, -1f);

                // negative freqs: half+1..N-1 -> +i
                for (int k = half + 1; k < N; k++) F[k] *= new Complex32(0f, +1f);
            }
            else
            {
                int half = N / 2; // floor
                                  // positive: 1..half -> -i
                for (int k = 1; k <= half; k++) F[k] *= new Complex32(0f, -1f);
                // negative: half+1..N-1 -> +i
                for (int k = half + 1; k < N; k++) F[k] *= new Complex32(0f, +1f);
            }
        }
        #endregion
    }
}
