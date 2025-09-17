using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast Hilbert transform.
    /// </summary>
    /// <remarks>
    /// NOT RECOMMENDED.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hilbert_transform
    /// </remarks>
    [Serializable]
    public class FastHilbertTransform : TransformBaseComplex32, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast Hilbert transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public FastHilbertTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(normalized, Direction.Both);
            this.Direction = direction;
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public override bool Normalized
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

        #region Fast Hilbert Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            var F = FFT.Forward(A);
            HilbertTransform.ApplyHilbertOperatorMaskInplace(F);
            return FFT.Backward(F);
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
    }
}
