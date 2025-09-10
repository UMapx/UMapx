using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines short-time Fourier transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    [Serializable]
    public class ShortTimeFourierTransform : TransformBaseComplex32, IWindowTransform, ITransform
    {
        #region Private data
        private readonly FourierTransform DFT;
        private IWindow window;
        private readonly float[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes short-time Fourier transform.
        /// </summary>
        /// <param name="function">Window function</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public ShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.DFT = new FourierTransform(normalized, direction);
            Direction = direction;
            Window = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-3f);
        }
        /// <summary>
        /// Normalized transform or not.
        /// </summary>
        public bool Normalized
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
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-time Fourier transform
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length, frame = coefs.Length;
            if (N % frame != 0)
                throw new ArgumentException("Length N must be a multiple of the frame (window) length");

            Complex32[] B = new Complex32[N];

            for (int i = 0; i < N; i += frame)
            {
                Complex32[] data = new Complex32[frame];

                for (int j = 0; j < frame; j++)
                    data[j] = A[i + j] * coefs[j];

                data = this.DFT.Forward(data);

                for (int j = 0; j < frame; j++)
                    B[i + j] = data[j];
            }

            return B;
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length, frame = coefs.Length;
            if (N % frame != 0)
                throw new ArgumentException("Length N must be a multiple of the frame (window) length");

            Complex32[] A = new Complex32[N];

            for (int i = 0; i < N; i += frame)
            {
                Complex32[] data = new Complex32[frame];

                for (int j = 0; j < frame; j++)
                    data[j] = B[i + j];

                data = this.DFT.Backward(data);

                for (int j = 0; j < frame; j++)
                    A[i + j] = data[j] / coefs[j];
            }

            return A;
        }
        #endregion
    }
}
