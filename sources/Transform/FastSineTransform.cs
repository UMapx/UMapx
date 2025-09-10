using System;
using System.Numerics;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast sine transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://sernam.ru/book_prett1.php?id=91
    /// </remarks>
    [Serializable]
    public class FastSineTransform : TransformBaseFloat, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Defines the fast sine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public FastSineTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(false, Direction.Both);
            this.Direction = direction;
        }
        #endregion

        #region Fast Sine Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int n = A.Length;
            int m = 2 * (n + 1);

            Complex32[] extended = new Complex32[m];
            extended[0] = Complex.Zero;

            for (int i = 0; i < n; i++)
            {
                extended[i + 1] = new Complex(-A[i], 0);
                extended[m - 1 - i] = new Complex(A[i], 0);
            }

            extended[n + 1] = Complex.Zero;
            extended = FFT.Forward(extended);

            float scale = Maths.Sqrt(1.0f / m);
            float[] output = new float[n];

            for (int i = 0; i < n; i++)
            {
                output[i] = scale * extended[i + 1].Imag;
            }

            return output;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            int n = B.Length;
            int m = 2 * (n + 1);

            Complex32[] spectrum = new Complex32[m];
            spectrum[0] = Complex.Zero;

            for (int i = 0; i < n; i++)
            {
                spectrum[i + 1] = new Complex(0, -B[i]);
                spectrum[m - 1 - i] = new Complex(0, B[i]);
            }

            spectrum[n + 1] = Complex.Zero;
            spectrum = FFT.Backward(spectrum);

            float scale = Maths.Sqrt(1.0f / m);
            float[] output = new float[n];

            for (int i = 0; i < n; i++)
            {
                output[i] = scale * spectrum[i + 1].Real;
            }

            return output;
        }
        #endregion
    }
}
