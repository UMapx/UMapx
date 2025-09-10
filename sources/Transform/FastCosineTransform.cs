using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the fast cosine transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Discrete_cosine_transform
    /// </remarks>
    [Serializable]
    public class FastCosineTransform : FloatToComplex32TransformBase, ITransform
    {
        #region Private data
        /// <summary>
        /// Fourier transform.
        /// </summary>
        private readonly FastFourierTransform FFT;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the fast cosine transform.
        /// </summary>
        /// <param name="direction">Processing direction</param>
        public FastCosineTransform(Direction direction = Direction.Vertical)
        {
            this.FFT = new FastFourierTransform(false, Direction.Both);
            this.Direction = direction;
        }
        #endregion

        #region Fast Cosine Transform
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N = A.Length, N2 = N / 2, i, k;
            Complex32[] B = new Complex32[N];

            for (i = 0; i < N2; i++)
            {
                var j = 2 * i;
                B[i     ] = A[j];
                B[i + N2] = A[N - j - 1];
            }

            B = FFT.Forward(B);

            float[] C = new float[N];

            Complex32 c = -Maths.I * Maths.Pi / (2f * N);

            for (k = 0; k < N; k++)
            {
                C[k] = 2.0f * (B[k] * Maths.Exp(c * k)).Real / Maths.Sqrt(2f * N);
            }

            C[0] = C[0] / Maths.Sqrt2;

            return C;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            int N = B.Length, N2 = N / 2, i, k;
            Complex32[] A = new Complex32[N];
            Complex32 c = Maths.I * Maths.Pi / (2f * N);

            for (k = 0; k < N; k++)
            {
                A[k] = B[k] * Maths.Exp(c * k) * Maths.Sqrt(2f * N);
            }

            A[0] /= Maths.Sqrt(2f);

            A = FFT.Backward(A);

            float[] C = new float[N];

            for (i = 0; i < N2; i++)
            {
                var j = 2 * i;
                C[j    ] = A[i        ].Real / N;
                C[j + 1] = A[N - i - 1].Real / N;
            }

            return C;
        }
        #endregion
    }
}
