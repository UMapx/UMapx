using System;
using UMapx.Core;

namespace UMapx.Response
{
    /// <summary>
    /// Defines a filter with a finite impulse response.
    /// <remarks>
    /// A filter with a finite impulse response (transverse filter, FIR filter or FIR filter) is one of the types of linear
    /// digital filters, a characteristic feature of which is the limited time of its impulse response
    /// (from some point in time it becomes exactly equal to zero). Such a filter is also called non-recursive due to the lack of feedback.
    /// The denominator of the transfer function of such a filter is a certain constant.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FIR : IResponse
    {
        #region Private data
        private float[] b;
        #endregion

        #region FIR Components
        /// <summary>
        /// Initializes a filter with a finite impulse response.
        /// </summary>
        public FIR() { }
        /// <summary>
        /// Initializes a filter with a finite impulse response.
        /// </summary>
        /// <param name="b">Array of signal coefficients</param>
        public FIR(float[] b)
        {
            B = b;
        }
        /// <summary>
        /// Gets or sets the array of signal coefficients.
        /// </summary>
        public float[] B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Returns an array of filter response values when a discrete function is supplied.
        /// </summary>
        /// <param name="u">Array</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        public float[] Reaction(float[] u)
        {
            int length = u.Length;
            float[] y = new float[length];

            float input;
            int t, P = b.Length;
            int n, i;

            for (n = 0; n < length; n++)
            {
                input = 0;

                for (i = 0; i < P; i++)
                {
                    t = n - i;
                    if (t < 0) continue;
                    input += b[i] * u[t];
                }

                y[n] = input;

            }
            return y;
        }
        /// <summary>
        /// Returns the frequency response of the filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / s)</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        public float[] Amplitude(float[] w)
        {
            int i, j, length = b.Length;
            Complex32 K1;
            float[] amplitude = new float[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex32.Zero;

                for (i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                amplitude[j] = K1.Abs;
            }
            return amplitude;
        }
        /// <summary>
        /// Returns the phase-frequency response of a filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / s)</param>
        /// <returns>Discrete function in a Cartesian coordinate system</returns>
        public float[] Phase(float[] w)
        {
            int j, i, length = b.Length;
            Complex32 K1;
            float[] phase = new float[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex32.Zero;

                for (i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                phase[j] = K1.Angle;
            }
            return phase;
        }
        /// <summary>
        /// Returns the amplitude value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / s)</param>
        /// <returns>Value</returns>
        public float Amplitude(float w)
        {
            int i;
            int length = b.Length;
            Complex32 K1 = new Complex32(0, 0);

            for (i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Abs;
        }
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / s)</param>
        /// <returns>Value</returns>
        public float Phase(float w)
        {
            int i;
            int length = b.Length;
            Complex32 K1 = new Complex32(0, 0);

            for (i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Angle;
        }
        #endregion

        #region Sample filters
        /// <summary>
        /// Gets the finished low pass filter.
        /// </summary>
        public static FIR LowPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new float[3] { 1.0f, 1.0f, 0.0f };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished high-pass filter.
        /// </summary>
        public static FIR HighPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new float[3] { 1.0f, -1.0f, 0.0f };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished bandpass filter.
        /// </summary>
        public static FIR BandPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new float[3] { 1.0f, 1.0f, -1.0f };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished notch filter.
        /// </summary>
        public static FIR Notch
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new float[3] { 1.0f, -1.0f, 1.0f };
                return cv;
            }
        }
        #endregion
    }
}
