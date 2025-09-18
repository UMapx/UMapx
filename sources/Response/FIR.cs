using System;
using UMapx.Core;

namespace UMapx.Response
{
    /// <summary>
    /// Defines a filter with a finite impulse response.
    /// </summary>
    /// <remarks>
    /// A filter with a finite impulse response, also called a transversal (FIR) filter, is one type of linear
    /// digital filter whose impulse response becomes exactly zero after a finite time.
    /// Such a filter is also called non-recursive due to the lack of feedback.
    /// The denominator of the transfer function of such a filter is a constant.
    /// </remarks>
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
            get { return this.b; }
            set { this.b = value; }
        }
        /// <summary>
        /// Returns an array of filter response values when a discrete function is supplied.
        /// </summary>
        /// <param name="u">Array</param>
        /// <returns>Array</returns>
        public float[] Reaction(float[] u)
        {
            int length = u.Length;
            float[] y = new float[length];

            int P = b.Length;

            for (int n = 0; n < length; n++)
            {
                float accX = 0f;

                for (int i = 0; i < P; i++)
                {
                    int t = n - i;
                    if (t < 0) break;
                    accX += b[i] * u[t];
                }

                y[n] = accX;
            }
            return y;
        }
        /// <summary>
        /// Returns an array of filter response values when a discrete function is supplied.
        /// </summary>
        /// <param name="u">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Reaction(Complex32[] u)
        {
            int length = u.Length;
            Complex32[] y = new Complex32[length];

            int P = b.Length;

            for (int n = 0; n < length; n++)
            {
                Complex32 accX = 0f;

                for (int i = 0; i < P; i++)
                {
                    int t = n - i;
                    if (t < 0) break;
                    accX += b[i] * u[t];
                }

                y[n] = accX;
            }
            return y;
        }
        /// <summary>
        /// Returns the frequency response of the filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / sample)</param>
        /// <returns>Array</returns>
        public float[] Amplitude(float[] w)
        {
            int length = b.Length;
            float[] amplitude = new float[w.Length];

            for (int j = 0; j < w.Length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                for (int i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                amplitude[j] = K1.Abs;
            }
            return amplitude;
        }
        /// <summary>
        /// Returns the frequency response of the filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / sample)</param>
        /// <returns>Array</returns>
        public Complex32[] Amplitude(Complex32[] w)
        {
            int length = b.Length;
            Complex32[] amplitude = new Complex32[w.Length];

            for (int j = 0; j < w.Length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                for (int i = 0; i < length; i++)
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
        /// <param name="w">Array of frequencies (rad / sample)</param>
        /// <returns>Array</returns>
        public float[] Phase(float[] w)
        {
            int length = b.Length;
            float[] phase = new float[w.Length];

            for (int j = 0; j < w.Length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                for (int i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                phase[j] = K1.Angle;
            }
            return phase;
        }
        /// <summary>
        /// Returns the phase-frequency response of a filter.
        /// </summary>
        /// <param name="w">Array of frequencies (rad / sample)</param>
        /// <returns>Array</returns>
        public Complex32[] Phase(Complex32[] w)
        {
            int length = b.Length;
            Complex32[] phase = new Complex32[w.Length];

            for (int j = 0; j < w.Length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                for (int i = 0; i < length; i++)
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
        /// <param name="w">Frequency (rad / sample)</param>
        /// <returns>Value</returns>
        public float Amplitude(float w)
        {
            Complex32 K1 = Complex32.Zero;
            int length = b.Length;

            for (int i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Abs;
        }
        /// <summary>
        /// Returns the amplitude value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / sample)</param>
        /// <returns>Value</returns>
        public Complex32 Amplitude(Complex32 w)
        {
            Complex32 K1 = Complex32.Zero;
            int length = b.Length;

            for (int i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Abs;
        }
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / sample)</param>
        /// <returns>Value</returns>
        public float Phase(float w)
        {
            Complex32 K1 = Complex32.Zero;
            int length = b.Length;

            for (int i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Angle;
        }
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / sample)</param>
        /// <returns>Value</returns>
        public Complex32 Phase(Complex32 w)
        {
            Complex32 K1 = Complex32.Zero;
            int length = b.Length;

            for (int i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Angle;
        }
        /// <summary>
        /// Checks if the specified filter is stable.
        /// </summary>
        public bool Stability
        {
            get
            {
                var l1Gain = 0.0f;

                foreach (var t in b)
                {
                    if (float.IsNaN(t) || float.IsInfinity(t)) return false;
                    l1Gain += Maths.Abs(t);
                    if (float.IsInfinity(l1Gain)) return false;
                }
                return true;
            }
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
