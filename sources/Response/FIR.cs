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
        private double[] b;
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
        public FIR(double[] b)
        {
            B = b;
        }
        /// <summary>
        /// Gets or sets the array of signal coefficients.
        /// </summary>
        public double[] B
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
        public double[] Reaction(double[] u)
        {
            int length = u.Length;
            double[] y = new double[length];

            double input;
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
        public double[] Amplitude(double[] w)
        {
            int i, j, length = b.Length;
            Complex K1;
            double[] amplitude = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;

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
        public double[] Phase(double[] w)
        {
            int j, i, length = b.Length;
            Complex K1;
            double[] phase = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;

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
        /// <returns>Double precision floating point number</returns>
        public double Amplitude(double w)
        {
            int i;
            int length = b.Length;
            Complex K1 = new Complex(0, 0);

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
        /// <returns>Double precision floating point number</returns>
        public double Phase(double w)
        {
            int i;
            int length = b.Length;
            Complex K1 = new Complex(0, 0);

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
                cv.B = new double[3] { 1.0, 1.0, 0.0 };
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
                cv.B = new double[3] { 1.0, -1.0, 0.0 };
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
                cv.B = new double[3] { 1.0, 1.0, -1.0 };
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
                cv.B = new double[3] { 1.0, -1.0, 1.0 };
                return cv;
            }
        }
        #endregion
    }
}
