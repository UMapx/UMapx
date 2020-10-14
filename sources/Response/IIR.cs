﻿using System;
using UMapx.Core;

namespace UMapx.Response
{
    /// <summary>
    /// Defines a filter with an infinite impulse response.
    /// <remarks>
    /// Filter with infinite impulse response (recursive filter, IIR filter or IIR filter) - a linear electronic filter,
    /// using one or more of its outputs as an input, i.e. forms a feedback. The main property of such filters
    /// is that their impulse response is of infinite length in the time domain, and the transfer function
    /// has a fractional rational look.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class IIR : IResponse
    {
        #region Private data
        private double[] a;
        private double[] b;
        #endregion

        #region IIR Components
        /// <summary>
        /// Initializes a filter with an infinite impulse response.
        /// </summary>
        public IIR() { }
        /// <summary>
        /// Initializes a filter with an infinite impulse response.
        /// </summary>
        /// <param name="b">Array of signal coefficients</param>
        /// <param name="a">Array of feedback coefficients</param>
        public IIR(double[] b, double[] a)
        {
            B = b; A = a;
        }
        /// <summary>
        /// Gets or sets the array of feedback coefficients.
        /// </summary>
        public double[] A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
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
            double input, output;
            int t, P = b.Length, Q = a.Length;
            int n, i, k;

            for (n = 0; n < length; n++)
            {
                input = 0; output = 0;

                for (i = 0; i < P; i++)
                {
                    t = n - i;
                    if (t < 0) continue;
                    input += b[i] * u[t];
                }

                for (k = 0; k < Q; k++)
                {
                    t = n - k;
                    if (t < 0) continue;
                    output += a[k] * y[t];
                }

                y[n] = input + output;

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
            int i, j;
            Complex K1;
            Complex K2;
            int length = w.Length;
            int P = b.Length, Q = a.Length;
            double[] amplitude = new double[length];

            for (j = 0; j < length; j++)
            {
                K1 = Complex.Zero;
                K2 = Complex.One;

                for (i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                amplitude[j] = K1.Abs / K2.Abs;
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
            int i, j;
            Complex K1;
            Complex K2;
            int length = w.Length;
            int P = b.Length, Q = a.Length;
            double[] phase = new double[length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;
                K2 = Complex.One;

                for (i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                phase[j] = K1.Angle - K2.Angle;
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
            Complex K1 = new Complex(0, 0);
            Complex K2 = new Complex(1, 0);

            for (i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Abs / K2.Abs;
        }
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / s)</param>
        /// <returns>Double precision floating point number</returns>
        public double Phase(double w)
        {
            int i;
            Complex K1 = new Complex(0, 0);
            Complex K2 = new Complex(1, 0);
            int P = b.Length, Q = a.Length;

            for (i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Angle - K2.Angle;
        }
        /// <summary>
        /// Checks if the specified filter is stable.
        /// </summary>
        public bool Stability
        {
            get
            {
                Complex sum = Complex.Zero;
                Complex exp = Maths.Exp(-Maths.I);
                int length = a.Length;

                for (int j = 0; j < length; j++)
                {
                    sum += a[j] * exp;
                }

                if (sum == Complex.Zero)
                {
                    return true;
                }
                return false;
            }
        }
        #endregion

        #region Sample filters
        /// <summary>
        /// Gets the finished low pass filter.
        /// </summary>
        public static IIR LowPass
        {
            get
            {
                IIR cv = new IIR();
                cv.B = new double[3] { 1.0, 1.0, 0.0 };
                cv.A = new double[3] { 0.0, 0.5, 0.5 };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished high-pass filter.
        /// </summary>
        public static IIR HighPass
        {
            get
            {
                IIR cv = new IIR();
                cv.B = new double[3] { 1.0, -1.0, 0.0 };
                cv.A = new double[3] { 0.0, 0.5, 0.5 };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished bandpass filter.
        /// </summary>
        public static IIR BandPass
        {
            get
            {
                IIR cv = new IIR();
                cv.B = new double[3] { 1.0, 1.0, -1.0 };
                cv.A = new double[3] { 0.0, 0.5, -0.5 };
                return cv;
            }
        }
        /// <summary>
        /// Gets the finished notch filter.
        /// </summary>
        public static IIR Notch
        {
            get
            {
                IIR cv = new IIR();
                cv.B = new double[3] { 1.0, -1.0, 1.0 };
                cv.A = new double[3] { 0.0, 0.5, 0.5 };
                return cv;
            }
        }
        #endregion
    }
}