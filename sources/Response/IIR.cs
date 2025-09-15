using System;
using UMapx.Analysis;
using UMapx.Core;

namespace UMapx.Response
{
    /// <summary>
    /// Defines a filter with an infinite impulse response.
    /// </summary>
    /// <remarks>
    /// Filter with infinite impulse response (recursive filter, IIR filter) — a linear digital filter with feedback.
    /// The impulse response is of infinite length in time, and the transfer function is rational.
    /// Convention used here:
    ///     H(e^{i w}) = (Σ_{i=0}^{P-1} b[i] e^{-i w i}) / (1 - Σ_{k=0}^{Q-1} a[k] e^{-i w k})
    /// which leads to the difference equation:
    ///     y[n] = ( Σ_i b[i] x[n-i] + Σ_{k=1}^{Q-1} a[k] y[n-k] ) / (1 - a[0])
    /// </remarks>
    [Serializable]
    public class IIR : IResponse
    {
        #region Private data
        private float[] a;
        private float[] b;
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
        public IIR(float[] b, float[] a)
        {
            B = b; A = a;
        }
        /// <summary>
        /// Gets or sets the array of feedback coefficients.
        /// </summary>
        public float[] A
        {
            get { return this.a; }
            set { this.a = value; }
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

            int P = (b != null) ? b.Length : 0;
            int Q = (a != null) ? a.Length : 0;

            float a0 = (Q > 0) ? a[0] : 0f;
            float denom = 1f - a0;

            if (denom == 0f)
                throw new InvalidOperationException("IIR.Reaction: a[0] = 1 → denominator is zero (pole on z=1)");

            for (int n = 0; n < length; n++)
            {
                float accX = 0f, accY = 0f;

                // feedforward (x-part)
                for (int i = 0; i < P; i++)
                {
                    int t = n - i;
                    if (t < 0) break;
                    accX += b[i] * u[t];
                }

                // feedback (y-part), starts from k = 1
                for (int k = 1; k < Q; k++)
                {
                    int t = n - k;
                    if (t < 0) break;
                    accY += a[k] * y[t];
                }

                y[n] = (accX + accY) / denom;
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
            int length = w.Length;
            int P = b.Length, Q = a.Length;
            float[] amplitude = new float[length];

            for (int j = 0; j < length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                Complex32 K2 = Complex32.One;

                for (int i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (int i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                amplitude[j] = K1.Abs / K2.Abs;
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
            int length = w.Length;
            int P = b.Length, Q = a.Length;
            float[] phase = new float[length];

            for (int j = 0; j < length; j++)
            {
                Complex32 K1 = Complex32.Zero;
                Complex32 K2 = Complex32.One;

                for (int i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (int i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                phase[j] = K1.Angle - K2.Angle;
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
            Complex32 K2 = Complex32.One;

            for (int i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (int i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Abs / K2.Abs;
        }
        /// <summary>
        /// Returns the phase value at the given frequency.
        /// </summary>
        /// <param name="w">Frequency (rad / sample)</param>
        /// <returns>Value</returns>
        public float Phase(float w)
        {
            Complex32 K1 = Complex32.Zero;
            Complex32 K2 = Complex32.One;

            int P = b.Length, Q = a.Length;

            for (int i = 0; i < P; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (int i = 0; i < Q; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Angle - K2.Angle;
        }
        /// <summary>
        /// Checks if the specified filter is stable (all poles inside the unit circle).
        /// </summary>
        public bool Stability
        {
            get
            {
                if (a == null || a.Length == 0) return true;

                int Q = a.Length;
                float[] p = new float[Q + 1];

                p[0] = 1f;
                for (int k = 0; k < Q; k++)
                    p[k + 1] = -a[k];

                var eps = 1e-8f;
                var rts = new Roots(eps).Compute(p);  // p(1)*x^n + ... + p(n+1) = 0

                for (int i = 0; i < rts.Length; i++)
                {
                    if (rts[i].Abs >= 1f - eps)
                        return false;
                }
                return true;
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
                // Demo 2nd-order low-pass style (stable): D(z) = 1 - 0.5 z^{-1} - 0.25 z^{-2}
                IIR cv = new IIR();
                cv.B = new float[3] { 1.0f, 1.0f, 0.0f };
                cv.A = new float[3] { 0.0f, 0.5f, 0.25f };
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
                // Demo 2nd-order high-pass style (stable): same poles as LowPass, different zeros
                IIR cv = new IIR();
                cv.B = new float[3] { 1.0f, -1.0f, 0.0f };
                cv.A = new float[3] { 0.0f, 0.5f, 0.25f };
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
                // Demo (stable): complex-conjugate poles with radius 0.5
                IIR cv = new IIR();
                cv.B = new float[3] { 1.0f, 1.0f, -1.0f };
                cv.A = new float[3] { 0.0f, 0.5f, -0.5f };
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
                // Demo notch, stable poles (reuse LowPass poles)
                IIR cv = new IIR();
                cv.B = new float[3] { 1.0f, -1.0f, 1.0f };
                cv.A = new float[3] { 0.0f, 0.5f, 0.25f };
                return cv;
            }
        }
        #endregion
    }
}
