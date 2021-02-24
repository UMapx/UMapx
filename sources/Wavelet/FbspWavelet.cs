using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the continuous complex frequency B-spline wavelet.
    /// </summary>
    [Serializable]
    public class FbspWavelet : IComplexWavelet
    {
        #region Private data
        private float m;
        private float fb;
        private float fc;
        #endregion

        #region Wavelet compoents
        /// <summary>
        /// Initializes the continuous complex frequency B-spline wavelet.
        /// </summary>
        /// <param name="m">Order</param>
        /// <param name="fb">Bandwidth</param>
        /// <param name="fc">Center frequency</param>
        public FbspWavelet(float m = 3, float fb = 1, float fc = 2)
        {
            M = m; Fb = fb; Fc = fc;
        }
        /// <summary>
        /// Gets or sets the value of the wavelet order.
        /// </summary>
        public float M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the bandwidth parameter.
        /// </summary>
        public float Fb
        {
            get
            {
                return this.fb;
            }
            set
            {
                this.fb = value;
            }
        }
        /// <summary>
        /// Gets or sets the center frequency value of the wavelet.
        /// </summary>
        public float Fc
        {
            get
            {
                return this.fc;
            }
            set
            {
                this.fc = value;
            }
        }
        /// <summary>
        /// Returns the value of the scaling function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex32 Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex32 Wavelet(float x)
        {
            float a = (float)Math.Sqrt(fb);
            float b = x / (float)Math.Pow(fb, m);
            float c = Special.Sinc(b, 1);
            float d = (float)Math.Pow(c, m);
            Complex32 e = Maths.Exp(Maths.I * 2 * Maths.Pi * fc * x);
            return a * d * e;
        }
        #endregion
    }
}
