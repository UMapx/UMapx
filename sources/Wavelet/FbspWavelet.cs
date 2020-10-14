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
        private double m;
        private double fb;
        private double fc;
        #endregion

        #region Wavelet compoents
        /// <summary>
        /// Initializes the continuous complex frequency B-spline wavelet.
        /// </summary>
        /// <param name="m">Order</param>
        /// <param name="fb">Bandwidth</param>
        /// <param name="fc">Center frequency</param>
        public FbspWavelet(double m = 3, double fb = 1, double fc = 2)
        {
            M = m; Fb = fb; Fc = fc;
        }
        /// <summary>
        /// Gets or sets the value of the wavelet order.
        /// </summary>
        public double M
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
        public double Fb
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
        public double Fc
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
        public Complex Scaling(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the wavelet function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Function</returns>
        public Complex Wavelet(double x)
        {
            double a = Math.Sqrt(fb);
            double b = x / Math.Pow(fb, m);
            double c = Special.Sinc(b, 1);
            double d = Math.Pow(c, m);
            Complex e = Maths.Exp(Maths.I * 2 * Maths.Pi * fc * x);
            return a * d * e;
        }
        #endregion
    }
}
