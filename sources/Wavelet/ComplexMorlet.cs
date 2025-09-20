using System;
using UMapx.Core;

namespace UMapx.Wavelet
{
    /// <summary>
    /// Defines the normalized continuous complex Morlet wavelet.
    /// </summary>
    [Serializable]
    public class ComplexMorlet : IWaveletComplex32
    {
        #region Private data
        private float fb;
        private float fc;
        #endregion

        #region Wavelet components
        /// <summary>
        /// Initializes the continuous complex Morlet wavelet.
        /// </summary>
        /// <param name="fb">Bandwidth</param>
        /// <param name="fc">Center frequency</param>
        public ComplexMorlet(float fb = 0.5f, float fc = 1)
        {
            Fb = fb; Fc = fc;
        }
        /// <summary>
        /// Gets or sets the bandwidth.
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
        /// Gets or sets the center frequency.
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
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public Complex32 Scaling(float x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Returns the value of the normalized wavelet function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>Function</returns>
        public Complex32 Wavelet(float x)
        {
            // Validate (optional but recommended)
            if (fb <= 0)
                throw new ArgumentException("Fb (variance) must be > 0");
            if (fc < 0)
                throw new ArgumentException("Fc (center frequency) must be >= 0");

            // Normalization: (π·fb)^(-1/4)
            float amplitude = Maths.Pow(Maths.Pi * fb, -0.25f);

            // Zero-mean correction for admissibility with Gaussian exp(-x^2 / (2·fb))
            float correction = Maths.Exp(-2f * Maths.Pi * Maths.Pi * fb * fc * fc);

            // Complex carrier exp(i·2π·fc·x)
            Complex32 carrier = Maths.Exp(2f * Maths.Pi * Maths.I * fc * x);

            // Gaussian envelope with variance fb (σ² = fb)
            float gaussian = Maths.Exp(-(x * x) / (2f * fb));

            return amplitude * (carrier - correction) * gaussian;
        }
        #endregion
    }
}
