using System;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the temperature correction filter.
    /// <remarks>
    /// The filter uses an approximation of the Planck curve.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class TemperatureCorrection : PhotoFilter, IBitmapFilter
    {
        #region Private data
        private double temperature;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the temperature correction filter.
        /// </summary>
        /// <param name="temperature">Temperature [1E3K, 1E4K]</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public TemperatureCorrection(double temperature, double strength = 0.5)
        {
            Temperature = temperature; Strength = strength;
        }
        /// <summary>
        /// Initializes the temperature correction filter.
        /// </summary>
        public TemperatureCorrection()
        {
            Temperature = 1000; Strength = 0.5;
        }
        /// <summary>
        /// Gets or sets the temperature [1E3K, 1E4K].
        /// </summary>
        public double Temperature
        {
            get
            {
                return this.temperature;
            }
            set
            {
                this.temperature = value;
                this.Color = RGB.Temp2RGB(this.temperature);
            }
        }
        #endregion
    }
}
