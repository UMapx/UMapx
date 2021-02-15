using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the gamma correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gamma_correction
    /// </remarks>
    /// </summary>
    [Serializable]
    public class GammaCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private float g;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the gamma correction filter.
        /// </summary>
        /// <param name="g">Gamma [0, 20]</param>
        /// <param name="space">Color space</param>
        public GammaCorrection(float g, Space space)
        {
            Gamma = g; Space = space;
        }
        /// <summary>
        /// Initializes the gamma correction filter.
        /// </summary>
        public GammaCorrection()
        {
            Gamma = 2.2f;
        }
        /// <summary>
        /// Gets or sets the gamma value [0, 20].
        /// </summary>
        public float Gamma
        {
            get
            {
                return g;
            }
            set
            {
                g = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Gamma(this.g, 256);
        }
        #endregion
    }
}
