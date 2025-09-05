using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the exposure correction filter.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Exposure_(photography)
    /// </remarks>
    [Serializable]
    public class ExposureCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private float average;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the exposure correction filter.
        /// </summary>
        /// <param name="average">Average [0, 2500]</param>
        /// <param name="space">Color space</param>
        public ExposureCorrection(float average, Space space)
        {
            Average = average; this.Space = space;
        }
        /// <summary>
        /// Initializes the exposure correction filter.
        /// </summary>
        public ExposureCorrection()
        {
            Average = 128;
        }
        /// <summary>
        /// Gets or sets the average [0, 2500].
        /// </summary>
        public float Average
        {
            get
            {
                return this.average;
            }
            set
            {
                this.average = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Exposure(this.average, 256);
        }
        #endregion
    }
}
