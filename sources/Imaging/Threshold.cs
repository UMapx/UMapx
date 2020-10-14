using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the threshold filter.
    /// </summary>
    [Serializable]
    public class Threshold : Correction, IBitmapFilter
    {
        #region Private data
        private double threshold;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <param name="space">Color space</param>
        public Threshold(double threshold, Space space)
        {
            this.Value = threshold;
            this.Space = space;
        }
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        public Threshold()
        {
            Value = 0.5;
        }
        /// <summary>
        /// Gets or sets the threshold value [0, 1].
        /// </summary>
        public double Value
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bin(this.threshold, 256);
        }
        #endregion
    }
}
