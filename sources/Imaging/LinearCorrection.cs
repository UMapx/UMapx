using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the linear correction filter.
    /// </summary>
    [Serializable]
    public class LinearCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble range;
        private double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        /// <param name="range">Range values</param>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public LinearCorrection(RangeDouble range, double delta, Space space)
        {
            Range = range; Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        /// <param name="delta">Delta [-100, 100]</param>
        /// <param name="space">Color space</param>
        public LinearCorrection(double delta, Space space)
        {
            Range = new RangeDouble(0, 1); Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        public LinearCorrection()
        {
            Range = new RangeDouble(0, 1); Delta = 0.5; this.Space = space;
        }
        /// <summary>
        /// Gets or sets range values.
        /// </summary>
        public RangeDouble Range
        {
            get
            {
                return this.range;
            }
            set
            {
                this.range = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public double Delta
        {
            get
            {
                return this.delta;
            }
            set
            {
                this.delta = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Linear(range, delta / 2.0, 256);
        }
        #endregion
    }
}
