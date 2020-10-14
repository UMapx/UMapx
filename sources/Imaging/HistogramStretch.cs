using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the global histogram stretch filter.
    /// </summary>
    [Serializable]
    public class HistogramStretch : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble range;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the global histogram stretch filter.
        /// </summary>
        /// <param name="min">Minimum intensity [0, 1]</param>
        /// <param name="max">Maximum intensity [0, 1]</param>
        /// <param name="space">Color space</param>
        public HistogramStretch(double min, double max, Space space)
        {
            Range = new RangeDouble(min, max);
            Space = space;
        }
        /// <summary>
        /// Initializes the global histogram stretch filter.
        /// </summary>
        /// <param name="range">Intensity range</param>
        /// <param name="space">Color space</param>
        public HistogramStretch(RangeDouble range, Space space)
        {
            Range = range;
            Space = space;
        }
        /// <summary>
        /// Gets or sets the intensity range.
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
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Equalize(range.Min, range.Max, 256);
        }
        #endregion
    }
}
