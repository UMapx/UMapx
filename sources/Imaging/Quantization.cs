using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the quantization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://en.wikipedia.org/wiki/Posterization
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Quantization : Correction, IBitmapFilter
    {
        #region Private data
        private int levels;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the quantization filter.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="space">Color space</param>
        public Quantization(int levels, Space space)
        {
            Levels = levels; Space = space;
        }
        /// <summary>
        /// Initializes the quantization filter.
        /// </summary>
        public Quantization()
        {
            Levels = 4;
        }
        /// <summary>
        /// Gets or sets number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return levels;
            }
            set
            {
                this.levels = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Quantize(this.levels, 256);
        }
        #endregion
    }
}
