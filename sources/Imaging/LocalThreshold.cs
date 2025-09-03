using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Bradley local threshold filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LocalThreshold : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private float difference;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(int radius, Space space, float difference = 0.15f)
        {
            gb = new BoxBlur(radius);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(int width, int height, Space space, float difference = 0.15f)
        {
            gb = new BoxBlur(width, height);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(SizeInt size, Space space, float difference = 0.15f)
        {
            gb = new BoxBlur(size);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Gets or sets the difference [0, 1].
        /// </summary>
        public float Difference
        {
            get
            {
                return difference;
            }
            set
            {
                difference = MathF.Float(value);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bradley(this.difference, 256);
        }
        #endregion
    }
}
