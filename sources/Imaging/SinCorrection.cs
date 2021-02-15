using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the sine correction filter.
    /// </summary>
    [Serializable]
    public class SinCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private float delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the sine correction filter.
        /// </summary>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public SinCorrection(float delta, Space space)
        {
            Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the sine correction filter.
        /// </summary>
        public SinCorrection()
        {
            Delta = 20;
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public float Delta
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
            this.values = Intensity.Sin(delta / 2.0f, 256);
        }
        #endregion
    }
}
