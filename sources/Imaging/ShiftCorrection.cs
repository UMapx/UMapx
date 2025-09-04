using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the shift correction filter.
    /// </summary>
    [Serializable]
    public class ShiftCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private float offset;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shift correction filter.
        /// </summary>
        /// <param name="offset">Offset (-0.5, 0.5)</param>
        /// <param name="space">Color space</param>
        public ShiftCorrection(float offset, Space space)
        {
            Offset = offset; Space = space;
        }
        /// <summary>
        /// Gets or sets the offset value (-0.5, 0.5).
        /// </summary>
        public float Offset
        {
            get
            {
                return this.offset;
            }
            set
            {
                this.offset = MathsF.Range(value, -0.49998f, 0.49998f);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Shift(this.offset, 256);
        }
        #endregion
    }
}
