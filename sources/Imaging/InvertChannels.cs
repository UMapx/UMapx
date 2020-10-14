using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the channels inversion filter.
    /// </summary>
    [Serializable]
    public class InvertChannels : Correction, IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the channels inversion filter.
        /// </summary>
        public InvertChannels(Space space)
        {
            this.values = Intensity.Invert(256);
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion
    }
}
