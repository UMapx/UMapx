using System;

namespace UMapx.Core
{
    /// <summary>
    /// Defines a morphology mode.
    /// </summary>
    [Serializable]
    public enum MorphologyMode
    {
        /// <summary>
        /// Erosion.
        /// </summary>
        Erosion,
        /// <summary>
        /// Median.
        /// </summary>
        Median,
        /// <summary>
        /// Dilatation.
        /// </summary>
        Dilatation
    }
}
