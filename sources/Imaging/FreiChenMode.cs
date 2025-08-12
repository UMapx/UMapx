using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Frei-Chen convolution filter mode.
    /// </summary>
    [Serializable]
    public enum FreiChenMode
    {
        /// <summary>
        /// Gradient components (F1, F2).
        /// </summary>
        Gradient,
        /// <summary>
        /// Edge components (F1-F4).
        /// </summary>
        Edge,
        /// <summary>
        /// Line components (F5-F8).
        /// </summary>
        Line,
        /// <summary>
        /// Average component (F9).
        /// </summary>
        Average
    }
}
