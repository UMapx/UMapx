using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Hough line.
    /// </summary>
    [Serializable]
    public struct HoughLine : IComparable
    {
        #region Struct components
        /// <summary>
        /// Slope of the line [0, 180).
        /// <remarks>
        /// It is the angle between the polar axis and the radius of the line.
        /// </remarks>
        /// </summary>
        public readonly double Theta;
        /// <summary>
        /// Line distance from the center of the image (-inf, +inf).
        /// <remarks>
        /// A negative radius line means that the line is at the bottom of the polar coordinate system. Therefore
        /// the angle θ should be increased by 180 degrees, and Radius should be positive.
        /// </remarks>
        /// </summary>
        public readonly short Radius;
        /// <summary>
        /// Absolute line intensity (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Relative line intensity (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Initializes the Hough line.
        /// </summary>
        /// <param name="theta">Slope of the line [0, 180)</param>
        /// <param name="radius">Radius (-inf, +inf)</param>
        /// <param name="intensity">Absolute line intensity (0, +inf)</param>
        /// <param name="relativeIntensity">Relative line intensity (0, 1]</param>
        public HoughLine(double theta, short radius, short intensity, double relativeIntensity)
        {
            Theta = theta;
            Radius = radius;
            Intensity = intensity;
            RelativeIntensity = relativeIntensity;
        }
        /// <summary>
        /// Compares object to another instance of this class.
        /// </summary>
        /// <param name="value">Object</param>
        /// <returns>Integer number</returns>
        public int CompareTo(object value)
        {
            return (-Intensity.CompareTo(((HoughLine)value).Intensity));
        }
        #endregion
    }
}
