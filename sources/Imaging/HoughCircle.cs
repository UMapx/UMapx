using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Hough circle.
    /// </summary>
    [Serializable]
    public struct HoughCircle : IComparable
    {
        #region Struct components
        /// <summary>
        /// Coordinate X.
        /// </summary>
        public readonly int X;
        /// <summary>
        /// Coordinate Y.
        /// </summary>
        public readonly int Y;
        /// <summary>
        /// Radius.
        /// </summary>
        public readonly int Radius;
        /// <summary>
        /// Absolute line intensity (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Relative line intensity (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Initializes the Hough circle.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="radius">Radius</param>
        /// <param name="intensity">Absolute line intensity (0, +inf)</param>
        /// <param name="relativeIntensity">Relative line intensity (0, 1]</param>
        public HoughCircle(int x, int y, int radius, short intensity, double relativeIntensity)
        {
            X = x;
            Y = y;
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
            return (-Intensity.CompareTo(((HoughCircle)value).Intensity));
        }
        #endregion
    }
}
