using System;
using SkiaDrawing;

namespace UMapx.Video
{
    /// <summary>
    ///  Capabilities of video device such as frame size and frame rate.
    /// </summary>
    public class VideoCapabilities
    {
        #region Properties

        /// <summary>
        /// Frame size supported by video device.
        /// </summary>
        public Size FrameSize { get; protected set; }

        /// <summary>
        /// Average frame rate of video device for corresponding <see cref="FrameSize">frame size</see>.
        /// </summary>
        public int AverageFrameRate { get; protected set; }

        /// <summary>
        /// Maximum frame rate of video device for corresponding <see cref="FrameSize">frame size</see>.
        /// </summary>
        public int MaximumFrameRate { get; protected set; }

        /// <summary>
        /// Number of bits per pixel provided by the camera.
        /// </summary>
        public int BitCount { get; protected set; }

        #endregion

        #region Constructors

        /// <summary>
        /// Initializes video capabilities.
        /// </summary>
        protected VideoCapabilities() { }

        /// <summary>
        /// Initializes video capabilities.
        /// </summary>
        /// <param name="frameSize">Frame size</param>
        /// <param name="averageFrameRate">Average frame rate</param>
        /// <param name="maximumFrameRate">Maximum frame rate</param>
        /// <param name="bitCount">Bit count</param>
        public VideoCapabilities(Size frameSize, int averageFrameRate, int maximumFrameRate, int bitCount)
        {
            FrameSize = frameSize;
            AverageFrameRate = averageFrameRate;
            MaximumFrameRate = maximumFrameRate;
            BitCount = bitCount;
        }

        #endregion

        #region Overrides

        /// <summary>
        /// Check if the video capability equals to the specified object.
        /// </summary>
        /// 
        /// <param name="obj">Object to compare with.</param>
        /// 
        /// <returns>Returns true if both are equal are equal or false otherwise.</returns>
        /// 
        public override bool Equals(object obj)
        {
            return Equals(obj as VideoCapabilities);
        }

        /// <summary>
        /// Check if two video capabilities are equal.
        /// </summary>
        /// 
        /// <param name="vc2">Second video capability to compare with.</param>
        /// 
        /// <returns>Returns true if both video capabilities are equal or false otherwise.</returns>
        /// 
        public bool Equals(VideoCapabilities vc2)
        {
            if (vc2 is null)
            {
                return false;
            }

            return ((FrameSize == vc2.FrameSize) && (BitCount == vc2.BitCount));
        }

        /// <summary>
        /// Get hash code of the object.
        /// </summary>
        /// 
        /// <returns>Returns hash code ot the object </returns>
        public override int GetHashCode()
        {
            return FrameSize.GetHashCode() ^ BitCount;
        }

        /// <summary>
        /// Equality operator.
        /// </summary>
        /// 
        /// <param name="a">First object to check.</param>
        /// <param name="b">Seconds object to check.</param>
        /// 
        /// <returns>Return true if both objects are equal or false otherwise.</returns>
        public static bool operator ==(VideoCapabilities a, VideoCapabilities b)
        {
            // if both are null, or both are same instance, return true.
            if (object.ReferenceEquals(a, b))
            {
                return true;
            }

            // if one is null, but not both, return false.
            if ((a is null) || (b is null))
            {
                return false;
            }

            return a.Equals(b);
        }

        /// <summary>
        /// Inequality operator.
        /// </summary>
        /// 
        /// <param name="a">First object to check.</param>
        /// <param name="b">Seconds object to check.</param>
        /// 
        /// <returns>Return true if both objects are not equal or false otherwise.</returns>
        public static bool operator !=(VideoCapabilities a, VideoCapabilities b)
        {
            return !(a == b);
        }

        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString()
        {
            return String.Format("{0}x{1}, {2} fps ({3} max fps), {4} bpp",
                FrameSize.Width, FrameSize.Height,
                AverageFrameRate, MaximumFrameRate,
                BitCount);
        }

        #endregion
    }
}
