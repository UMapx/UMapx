using System;

namespace UMapx.Imaging
{
    /// <summary>
    /// Using for epsilon extensions.
    /// </summary>
    internal static class EpsilonExtensions
    {
        /// <summary>
        /// Epsilon equals method.
        /// </summary>
        /// <param name="f">Value</param>
        /// <param name="other">Value</param>
        /// <param name="epsilon">Value</param>
        /// <returns>True or false</returns>
        public static bool EpsilonEquals(this float f, float other, float epsilon)
        {
            return Math.Abs(f - other) <= epsilon;
        }
    }
}