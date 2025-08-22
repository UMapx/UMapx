using System;
using System.Linq;
using UMapx.Core;

namespace UMapx.Visualization
{
    /// <summary>
    /// Used to work with points.
    /// </summary>
    internal static class Points
    {
        #region Points options
        /// <summary>
        /// Maps a horizontal pixel coordinate to the corresponding data-space X value.
        /// </summary>
        /// <param name="a">
        /// Pixel position along X (0 at the left edge, increasing to the right).
        /// </param>
        /// <param name="amin">Data-space minimum (left bound) of the X axis</param>
        /// <param name="amax">Data-space maximum (right bound) of the X axis</param>
        /// <param name="width">Total drawable width in pixels (must be &gt; 0)</param>
        /// <returns>
        /// The data-space value X that corresponds to pixel <paramref name="a"/>.
        /// </returns>
        /// <remarks>
        /// Assumes a linear mapping from pixel space [0, <paramref name="width"/>]
        /// to data space [<paramref name="amin"/>, <paramref name="amax"/>].
        /// </remarks>
        public static float Point2X(float a, float amin, float amax, float width)
        {
            float dx = (amax - amin) / width;
            return a * dx + amin;
        }
        /// <summary>
        /// Maps a vertical pixel coordinate to the corresponding data-space Y value.
        /// </summary>
        /// <param name="a">
        /// Pixel position along Y (0 at the top edge, increasing downward — screen coordinates).
        /// </param>
        /// <param name="amin">Data-space minimum (bottom bound) of the Y axis</param>
        /// <param name="amax">Data-space maximum (top bound) of the Y axis</param>
        /// <param name="height">Total drawable height in pixels (must be &gt; 0)</param>
        /// <returns>
        /// The data-space value Y that corresponds to pixel <paramref name="a"/>.
        /// </returns>
        /// <remarks>
        /// Uses an inverted Y axis (top-left origin), common for raster graphics:
        /// pixel 0 maps to <paramref name="amax"/>, pixel <paramref name="height"/>
        /// maps to <paramref name="amin"/>.
        /// </remarks>
        public static float Point2Y(float a, float amin, float amax, float height)
        {
            float dy = (amax - amin) / height;
            return (height - a) * dy + amin;
        }
        /// <summary>
        /// Maps a data-space X value to the corresponding horizontal pixel coordinate.
        /// </summary>
        /// <param name="a">Data-space X value to convert.</param>
        /// <param name="amin">Data-space minimum (left bound) of the X axis</param>
        /// <param name="amax">Data-space maximum (right bound) of the X axis</param>
        /// <param name="width">Total drawable width in pixels (must be &gt; 0)</param>
        /// <returns>
        /// Pixel position along X in the range [0, <paramref name="width"/>].
        /// </returns>
        public static float X2Point(float a, float amin, float amax, float width)
        {
            float dx = (amax - amin) / width;
            return (a - amin) / dx;
        }
        /// <summary>
        /// Maps a data-space Y value to the corresponding vertical pixel coordinate.
        /// </summary>
        /// <param name="a">Data-space Y value to convert</param>
        /// <param name="amin">Data-space minimum (bottom bound) of the Y axis</param>
        /// <param name="amax">Data-space maximum (top bound) of the Y axis</param>
        /// <param name="height">Total drawable height in pixels (must be &gt; 0)</param>
        /// <returns>
        /// Pixel position along Y in the range [0, <paramref name="height"/>],
        /// where 0 is the top edge.
        /// </returns>
        /// <remarks>
        /// Inverts the Y axis (top-left origin) to match typical screen coordinates.
        /// </remarks>
        public static float Y2Point(float a, float amin, float amax, float height)
        {
            float dy = (amax - amin) / height;
            return height - (a - amin) / dy;
        }
        /// <summary>
        /// Tests whether a scalar value is singular (not a finite real number).
        /// </summary>
        /// <param name="a">Value to test</param>
        /// <returns>
        /// <see langword="true"/> if <paramref name="a"/> is <see cref="float.NaN"/>,
        /// <see cref="float.PositiveInfinity"/>, or <see cref="float.NegativeInfinity"/>;
        /// otherwise <see langword="false"/>.
        /// </returns>
        /// <remarks>
        /// This is a convenience wrapper for detecting non-finite values.
        /// </remarks>
        public static bool IsSingularPoint(float a)
        {
            return Maths.IsSquare(a);
        }
        /// <summary>
        /// Clips a value to a half-open range, returning a sentinel just outside if out of bounds.
        /// </summary>
        /// <param name="a">Value to clip</param>
        /// <param name="amin">Inclusive lower bound</param>
        /// <param name="amax">Inclusive upper bound</param>
        /// <returns>
        /// <para>
        /// If <paramref name="a"/> is within [<paramref name="amin"/>, <paramref name="amax"/>],
        /// returns <paramref name="a"/>.
        /// </para>
        /// <para>
        /// If below the range, returns <paramref name="amin"/> - 1; if above, returns
        /// <paramref name="amax"/> + 1. This sentinel behavior is useful to mark points
        /// as “just outside” the drawable area without losing ordering information.
        /// </para>
        /// </returns>
        public static float ClipPoint(float a, float amin, float amax)
        {
            if (a < amin)
            {
                return amin - 1;
            }
            else if (a > amax)
            {
                return amax + 1;
            }
            return a;
        }
        /// <summary>
        /// Returns the minimum finite value in an array, ignoring singular entries.
        /// </summary>
        /// <param name="v">Input array</param>
        /// <returns>
        /// The minimum finite value if present; otherwise <see langword="null"/> when all
        /// entries are singular (<see cref="float.NaN"/> or infinities) or the array is empty.
        /// </returns>
        /// <remarks>
        /// Singular values are filtered using <c>Maths.IsSingular</c>.
        /// </remarks>
        public static float? GetMin(this float[] v)
        {
            var w = v.Where(x => !Maths.IsSingular(x)).ToArray();
            if (w.Length > 0) return w.Min();
            return null;
        }
        /// <summary>
        /// Returns the maximum finite value in an array, ignoring singular entries.
        /// </summary>
        /// <param name="v">Input array</param>
        /// <returns>
        /// The maximum finite value if present; otherwise <see langword="null"/> when all
        /// entries are singular (<see cref="float.NaN"/> or infinities) or the array is empty.
        /// </returns>
        /// <remarks>
        /// Singular values are filtered using <c>Maths.IsSingular</c>.
        /// </remarks>
        public static float? GetMax(this float[] v)
        {
            var w = v.Where(x => !Maths.IsSingular(x)).ToArray();
            if (w.Length > 0) return w.Max();
            return null;
        }
        /// <summary>
        /// Generates evenly spaced tick marks between <paramref name="min"/> and <paramref name="max"/> (inclusive).
        /// </summary>
        /// <param name="min">Lower bound of the axis</param>
        /// <param name="max">Upper bound of the axis</param>
        /// <param name="points">
        /// Number of intervals to split the range into. The resulting array has length <c>points + 1</c>.
        /// </param>
        /// <returns>
        /// An array of size <c>points + 1</c> containing tick values from
        /// <paramref name="min"/> to <paramref name="max"/> with a uniform step.
        /// </returns>
        /// <remarks>
        /// The step is rounded to 2 decimal places, which may cause the last tick
        /// to deviate slightly from <paramref name="max"/> for certain ranges.
        /// </remarks>
        public static float[] GetPoints(float min, float max, int points)
        {
            float delta = (float)Math.Round((max - min) / points, 2);
            float[] marks = new float[points + 1];
            float c = min;

            for (int i = 0; i < points + 1; i++)
            {
                marks[i] = c;
                c += delta;
            }
            return marks;
        }
        #endregion
    }
}
