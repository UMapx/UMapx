using System;

namespace UMapx.Visualization
{
    /// <summary>
    /// Used to work with points.
    /// </summary>
    internal static class Points
    {
        #region Points options
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="amin"></param>
        /// <param name="amax"></param>
        /// <param name="width"></param>
        /// <returns></returns>
        public static float Point2X(float a, float amin, float amax, float width)
        {
            float dx = (amax - amin) / width;
            return a * dx + amin;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="amin"></param>
        /// <param name="amax"></param>
        /// <param name="height"></param>
        /// <returns></returns>
        public static float Point2Y(float a, float amin, float amax, float height)
        {
            float dy = (amax - amin) / height;
            return (height - a) * dy + amin;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="amin"></param>
        /// <param name="amax"></param>
        /// <param name="width"></param>
        /// <returns></returns>
        public static float X2Point(float a, float amin, float amax, float width)
        {
            float dx = (amax - amin) / width;
            return (a - amin) / dx;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="amin"></param>
        /// <param name="amax"></param>
        /// <param name="height"></param>
        /// <returns></returns>
        public static float Y2Point(float a, float amin, float amax, float height)
        {
            float dy = (amax - amin) / height;
            return height - (a - amin) / dy;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static bool IsSingularPoint(float a)
        {
            if (float.IsNaN(a))
            {
                return true;
            }
            else if (float.IsNegativeInfinity(a) || float.IsPositiveInfinity(a))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="amin"></param>
        /// <param name="amax"></param>
        /// <returns></returns>
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
        /// 
        /// </summary>
        /// <param name="min"></param>
        /// <param name="max"></param>
        /// <param name="points"></param>
        /// <returns></returns>
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
