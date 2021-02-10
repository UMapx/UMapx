using System;

namespace UMapx.Visualization
{
    /// <summary>
    /// Uses for points operations.
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
        public static double Point2X(double a, double amin, double amax, double width)
        {
            double dx = (amax - amin) / width;
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
        public static double Point2Y(double a, double amin, double amax, double height)
        {
            double dy = (amax - amin) / height;
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
        public static double X2Point(double a, double amin, double amax, double width)
        {
            double dx = (amax - amin) / width;
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
        public static double Y2Point(double a, double amin, double amax, double height)
        {
            double dy = (amax - amin) / height;
            return height - (a - amin) / dy;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="a"></param>
        /// <returns></returns>
        public static bool IsSingularPoint(double a)
        {
            if (double.IsNaN(a))
            {
                return true;
            }
            else if (double.IsNegativeInfinity(a) || double.IsPositiveInfinity(a))
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
        public static double ClipPoint(double a, double amin, double amax)
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
        public static double[] GetPoints(double min, double max, int points)
        {
            double delta = Math.Round((max - min) / points, 2);
            double[] marks = new double[points + 1];
            double c = min;

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
