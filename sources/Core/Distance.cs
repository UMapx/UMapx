using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to calculate distances.
    /// </summary>
    public static class Distance
    {
        #region Euclidean distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Euclidean(double[] p, double[] q)
        {
            double sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Pow(p[k] - q[k], 2);
            }

            return Math.Sqrt(sum);
        }
        #endregion

        #region Chebyshev distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Chebyshev(double[] p, double[] q)
        {
            int n = p.Length;
            double max = Math.Abs(p[0] - q[0]);
            double tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = Math.Abs(p[k] - q[k]);
                max = tmp > max ? tmp : max;
            }

            return max;
        }
        #endregion

        #region Manhattan distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Manhattan(double[] p, double[] q)
        {
            double sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Abs(p[k] - q[k]);
            }

            return sum;
        }
        #endregion

        #region Angular distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Angular(double[] p, double[] q)
        {
            int n = p.Length;
            double s = 0;
            double x = 0;
            double y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            double den = Math.Sqrt(x) * Math.Sqrt(y);
            double similarity = s == 0 ? 1.0 : 1.0 - (s / den);

            return Math.Acos(similarity);
        }
        #endregion

        #region Bray-Curtis distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double BrayCurtis(double[] p, double[] q)
        {
            int n = p.Length;
            double x = 0;
            double y = 0;

            for (int i = 0; i < n; i++)
            {
                y += Math.Abs(p[i] - q[i]);
                x += Math.Abs(p[i] + q[i]);
            }

            return y / x;
        }
        #endregion

        #region Canberra distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Canberra(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Abs(p[i] - q[i]) / (Math.Abs(p[i]) + Math.Abs(q[i]));
            }
            return sum;
        }
        #endregion

        #region Dice distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Dice(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (tf + ft) / (double)(2 * tt + ft + tf);
        }
        #endregion

        #region Hellinger distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Hellinger(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(Math.Sqrt(p[i]) - Math.Sqrt(q[i]), 2);
            }

            return sum / Math.Sqrt(2);
        }
        #endregion

        #region Jaccard distance
        /// <summary>
        /// Returns distance value".
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Jaccard(double[] p, double[] q)
        {
            int n = p.Length;
            int inter = 0;
            int union = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 || q[i] != 0)
                {
                    if (p[i] == q[i])
                        inter++;
                    union++;
                }
            }

            return (union == 0) ? 0 : 1.0 - (inter / (double)union);
        }
        #endregion

        #region Kulczynski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Kulczynski(double[] p, double[] q)
        {
            // TODO: Rewrite the integer dissimilarities (Yule, Russel-Rao,...)
            // using generics
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            double num = tf + ft - tt + n;
            double den = ft + tf + n;
            return num / den;
        }
        #endregion

        #region Minkowski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <param name="order">Order</param>
        /// <returns>Double precision floating point number</returns>
        public static double Minkowski(double[] p, double[] q, double order)
        {
            int n = p.Length;
            double sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Pow(Math.Abs(p[i] - q[i]), order);
            }
            return Math.Pow(sum, 1 / order);
        }
        #endregion

        #region Russel-Rao distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double RusselRao(double[] p, double[] q)
        {
            int n = p.Length;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (n - tt) / (double)(n);
        }
        #endregion

        #region Sokal-Michener distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double SokalMichener(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] == 1 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] == 1) ft++;
                if (p[i] == 1 && q[i] == 1) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            int r = 2 * (tf + ft);
            return r / (double)(ff + tt + r);
        }
        #endregion

        #region Sokal-Sneath distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double SokalSneath(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            int r = 2 * (tf + ft);
            return r / (double)(tt + r);
        }
        #endregion

        #region Yule distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Yule(double[] p, double[] q)
        {
            int n = p.Length;
            int tf = 0;
            int ft = 0;
            int tt = 0;
            int ff = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] == 0) tf++;
                if (p[i] == 0 && q[i] != 0) ft++;
                if (p[i] != 0 && q[i] != 0) tt++;
                if (p[i] == 0 && q[i] == 0) ff++;
            }

            double r = 2 * (tf + ft);
            return r / (tt + ff + r / 2);
        }
        #endregion

        #region Square-Euclidian distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double SquareEuclidian(double[] p, double[] q)
        {
            int n = p.Length;
            double sum = 0.0;
            double u;

            for (int i = 0; i < n; i++)
            {
                u = p[i] - q[i];
                sum += u * u;
            }

            return sum;
        }
        #endregion
    }
}
