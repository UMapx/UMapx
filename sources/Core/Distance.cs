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
        /// <returns>float precision floating point number</returns>
        public static float Euclidean(float[] p, float[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += (float)Math.Pow(p[k] - q[k], 2);
            }

            return (float)Math.Sqrt(sum);
        }
        #endregion

        #region Chebyshev distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Chebyshev(float[] p, float[] q)
        {
            int n = p.Length;
            float max = Math.Abs(p[0] - q[0]);
            float tmp;

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
        /// <returns>float precision floating point number</returns>
        public static float Manhattan(float[] p, float[] q)
        {
            float sum = 0;
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
        /// <returns>float precision floating point number</returns>
        public static float Angular(float[] p, float[] q)
        {
            int n = p.Length;
            float s = 0;
            float x = 0;
            float y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            float den = (float)Math.Sqrt(x) * (float)Math.Sqrt(y);
            float similarity = s == 0 ? 1.0f : 1.0f - (s / den);

            return (float)Math.Acos(similarity);
        }
        #endregion

        #region Bray-Curtis distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float BrayCurtis(float[] p, float[] q)
        {
            int n = p.Length;
            float x = 0;
            float y = 0;

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
        /// <returns>float precision floating point number</returns>
        public static float Canberra(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

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
        /// <returns>float precision floating point number</returns>
        public static float Dice(float[] p, float[] q)
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

            return (tf + ft) / (float)(2 * tt + ft + tf);
        }
        #endregion

        #region Hellinger distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Hellinger(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += (float)Math.Pow(Math.Sqrt(p[i]) - Math.Sqrt(q[i]), 2);
            }

            return sum / Maths.Sqrt2;
        }
        #endregion

        #region Jaccard distance
        /// <summary>
        /// Returns distance value".
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Jaccard(float[] p, float[] q)
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

            return (union == 0) ? 0 : 1.0f - (inter / (float)union);
        }
        #endregion

        #region Kulczynski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Kulczynski(float[] p, float[] q)
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

            float num = tf + ft - tt + n;
            float den = ft + tf + n;
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
        /// <returns>float precision floating point number</returns>
        public static float Minkowski(float[] p, float[] q, float order)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += (float)Math.Pow(Math.Abs(p[i] - q[i]), order);
            }
            return (float)Math.Pow(sum, 1 / order);
        }
        #endregion

        #region Russel-Rao distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float RusselRao(float[] p, float[] q)
        {
            int n = p.Length;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (n - tt) / (float)(n);
        }
        #endregion

        #region Sokal-Michener distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float SokalMichener(float[] p, float[] q)
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
            return r / (float)(ff + tt + r);
        }
        #endregion

        #region Sokal-Sneath distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float SokalSneath(float[] p, float[] q)
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
            return r / (float)(tt + r);
        }
        #endregion

        #region Yule distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Yule(float[] p, float[] q)
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

            float r = 2 * (tf + ft);
            return r / (tt + ff + r / 2);
        }
        #endregion

        #region Square-Euclidian distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float SquareEuclidian(float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0.0f;
            float u;

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
