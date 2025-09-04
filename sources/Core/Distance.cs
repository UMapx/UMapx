using System;

namespace UMapx.Core
{
    /// <summary>
    /// Used to calculate distances.
    /// </summary>
    public static class Distance
    {
        #region Euclidean distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Euclidean(this float[] p, float[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += (float)Math.Pow(p[k] - q[k], 2);
            }

            return (float)Math.Sqrt(sum);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Euclidean(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Euclidean(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Euclidean(this Complex32[] p, Complex32[] q)
        {
            Complex32 sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Maths.Pow(p[k] - q[k], 2);
            }

            return Maths.Sqrt(sum);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Euclidean(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Euclidean(u);
            }

            return v;
        }
        #endregion

        #region Chebyshev distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Chebyshev(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Chebyshev(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Chebyshev(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Chebyshev(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            float max = Maths.Abs(p[0] - q[0]);
            float tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = Maths.Abs(p[k] - q[k]);
                max = tmp > max ? tmp : max;
            }

            return max;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Chebyshev(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Chebyshev(u);
            }

            return v;
        }
        #endregion

        #region Manhattan distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Manhattan(this float[] p, float[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Math.Abs(p[k] - q[k]);
            }

            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Manhattan(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Manhattan(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Manhattan(this Complex32[] p, Complex32[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += Maths.Abs(p[k] - q[k]);
            }

            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Manhattan(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Manhattan(u);
            }

            return v;
        }
        #endregion

        #region Angular distance
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Angular(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Angular(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Angular(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Angular(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            Complex32 s = 0;
            Complex32 x = 0;
            Complex32 y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            Complex32 den = Maths.Sqrt(x) * Maths.Sqrt(y);
            Complex32 similarity = s == 0 ? 1.0f : 1.0f - (s / den);

            return Maths.Acos(similarity);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Angular(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Angular(u);
            }

            return v;
        }
        #endregion

        #region Bray-Curtis distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float BrayCurtis(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] BrayCurtis(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.BrayCurtis(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 BrayCurtis(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            Complex32 x = 0;
            Complex32 y = 0;

            for (int i = 0; i < n; i++)
            {
                y += Maths.Abs(p[i] - q[i]);
                x += Maths.Abs(p[i] + q[i]);
            }

            return y / x;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] BrayCurtis(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.BrayCurtis(u);
            }

            return v;
        }
        #endregion

        #region Canberra distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Canberra(this float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Math.Abs(p[i] - q[i]) / (Math.Abs(p[i]) + Math.Abs(q[i]));
            }
            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Canberra(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Canberra(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Canberra(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            Complex32 sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Maths.Abs(p[i] - q[i]) / (Maths.Abs(p[i]) + Maths.Abs(q[i]));
            }
            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Canberra(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Canberra(u);
            }

            return v;
        }
        #endregion

        #region Dice distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Dice(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Dice(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Dice(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Dice(this Complex32[] p, Complex32[] q)
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

            return (tf + ft) / (float)(2.0f * tt + ft + tf);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Dice(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Dice(u);
            }

            return v;
        }
        #endregion

        #region Hellinger distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Hellinger(this float[] p, float[] q)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += (float)Math.Pow(Math.Sqrt(p[i]) - Math.Sqrt(q[i]), 2);
            }

            return sum / Maths.Sqrt2;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Hellinger(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Hellinger(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Hellinger(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            Complex32 sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Maths.Pow(Maths.Sqrt(p[i]) - Maths.Sqrt(q[i]), 2);
            }

            return sum / Maths.Sqrt2;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Hellinger(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Hellinger(u);
            }

            return v;
        }
        #endregion

        #region Jaccard distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Jaccard(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Jaccard(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Jaccard(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Jaccard(this Complex32[] p, Complex32[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Jaccard(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Jaccard(u);
            }

            return v;
        }
        #endregion

        #region Kulczynski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Kulczynski(this float[] p, float[] q)
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

            float num = tf + ft - tt + n;
            float den = ft + tf + n;
            return num / den;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Kulczynski(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Kulczynski(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Kulczynski(this Complex32[] p, Complex32[] q)
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

            float num = tf + ft - tt + n;
            float den = ft + tf + n;
            return num / den;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Kulczynski(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Kulczynski(u);
            }

            return v;
        }
        #endregion

        #region Minkowski distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <param name="order">Order</param>
        /// <returns>Value</returns>
        public static float Minkowski(this float[] p, float[] q, float order)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += (float)Math.Pow(Math.Abs(p[i] - q[i]), order);
            }
            return (float)Math.Pow(sum, 1 / order);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <param name="order">Order</param>
        /// <returns>Vector</returns>
        public static float[] Minkowski(this float[,] p, float[,] q, float order)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Minkowski(u, order);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <param name="order">Order</param>
        /// <returns>Value</returns>
        public static Complex32 Minkowski(this Complex32[] p, Complex32[] q, float order)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += Maths.Pow(Maths.Abs(p[i] - q[i]), order);
            }
            return Maths.Pow(sum, 1 / order);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <param name="order">Order</param>
        /// <returns>Vector</returns>
        public static Complex32[] Minkowski(this Complex32[,] p, Complex32[,] q, float order)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Minkowski(u, order);
            }

            return v;
        }
        #endregion

        #region Russel-Rao distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float RusselRao(this float[] p, float[] q)
        {
            int n = p.Length;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (n - tt) / (float)n;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] RusselRao(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.RusselRao(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 RusselRao(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            int tt = 0;

            for (int i = 0; i < n; i++)
            {
                if (p[i] != 0 && q[i] != 0) tt++;
            }

            return (n - tt) / (float)n;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] RusselRao(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.RusselRao(u);
            }

            return v;
        }
        #endregion

        #region Sokal-Michener distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float SokalMichener(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] SokalMichener(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SokalMichener(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 SokalMichener(this Complex32[] p, Complex32[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] SokalMichener(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SokalMichener(u);
            }

            return v;
        }
        #endregion

        #region Sokal-Sneath distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float SokalSneath(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] SokalSneath(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SokalSneath(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 SokalSneath(this Complex32[] p, Complex32[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] SokalSneath(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SokalSneath(u);
            }

            return v;
        }
        #endregion

        #region Yule distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float Yule(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Yule(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Yule(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Yule(this Complex32[] p, Complex32[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Yule(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Yule(u);
            }

            return v;
        }
        #endregion

        #region Square-Euclidian distance
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static float SquareEuclidian(this float[] p, float[] q)
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
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] SquareEuclidian(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SquareEuclidian(u);
            }

            return v;
        }
        /// <summary>
        /// Returns distance value.
        /// </summary>
        /// <param name="p">Array</param>
        /// <param name="q">Array</param>
        /// <returns>Value</returns>
        public static Complex32 SquareEuclidian(this Complex32[] p, Complex32[] q)
        {
            int n = p.Length;
            Complex32 sum = 0.0f;
            Complex32 u;

            for (int i = 0; i < n; i++)
            {
                u = p[i] - q[i];
                sum += u * u;
            }

            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] SquareEuclidian(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.SquareEuclidian(u);
            }

            return v;
        }
        #endregion

        #region Cosine distance
        /// <summary>
        /// Returns similarity function of two vectors.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value</returns>
        public static float Cosine(this float[] p, float[] b)
        {
            int length = p.Length;
            float A = Matrix.Abs(p, false);
            float B = Matrix.Abs(b, false);
            float s = 0;

            for (int i = 0; i < length; i++)
                s += p[i] * b[i];

            return s / (A * B);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Cosine(this float[,] p, float[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];
                float[] u = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Cosine(u);
            }

            return v;
        }
        /// <summary>
        /// Returns similarity function of two vectors.
        /// </summary>
        /// <param name="p">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Value</returns>
        public static Complex32 Cosine(this Complex32[] p, Complex32[] b)
        {
            int length = p.Length;
            Complex32 A = Matrix.Abs(p, false);
            Complex32 B = Matrix.Abs(b, false);
            Complex32 s = 0;

            for (int i = 0; i < length; i++)
                s += p[i] * b[i];

            return s / (A * B);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Cosine(this Complex32[,] p, Complex32[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];
                Complex32[] u = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = p[i, j];
                    u[j] = q[i, j];
                }

                v[i] = t.Cosine(u);
            }

            return v;
        }
        #endregion
    }
}
