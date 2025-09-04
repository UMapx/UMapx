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
        public static ComplexF Euclidean(this ComplexF[] p, ComplexF[] q)
        {
            ComplexF sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += MathsF.Pow(p[k] - q[k], 2);
            }

            return MathsF.Sqrt(sum);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Euclidean(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Chebyshev(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            float max = MathsF.Abs(p[0] - q[0]);
            float tmp;

            for (int k = 1; k < n; k++)
            {
                tmp = MathsF.Abs(p[k] - q[k]);
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
        public static ComplexF[] Chebyshev(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Manhattan(this ComplexF[] p, ComplexF[] q)
        {
            float sum = 0;
            int n = p.Length;

            for (int k = 0; k < n; k++)
            {
                sum += MathsF.Abs(p[k] - q[k]);
            }

            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Manhattan(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Angular(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            ComplexF s = 0;
            ComplexF x = 0;
            ComplexF y = 0;

            for (int i = 0; i < n; i++)
            {
                s += p[i] * q[i];
                x += p[i] * p[i];
                y += q[i] * q[i];
            }

            ComplexF den = MathsF.Sqrt(x) * MathsF.Sqrt(y);
            ComplexF similarity = s == 0 ? 1.0f : 1.0f - (s / den);

            return MathsF.Acos(similarity);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Angular(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF BrayCurtis(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            ComplexF x = 0;
            ComplexF y = 0;

            for (int i = 0; i < n; i++)
            {
                y += MathsF.Abs(p[i] - q[i]);
                x += MathsF.Abs(p[i] + q[i]);
            }

            return y / x;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] BrayCurtis(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Canberra(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            ComplexF sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += MathsF.Abs(p[i] - q[i]) / (MathsF.Abs(p[i]) + MathsF.Abs(q[i]));
            }
            return sum;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Canberra(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Dice(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] Dice(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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

            return sum / MathsF.Sqrt2;
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
        public static ComplexF Hellinger(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            ComplexF sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += MathsF.Pow(MathsF.Sqrt(p[i]) - MathsF.Sqrt(q[i]), 2);
            }

            return sum / MathsF.Sqrt2;
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Hellinger(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Jaccard(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] Jaccard(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Kulczynski(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] Kulczynski(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Minkowski(this ComplexF[] p, ComplexF[] q, float order)
        {
            int n = p.Length;
            float sum = 0;

            for (int i = 0; i < n; i++)
            {
                sum += MathsF.Pow(MathsF.Abs(p[i] - q[i]), order);
            }
            return MathsF.Pow(sum, 1 / order);
        }
        /// <summary>
        /// Returns distance value. 
        /// </summary>
        /// <param name="p">Matrix</param>
        /// <param name="q">Matrix</param>
        /// <param name="order">Order</param>
        /// <returns>Vector</returns>
        public static ComplexF[] Minkowski(this ComplexF[,] p, ComplexF[,] q, float order)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF RusselRao(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] RusselRao(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF SokalMichener(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] SokalMichener(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF SokalSneath(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] SokalSneath(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF Yule(this ComplexF[] p, ComplexF[] q)
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
        public static ComplexF[] Yule(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
        public static ComplexF SquareEuclidian(this ComplexF[] p, ComplexF[] q)
        {
            int n = p.Length;
            ComplexF sum = 0.0f;
            ComplexF u;

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
        public static ComplexF[] SquareEuclidian(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
            float A = MatrixF.Abs(p, false);
            float B = MatrixF.Abs(b, false);
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
        public static ComplexF Cosine(this ComplexF[] p, ComplexF[] b)
        {
            int length = p.Length;
            ComplexF A = MatrixF.Abs(p, false);
            ComplexF B = MatrixF.Abs(b, false);
            ComplexF s = 0;

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
        public static ComplexF[] Cosine(this ComplexF[,] p, ComplexF[,] q)
        {
            int r = p.GetLength(0);
            int c = p.GetLength(1);
            ComplexF[] v = new ComplexF[r];

            for (int i = 0; i < r; i++)
            {
                ComplexF[] t = new ComplexF[c];
                ComplexF[] u = new ComplexF[c];

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
