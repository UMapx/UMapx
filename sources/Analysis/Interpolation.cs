using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements interpolation.
    /// <remarks>
    /// This class is a solution to the problem of finding an intermediate value of the function F(x).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Interpolation
    {
        #region Private data
        private InterpolationMethod method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements interpolation.
        /// </summary>
        /// <param name="method">Interpolation method</param>
        public Interpolation(InterpolationMethod method = InterpolationMethod.Lagrange)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the interpolation method.
        /// </summary>
        public InterpolationMethod MethodType
        {
            get
            {
                return this.method;
            }
            set
            {
                this.method = value;
            }
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// <remarks>
        /// In this case, only bilinear interpolation is used.
        /// </remarks>
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="z">Function matrix</param>
        /// <param name="xl">The value of the first argument to calculate</param>
        /// <param name="yl">The value of the second argument to calculate</param>
        /// <returns>Value</returns>
        public float Compute(float[] x, float[] y, float[,] z, float xl, float yl)
        {
            return bilinear(x, y, z, xl, yl);
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Value</returns>
        public float Compute(float[] x, float[] y, float xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case InterpolationMethod.Newton:
                    return Interpolation.newto(x, y, xl);

                case InterpolationMethod.Barycentric:
                    return Interpolation.baryc(x, y, xl);

                case InterpolationMethod.Lagrange:
                    return Interpolation.lagra(x, y, xl);

                case InterpolationMethod.Linear:
                default:
                    return Interpolation.linear(x, y, xl);
            }
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Complex number</returns>
        public Complex32 Compute(Complex32[] x, Complex32[] y, Complex32 xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case InterpolationMethod.Newton:
                    return Interpolation.newto(x, y, xl);

                case InterpolationMethod.Barycentric:
                    return Interpolation.baryc(x, y, xl);

                case InterpolationMethod.Lagrange:
                default:
                    return Interpolation.lagra(x, y, xl);

                case InterpolationMethod.Linear:
                    throw new NotSupportedException("Linear interpolation is not defined for complex-valued functions");
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xl"></param>
        /// <returns></returns>
        private static float linear(float[] x, float[] y, float xl)
        {
            int n = x?.Length ?? 0;
            if (n == 0 || y == null || y.Length != n) throw new ArgumentException();
            if (xl <= x[0]) return y[0];
            if (xl >= x[n - 1]) return y[n - 1];

            int lo = 0, hi = n - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) >> 1;
                if (xl >= x[mid]) lo = mid; else hi = mid;
            }
            float t = (xl - x[lo]) / (x[lo + 1] - x[lo]);
            return y[lo] + t * (y[lo + 1] - y[lo]);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        /// <param name="xval"></param>
        /// <param name="yval"></param>
        /// <returns></returns>
        private static float bilinear(float[] x, float[] y, float[,] z, float xval, float yval)
        {
            int nx = x?.Length ?? 0, ny = y?.Length ?? 0;
            if (nx < 2 || ny < 2 || z == null || z.GetLength(0) != nx || z.GetLength(1) != ny)
                throw new ArgumentException();

            if (xval <= x[0]) return (yval <= y[0]) ? z[0, 0] : (yval >= y[ny - 1] ? z[0, ny - 1] : z[0, LowerIndex(y, yval)]);
            if (xval >= x[nx - 1]) return (yval <= y[0]) ? z[nx - 1, 0] : (yval >= y[ny - 1] ? z[nx - 1, ny - 1] : z[nx - 1, LowerIndex(y, yval)]);

            int ix = LowerIndex(x, xval); // x[ix] <= xval < x[ix+1]
            int iy = LowerIndex(y, yval); // y[iy] <= yval < y[iy+1]

            float x0 = x[ix], x1 = x[ix + 1];
            float y0 = y[iy], y1 = y[iy + 1];
            float tx = (xval - x0) / (x1 - x0);
            float ty = (yval - y0) / (y1 - y0);

            float z00 = z[ix, iy], z10 = z[ix + 1, iy], z01 = z[ix, iy + 1], z11 = z[ix + 1, iy + 1];
            float z0 = z00 + tx * (z10 - z00);
            float z1 = z01 + tx * (z11 - z01);
            return z0 + ty * (z1 - z0);

            static int LowerIndex(float[] a, float v)
            {
                int lo = 0, hi = a.Length - 1;
                if (v <= a[0]) return 0;
                if (v >= a[hi]) return hi - 1;
                while (hi - lo > 1)
                {
                    int mid = (lo + hi) >> 1;
                    if (v >= a[mid]) lo = mid; else hi = mid;
                }
                return lo;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static float lagra(float[] x, float[] y, float xval)
        {
            float yval = 0.0f;
            float Products = y[0];
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                Products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        Products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += Products;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static float newto(float[] x, float[] y, float xval)
        {
            int n = x.Length;
            var a = (float[])y.Clone();
            for (int i = 1; i < n; i++)
                for (int j = n - 1; j >= i; j--)
                    a[j] = (a[j] - a[j - 1]) / (x[j] - x[j - i]);

            float res = a[n - 1];
            for (int i = n - 2; i >= 0; i--)
                res = a[i] + (xval - x[i]) * res;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static float baryc(float[] x, float[] y, float xval)
        {
            int n = x.Length;
            for (int i = 0; i < n; i++) if (xval == x[i]) return y[i];

            var w = new float[n];
            for (int i = 0; i < n; i++)
            {
                float prod = 1f;
                for (int j = 0; j < n; j++) if (i != j) prod *= (x[i] - x[j]);
                if (prod == 0f) throw new ArgumentException("Duplicate nodes");
                w[i] = 1f / prod;
            }
            float num = 0f, den = 0f;
            for (int i = 0; i < n; i++)
            {
                float d = xval - x[i];
                float wi = w[i] / d;
                num += wi * y[i];
                den += wi;
            }
            return num / den;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex32 lagra(Complex32[] x, Complex32[] y, Complex32 xval)
        {
            Complex32 yval = 0.0;
            Complex32 Products = y[0];
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                Products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        Products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += Products;
            }
            return yval;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex32 newto(Complex32[] x, Complex32[] y, Complex32 xval)
        {
            int n = x.Length;
            var a = (Complex32[])y.Clone();
            for (int i = 1; i < n; i++)
                for (int j = n - 1; j >= i; j--)
                    a[j] = (a[j] - a[j - 1]) / (x[j] - x[j - i]);

            Complex32 res = a[n - 1];
            for (int i = n - 2; i >= 0; i--)
                res = a[i] + (xval - x[i]) * res;
            return res;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="xval"></param>
        /// <returns></returns>
        private static Complex32 baryc(Complex32[] x, Complex32[] y, Complex32 xval)
        {
            int n = x.Length;
            for (int i = 0; i < n; i++) if (xval == x[i]) return y[i];

            var w = new Complex32[n];
            for (int i = 0; i < n; i++)
            {
                Complex32 prod = 1f;
                for (int j = 0; j < n; j++) if (i != j) prod *= (x[i] - x[j]);
                if (prod == 0f) throw new ArgumentException("Duplicate nodes");
                w[i] = 1f / prod;
            }
            Complex32 num = 0f, den = 0f;
            for (int i = 0; i < n; i++)
            {
                Complex32 d = xval - x[i];
                Complex32 wi = w[i] / d;
                num += wi * y[i];
                den += wi;
            }
            return num / den;
        }
        #endregion
    }
}
