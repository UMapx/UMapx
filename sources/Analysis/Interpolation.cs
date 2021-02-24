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
        private Interpolation.Method method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements interpolation.
        /// </summary>
        /// <param name="method">Interpolation method</param>
        public Interpolation(Interpolation.Method method = Method.Lagrange)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the interpolation method.
        /// </summary>
        public Interpolation.Method MethodType
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public float Compute(float[] x, float[] y, float xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case Method.Lagrange:
                    return Interpolation.lagra(x, y, xl);

                case Method.Newton:
                    return Interpolation.newto(x, y, xl);

                case Method.Barycentric:
                    return Interpolation.baryc(x, y, xl);

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
                case Method.Newton:
                    return Interpolation.newto(x, y, xl);

                case Method.Barycentric:
                    return Interpolation.baryc(x, y, xl);

                default:
                    return Interpolation.lagra(x, y, xl);
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
            float yval = 0.0f;
            int length = x.Length - 1;

            for (int i = 0; i < length; i++)
            {
                if (xl >= x[i] && xl < x[i + 1])
                {
                    yval = y[i] + (xl - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
                }
            }
            return yval;
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
            float zval = 0.0f;
            int xlength = x.Length - 1;
            int ylength = y.Length - 1;

            for (int i = 0; i < xlength; i++)
            {
                for (int j = 0; j < ylength; j++)
                {
                    if (xval >= x[i] && xval < x[i + 1] && yval >= y[j] && yval < y[j + 1])
                    {
                        zval = z[i, j] * (x[i + 1] - xval) * (y[j + 1] - yval) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i + 1, j] * (xval - x[i]) * (y[j + 1] - yval) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i, j + 1] * (x[i + 1] - xval) * (yval - y[j]) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]) +
                        z[i + 1, j + 1] * (xval - x[i]) * (yval - y[j]) / (x[i + 1] - x[i]) / (y[j + 1] - y[j]);
                    }
                }
            }
            return zval;
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
            float yval;
            int length = x.Length;
            float[] tarray = new float[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                tarray[i] = y[i];
            }
            for (i = 0; i < length - 1; i++)
            {
                for (j = length - 1; j > i; j--)
                {
                    tarray[j] = (tarray[j - 1] - tarray[j]) / (x[j - 1 - i] - x[j]);
                }
            }
            yval = tarray[length - 1];
            for (i = length - 2; i >= 0; i--)
            {
                yval = tarray[i] + (xval - x[i]) * yval;
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
        private static float baryc(float[] x, float[] y, float xval)
        {
            float product;
            float deltaX;
            float bc1 = 0;
            float bc2 = 0;
            int length = x.Length;
            float[] weights = new float[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                product = 1;
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        product *= (x[i] - x[j]);
                        weights[i] = 1.0f / product;
                    }
                }
            }

            for (i = 0; i < length; i++)
            {
                deltaX = weights[i] / (xval - x[i]);
                bc1 += y[i] * deltaX;
                bc2 += deltaX;
            }
            return bc1 / bc2;
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
            Complex32 yval;
            int length = x.Length;
            Complex32[] tarray = new Complex32[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                tarray[i] = y[i];
            }
            for (i = 0; i < length - 1; i++)
            {
                for (j = length - 1; j > i; j--)
                {
                    tarray[j] = (tarray[j - 1] - tarray[j]) / (x[j - 1 - i] - x[j]);
                }
            }
            yval = tarray[length - 1];
            for (i = length - 2; i >= 0; i--)
            {
                yval = tarray[i] + (xval - x[i]) * yval;
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
        private static Complex32 baryc(Complex32[] x, Complex32[] y, Complex32 xval)
        {
            Complex32 product;
            Complex32 deltaX;
            Complex32 bc1 = 0;
            Complex32 bc2 = 0;
            int length = x.Length;
            Complex32[] weights = new Complex32[length];
            int i, j;

            for (i = 0; i < length; i++)
            {
                product = 1;
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        product *= (x[i] - x[j]);
                        weights[i] = 1.0 / product;
                    }
                }
            }

            for (i = 0; i < length; i++)
            {
                deltaX = weights[i] / (xval - x[i]);
                bc1 += y[i] * deltaX;
                bc2 += deltaX;
            }
            return bc1 / bc2;
        }
        #endregion

        #region Enums
        /// <summary>
        /// Interpolation method.
        /// </summary>
        public enum Method
        {
            #region Methods
            /// <summary>
            /// Linear method.
            /// </summary>
            Linear,
            /// <summary>
            /// Lagrange's method.
            /// </summary>
            Lagrange,
            /// <summary>
            /// Newton's method.
            /// </summary>
            Newton,
            /// <summary>
            /// Barycentric method.
            /// </summary>
            Barycentric,
            #endregion
        }
        #endregion
    }
}
