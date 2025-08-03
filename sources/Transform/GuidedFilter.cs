using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the guided filter.
    /// <remarks>
    /// This filter is a computationally effective analogue of a bilateral filter.
    /// More information can be found on the website:
    /// http://kaiminghe.com/eccv10/index.html
    /// </remarks>
    /// </summary>
    [Serializable]
    public class GuidedFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Epsilon.
        /// </summary>
        private float eps;
        /// <summary>
        /// Factor.
        /// </summary>
        private float factor;
        /// <summary>
        /// Radius.
        /// </summary>
        private int radius;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the guided filter.
        /// </summary>
        /// <param name="radius">Radius (>1)</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public GuidedFilter(int radius, float eps = 0.025f, float factor = -1.0f)
        {
            this.radius = radius;
            this.Eps = eps;
            this.factor = factor;
        }
        /// <summary>
        /// Gets or sets the radius.
        /// </summary>
        public int Radius
        {
            get
            {
                return this.radius;
            }
            set
            {
                this.radius = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of the epsilon (0, 1).
        /// <remarks>
        /// Optimal value ε = 0.025.
        /// </remarks>
        /// </summary>
        public float Eps
        {
            get
            {
                return this.eps;
            }
            set
            {
                this.eps = Maths.Float(value);
            }
        }
        /// <summary>
        /// Gets or sets the factor [-1, 1].
        /// </summary>
        public float Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Creates a guided filter with the specified parameters for a bilateral filter.
        /// </summary>
        /// <param name="s">σs</param>
        /// <param name="r">σr</param>
        /// <returns>Guided filter</returns>
        public static GuidedFilter FromBilateral(int s, float r)
        {
            return new GuidedFilter(s, r * r);
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                float[] copy = (float[])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0f + this.factor) * (data[i] - copy[i]) + copy[i];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                float[,] copy = (float[,])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0f + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex32[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                Complex32[] copy = (Complex32[])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex32[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                Complex32[,] copy = (Complex32[,])data.Clone();
                GuidedFilter.guidedfilter(copy, this.radius, this.eps);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
        }
        #endregion

        #region Private voids
        // **************************************************
        //                   GUIDED FILTER
        // **************************************************
        // ORIGINALS: Kaiming He, Jian Sun, and Xiaoou Tang.
        // Designed by Valery Asiryan (c), 2015-2018
        // **************************************************

        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(float[,] input, int r, float eps)
        {
            // Input signal properties:
            int l0 = input.GetLength(0), l1 = input.GetLength(1), i, j;

            // Calculating μ(I) and μ(I^2):
            float[,] x = (float[,])input.Clone();
            float[,] y = Matrice.Pow(input, 2.0f);

            // Applying fast box filter:
            x = x.Mean(r, r);
            y = y.Mean(r, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            float[,] c = new float[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    c[i, j] = y[i, j] - x[i, j] * x[i, j];

            // Calculating μ(a) and μ(b):
            float[,] a = new float[l0, l1];
            float[,] b = new float[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                {
                    a[i, j] = c[i, j] / (c[i, j] + eps);
                    b[i, j] = x[i, j] - a[i, j] * x[i, j];
                }

            // Applying fast box filter:
            a = a.Mean(r, r);
            b = b.Mean(r, r);

            // Calculating μ(a) * I + μ(b):
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    input[i, j] = a[i, j] * input[i, j] + b[i, j];

            return;
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(Complex32[,] input, int r, float eps)
        {
            // Input signal properties:
            int l0 = input.GetLength(0), l1 = input.GetLength(1), i, j;

            // Calculating μ(I) and μ(I^2):
            Complex32[,] x = (Complex32[,])input.Clone();
            Complex32[,] y = Matrice.Pow(input, 2.0f);

            // Applying fast box filter:
            x = x.Mean(r, r);
            y = y.Mean(r, r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            Complex32[,] c = new Complex32[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    c[i, j] = y[i, j] - x[i, j] * x[i, j];

            // Calculating μ(a) and μ(b):
            Complex32[,] a = new Complex32[l0, l1];
            Complex32[,] b = new Complex32[l0, l1];
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                {
                    a[i, j] = c[i, j] / (c[i, j] + eps);
                    b[i, j] = x[i, j] - a[i, j] * x[i, j];
                }

            // Applying fast box filter:
            a = a.Mean(r, r);
            b = b.Mean(r, r);

            // Calculating μ(a) * I + μ(b):
            for (i = 0; i < l0; i++)
                for (j = 0; j < l1; j++)
                    input[i, j] = a[i, j] * input[i, j] + b[i, j];
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(float[] input, int r, float eps)
        {
            // Input signal properties:
            int length = input.Length, i;

            // Calculating μ(I) and μ(I^2):
            float[] x = (float[])input.Clone();
            float[] y = Matrice.Pow(input, 2.0f);

            // Applying fast box filter:
            x = x.Mean(r);
            y = y.Mean(r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            float[] c = new float[length];
            for (i = 0; i < length; i++)
                c[i] = y[i] - x[i] * x[i];

            // Calculating μ(a) and μ(b):
            float[] a = new float[length];
            float[] b = new float[length];
            for (i = 0; i < length; i++)
            {
                a[i] = c[i] / (c[i] + eps);
                b[i] = x[i] - a[i] * x[i];
            }

            // Applying fast box filter:
            a = a.Mean(r);
            b = b.Mean(r);

            // Calculating μ(a) * I + μ(b):
            for (i = 0; i < length; i++)
                input[i] = a[i] * input[i] + b[i];
        }
        /// <summary>
        /// Guided filer function.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Filter size</param>
        /// <param name="eps">Epsilon (0, 1)</param>
        internal static void guidedfilter(Complex32[] input, int r, float eps)
        {
            // Input signal properties:
            int length = input.Length, i;

            // Calculating μ(I) and μ(I^2):
            Complex32[] x = (Complex32[])input.Clone();
            Complex32[] y = Matrice.Pow(input, 2.0f);

            // Applying fast box filter:
            x = x.Mean(r);
            y = y.Mean(r);

            // Calculating cov(I):
            // This is the covariance of input in each local patch:
            Complex32[] c = new Complex32[length];
            for (i = 0; i < length; i++)
                c[i] = y[i] - x[i] * x[i];

            // Calculating μ(a) and μ(b):
            Complex32[] a = new Complex32[length];
            Complex32[] b = new Complex32[length];
            for (i = 0; i < length; i++)
            {
                a[i] = c[i] / (c[i] + eps);
                b[i] = x[i] - a[i] * x[i];
            }

            // Applying fast box filter:
            a = a.Mean(r);
            b = b.Mean(r);

            // Calculating μ(a) * I + μ(b):
            for (i = 0; i < length; i++)
                input[i] = a[i] * input[i] + b[i];

            return;
        }
        #endregion
    }
}
