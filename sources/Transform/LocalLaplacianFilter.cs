using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the local Laplace pyramid filter.
    /// </summary>
    /// <remarks>
    /// Intensity remapping is sampled on [0, 1] and interpolated from 256-entry
    /// tables. Only detail levels are modified; the coarsest Gaussian level is
    /// preserved. Zero width, zero sampling intervals, or no detail levels leave
    /// the input unchanged. Complex filtering is not supported.
    /// More information can be found on the website:
    /// <see href="https://people.csail.mit.edu/sparis/publi/2011/siggraph/"/>.
    /// </remarks>
    [Serializable]
    public class LocalLaplacianFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Sigma.
        /// </summary>
        private float sigma;
        /// <summary>
        /// Factor.
        /// </summary>
        private float factor;
        /// <summary>
        /// Number of samples.
        /// </summary>
        private int n;
        /// <summary>
        /// Number of levels.
        /// </summary>
        private int levels;
        /// <summary>
        /// Radius.
        /// </summary>
        private int radius;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local Laplace pyramid filter.
        /// </summary>
        /// <param name="radius">Radius.</param>
        /// <param name="sigma">σ-parameter.</param>
        /// <param name="n">Number of intensity sampling intervals; zero disables filtering.</param>
        /// <param name="levels">Number of levels.</param>
        /// <param name="factor">Factor [-1, 1].</param>
        public LocalLaplacianFilter(int radius = 2, float sigma = 0.05f, int n = 10, int levels = 10, float factor = -1.0f)
        {
            this.Radius = radius;
            this.Sigma = sigma;
            this.N = n;
            this.Levels = levels;
            this.Factor = factor;
        }
        /// <summary>
        /// Gets or sets radius.
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
        /// Gets or sets the value of σ-parameter.
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                this.sigma = Maths.Float(value);
            }
        }
        /// <summary>
        /// Gets or sets the factor.
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
        /// <summary>
        /// Gets or sets the number of sampling intervals; zero disables the filter.
        /// </summary>
        public int N
        {
            get
            {
                return this.n;
            }
            set
            {
                this.n = Math.Max(value, 0);
            }
        }
        /// <summary>
        /// Gets or sets the number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                this.levels = value;
            }
        }
        #endregion

        #region Apply voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix.</param>
        public void Apply(float[,] data)
        {
            Llfilter(data, this.radius, this.sigma, this.factor, this.n, this.levels);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix.</param>
        public void Apply(float[] data)
        {
            Llfilter(data, this.radius, this.sigma, this.factor, this.n, this.levels);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix.</param>
        public void Apply(Complex32[,] data)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix.</param>
        public void Apply(Complex32[] data)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        // **************************************************
        //            Local Laplacian Filter
        // **************************************************
        // This function implements edge-aware detail and 
        // tone manipulation as described in:
        // "Fast and Robust Pyramid-based Image Processing"
        // Mathieu Aubry, Sylvain Paris, Samuel W. Hasinoff, 
        // Jan Kautz, and Fredo Durand.
        // MIT technical report, November 2011.
        // 
        // Designed by Valery Asiryan (c), 2015-2020
        // **************************************************


        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="input">Input data.</param>
        /// <param name="radius">Radius.</param>
        /// <param name="sigma">Sigma.</param>
        /// <param name="factor">Factor.</param>
        /// <param name="n">Number of steps.</param>
        /// <param name="levels">Levels.</param>
        private static void Llfilter(float[,] input, int radius, float sigma, float factor, int n, int levels)
        {
            // exception
            if (factor == 0 || sigma == 0 || n == 0 || levels <= 1)
                return;

            // data
            int height = input.GetLength(0);
            int width = input.GetLength(1);
            int y, x, level, length = 256;
            float step = 1.0f / n;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(Math.Min(height, width)) / Math.Log(2)), levels);
            if (n_levels < 2) return;
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels, radius);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels, radius);

            float[][,] input_gaussian_pyr = gpt.Forward(input);
            float[][,] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            float[][,] temp_laplace_pyr;
            float[,] I_temp, I_gaus, I_outp;
            float[] T;

            // do job
            // Include both endpoints without accumulated floating-point step error.
            for (int sample = 0; sample <= n; sample++)
            {
                float i = sample * step;
                height = input.GetLength(0); width = input.GetLength(1);
                I_temp = new float[height, width];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++)
                    {
                        I_temp[y, x] = SampleTable(T, input[y, x]);
                    }
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // Modify detail bands only; preserve the coarsest Gaussian level.
                for (level = 0; level < n_levels - 1; level++)
                {
                    I_gaus = input_gaussian_pyr[level];
                    I_temp = temp_laplace_pyr[level];
                    I_outp = output_laplace_pyr[level];
                    height = I_outp.GetLength(0);
                    width = I_outp.GetLength(1);

                    for (y = 0; y < height; y++)
                    {
                        for (x = 0; x < width; x++)
                        {
                            I_outp[y, x] += SampleTable(T, I_gaus[y, x]) * I_temp[y, x];
                        }
                    }

                    output_laplace_pyr[level] = I_outp;
                }
            }

            // backward transform
            I_outp = lpt.Backward(output_laplace_pyr);
            height = input.GetLength(0);
            width = input.GetLength(1);

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++)
                {
                    input[y, x] = I_outp[y, x];
                }
            }
        }
        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="input">Input data.</param>
        /// <param name="radius">Radius.</param>
        /// <param name="sigma">Sigma.</param>
        /// <param name="factor">Factor.</param>
        /// <param name="n">Number of steps.</param>
        /// <param name="levels">Levels.</param>
        private static void Llfilter(float[] input, int radius, float sigma, float factor, int n, int levels)
        {
            // exception
            if (factor == 0 || sigma == 0 || n == 0 || levels <= 1)
                return;

            // data
            int height = input.GetLength(0);
            int y, level, length = 256;
            float step = 1.0f / n;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(height) / Math.Log(2)), levels);
            if (n_levels < 2) return;
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels, radius);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels, radius);

            float[][] input_gaussian_pyr = gpt.Forward(input);
            float[][] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            float[][] temp_laplace_pyr;
            float[] I_temp, I_gaus, I_outp;
            float[] T;

            // do job
            // Include both endpoints without accumulated floating-point step error.
            for (int sample = 0; sample <= n; sample++)
            {
                float i = sample * step;
                height = input.GetLength(0);
                I_temp = new float[height];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    I_temp[y] = SampleTable(T, input[y]);
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // Modify detail bands only; preserve the coarsest Gaussian level.
                for (level = 0; level < n_levels - 1; level++)
                {
                    I_gaus = input_gaussian_pyr[level];
                    I_temp = temp_laplace_pyr[level];
                    I_outp = output_laplace_pyr[level];
                    height = I_outp.GetLength(0);

                    for (y = 0; y < height; y++)
                    {
                        I_outp[y] += SampleTable(T, I_gaus[y]) * I_temp[y];
                    }

                    output_laplace_pyr[level] = I_outp;
                }
            }

            // backward transform
            I_outp = lpt.Backward(output_laplace_pyr);
            height = input.GetLength(0);

            for (y = 0; y < height; y++)
            {
                input[y] = I_outp[y];
            }
        }

        /// <summary>
        /// Interpolates a lookup table sampled uniformly over the unit interval.
        /// </summary>
        /// <param name="table">At least two samples, including both interval endpoints.</param>
        /// <param name="value">Intensity, clamped to [0, 1] before lookup.</param>
        /// <returns>The linear interpolation of adjacent table samples.</returns>
        private static float SampleTable(float[] table, float value)
        {
            float position = Maths.Float(value) * (table.Length - 1);
            int index = Math.Min((int)position, table.Length - 2);
            float fraction = position - index;
            return table[index] + fraction * (table[index + 1] - table[index]);
        }

        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="i">Increment.</param>
        /// <param name="step">Step.</param>
        /// <returns>Function.</returns>
        private static float Rec(float x, float i, float step)
        {
            float y = Math.Abs(x - i);
            return y < step ? (1.0f - y / step) : 0;
        }
        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="i">Increment.</param>
        /// <param name="step">Step.</param>
        /// <param name="length">Length of table.</param>
        /// <returns>Table.</returns>
        private static float[] Rec(float i, float step, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rec(x / (float)(length - 1), i, step);
            }
            return table;
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="x">Value.</param>
        /// <param name="sigma">Sigma.</param>
        /// <param name="factor">Factor.</param>
        /// <param name="i">Increment.</param>
        /// <returns>Function.</returns>
        private static float Rem(float x, float sigma, float factor, float i)
        {
            float z = 2 * sigma * sigma;
            float y = x - i;
            return factor * y * Maths.Exp(-y * y / z);
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="sigma">Sigma.</param>
        /// <param name="factor">Factor.</param>
        /// <param name="i">Increment.</param>
        /// <param name="length">Length of table.</param>
        /// <returns>Table.</returns>
        private static float[] Rem(float sigma, float factor, float i, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rem(x / (float)(length - 1), sigma, factor, i);
            }
            return table;
        }
        #endregion
    }
}
