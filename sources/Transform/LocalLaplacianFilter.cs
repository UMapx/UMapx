using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the local Laplace pyramid filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://people.csail.mit.edu/sparis/publi/2011/siggraph/
    /// </remarks>
    /// </summary>
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
        /// <param name="radius">Radius</param>
        /// <param name="sigma">σ-parameter</param>
        /// <param name="n">Number of samples</param>
        /// <param name="levels">Number of levels</param>
        /// <param name="factor">Factor [-1, 1]</param>
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
        /// Gets or sets the number of samples.
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
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            llfilter(data, this.radius, this.sigma, this.factor, this.n, this.levels);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[] data)
        {
            llfilter(data, this.radius, this.sigma, this.factor, this.n, this.levels);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex32[,] data)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
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
        // Moscow, Russia.
        // **************************************************


        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="input">Input data</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="n">Number of steps</param>
        /// <param name="levels">Levels</param>
        /// <returns>Output data</returns>
        internal static void llfilter(float[,] input, int radius, float sigma, float factor, int n, int levels)
        {
            // exception
            if (factor == 0)
                return;

            // data
            int height = input.GetLength(0);
            int width = input.GetLength(1);
            int y, x, level, length = 256;
            float step = 1.0f / n;
            float min = 0.0f, max = 1.0f;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(Math.Min(height, width)) / Math.Log(2)), levels);
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels, radius);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels, radius);

            float[][,] input_gaussian_pyr = gpt.Forward(input);
            float[][,] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            float[][,] temp_laplace_pyr;
            float[,] I_temp, I_gaus, I_outp;
            float[] T;

            // do job
            for (float i = min; i <= max; i += step)
            {
                height = input.GetLength(0); width = input.GetLength(1);
                I_temp = new float[height, width];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++)
                    {
                        I_temp[y, x] = T[Maths.Byte(input[y, x] * (length - 1))];
                    }
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // pyramid reconstruction
                for (level = 0; level < n_levels; level++)
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
                            I_outp[y, x] += T[Maths.Byte(I_gaus[y, x] * (length - 1))] * I_temp[y, x];
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

            return;
        }
        /// <summary>
        /// Local laplacian filter.
        /// </summary>
        /// <param name="input">Input data</param>
        /// <param name="radius">Radius</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="n">Number of steps</param>
        /// <param name="levels">Levels</param>
        /// <returns>Output data</returns>
        internal static void llfilter(float[] input, int radius, float sigma, float factor, int n, int levels)
        {
            // exception
            if (factor == 0)
                return;

            // data
            int height = input.GetLength(0);
            int y, level, length = 256;
            float step = 1.0f / n;
            float min = 0.0f, max = 1.0f;

            // pyramids
            int n_levels = (int)Math.Min((Math.Log(height) / Math.Log(2)), levels);
            LaplacianPyramidTransform lpt = new LaplacianPyramidTransform(n_levels, radius);
            GaussianPyramidTransform gpt = new GaussianPyramidTransform(n_levels, radius);

            float[][] input_gaussian_pyr = gpt.Forward(input);
            float[][] output_laplace_pyr = lpt.Forward(input_gaussian_pyr);
            float[][] temp_laplace_pyr;
            float[] I_temp, I_gaus, I_outp;
            float[] T;

            // do job
            for (float i = min; i <= max; i += step)
            {
                height = input.GetLength(0);
                I_temp = new float[height];
                T = Rem(sigma, factor, i, length);

                // remapping function
                for (y = 0; y < height; y++)
                {
                    I_temp[y] = T[Maths.Byte(input[y] * (length - 1))];
                }

                temp_laplace_pyr = lpt.Forward(I_temp);
                T = Rec(i, step, length);

                // pyramid reconstruction
                for (level = 0; level < n_levels; level++)
                {
                    I_gaus = input_gaussian_pyr[level];
                    I_temp = temp_laplace_pyr[level];
                    I_outp = output_laplace_pyr[level];
                    height = I_outp.GetLength(0);

                    for (y = 0; y < height; y++)
                    {
                        I_outp[y] += T[Maths.Byte(I_gaus[y] * (length - 1))] * I_temp[y];
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

            return;
        }

        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="i">Increment</param>
        /// <param name="step">Step</param>
        /// <returns>Function</returns>
        internal static float Rec(float x, float i, float step)
        {
            float y = Math.Abs(x - i);
            return y < step ? (1.0f - y / step) : 0;
        }
        /// <summary>
        /// Reconstruct function.
        /// </summary>
        /// <param name="i">Increment</param>
        /// <param name="step">Step</param>
        /// <param name="length">Length of table</param>
        /// <returns>Table</returns>
        internal static float[] Rec(float i, float step, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rec(x / (float)length, i, step);
            }
            return table;
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="i">Increment</param>
        /// <returns>Function</returns>
        internal static float Rem(float x, float sigma, float factor, float i)
        {
            float z = 2 * sigma * sigma;
            float y = x - i;
            return factor * y * Maths.Exp(-y * y / z);
        }
        /// <summary>
        /// Remapping function.
        /// </summary>
        /// <param name="sigma">Sigma</param>
        /// <param name="factor">Factor</param>
        /// <param name="i">Increment</param>
        /// <param name="length">Length of table</param>
        /// <returns>Table</returns>
        internal static float[] Rem(float sigma, float factor, float i, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LocalLaplacianFilter.Rem(x / (float)length, sigma, factor, i);
            }
            return table;
        }
        #endregion
    }
}
