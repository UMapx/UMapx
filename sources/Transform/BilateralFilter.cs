using System;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the bilateral filter.
    /// <remarks>
    /// Optimized implementation.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Bilateral_filter
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BilateralFilter : IFilter
    {
        #region Private data
        private int radius;
        private float sigma;
        private int levels;
        private float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bilateral filter.
        /// </summary>
        /// <param name="radius">Radius (>1)</param>
        /// <param name="sigma">Range Gaussian sigma (>0)</param>
        /// <param name="levels">Number of quantization levels for intensity</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public BilateralFilter(int radius, float sigma = 0.1f, int levels = 32, float factor = -1.0f)
        {
            this.Radius = radius;
            this.Sigma = sigma;
            this.Levels = levels;
            this.Factor = factor;
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
        /// Gets or sets the value of the range Gaussian sigma (>0).
        /// </summary>
        public float Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                this.sigma = value;
            }
        }
        /// <summary>
        /// Gets or sets the number of quantization levels for intensity.
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
                BilateralFilter.Bilateralfilter(copy, this.radius, this.sigma, this.levels);

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
                BilateralFilter.Bilateralfilter(copy, this.radius, this.sigma, this.levels);

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
                BilateralFilter.Bilateralfilter(copy, this.radius, this.sigma, this.levels);

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
                BilateralFilter.Bilateralfilter(copy, this.radius, this.sigma, this.levels);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Applies fast bilateral filter.
        /// </summary>
        /// <param name="array">Input signal</param>
        /// <param name="r">Radius</param>
        /// <param name="s">Range Gaussian sigma</param>
        /// <param name="samples">Number of quantization levels for intensity</param>
        /// <returns>Filtered image</returns>
        private static void Bilateralfilter(float[,] array, int r, float s = 0.1f, int samples = 32)
        {
            int height = array.GetLength(0);
            int width = array.GetLength(1);

            float[] levels = new float[samples];

            for (int i = 0; i < samples; i++)
            {
                levels[i] = i / (float)(samples - 1);
            }

            float[,] output = new float[height, width];
            float[,] norm = new float[height, width];
            float[][,] outputLocals = new float[samples][,];
            float[][,] normLocals = new float[samples][,];

            Parallel.For(0, samples, l =>
            {
                float[,] rangeWeight = new float[height, width];
                float level = levels[l];

                float rangeSigma2 = 2 * s * s;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        float diff = array[y, x] - level;
                        rangeWeight[y, x] = Maths.Exp(-diff * diff / rangeSigma2);
                    }
                }

                float[,] blurred = Matrix.Mean(array, rangeWeight, r, r);

                float[,] outputLocal = new float[height, width];
                float[,] normLocal = new float[height, width];

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        float w = rangeWeight[y, x];
                        outputLocal[y, x] = blurred[y, x] * w;
                        normLocal[y, x] = w;
                    }
                }

                outputLocals[l] = outputLocal;
                normLocals[l] = normLocal;
            });

            for (int l = 0; l < samples; l++)
            {
                var o = outputLocals[l];
                var n = normLocals[l];

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        output[y, x] += o[y, x];
                        norm[y, x] += n[y, x];
                    }
                }
            }

            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    array[y, x] = output[y, x] / (norm[y, x] + 1e-8f);
                }
            });
        }
        /// <summary>
        /// Applies fast bilateral filter.
        /// </summary>
        /// <param name="array">Input signal</param>
        /// <param name="r">Radius</param>
        /// <param name="s">Range Gaussian sigma</param>
        /// <param name="samples">Number of quantization levels for intensity</param>
        /// <returns>Filtered image</returns>
        private static void Bilateralfilter(Complex32[,] array, int r, float s = 0.1f, int samples = 32)
        {
            int height = array.GetLength(0);
            int width = array.GetLength(1);

            Complex32[] levels = new Complex32[samples];

            for (int i = 0; i < samples; i++)
            {
                levels[i] = i / (float)(samples - 1);
            }

            Complex32[,] output = new Complex32[height, width];
            Complex32[,] norm = new Complex32[height, width];
            Complex32[][,] outputLocals = new Complex32[samples][,];
            Complex32[][,] normLocals = new Complex32[samples][,];

            Parallel.For(0, samples, l =>
            {
                Complex32[,] rangeWeight = new Complex32[height, width];
                Complex32 level = levels[l];

                float rangeSigma2 = 2 * s * s;

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Complex32 diff = array[y, x] - level;
                        rangeWeight[y, x] = Maths.Exp(-diff * diff / rangeSigma2);
                    }
                }

                Complex32[,] blurred = Matrix.Mean(array, rangeWeight, r, r);

                Complex32[,] outputLocal = new Complex32[height, width];
                Complex32[,] normLocal = new Complex32[height, width];

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Complex32 w = rangeWeight[y, x];
                        outputLocal[y, x] = blurred[y, x] * w;
                        normLocal[y, x] = w;
                    }
                }

                outputLocals[l] = outputLocal;
                normLocals[l] = normLocal;
            });

            for (int l = 0; l < samples; l++)
            {
                var o = outputLocals[l];
                var n = normLocals[l];

                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        output[y, x] += o[y, x];
                        norm[y, x] += n[y, x];
                    }
                }
            }

            Parallel.For(0, height, y =>
            {
                for (int x = 0; x < width; x++)
                {
                    array[y, x] = output[y, x] / (norm[y, x] + 1e-8f);
                }
            });
        }
        /// <summary>
        /// Applies fast bilateral filter.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Radius</param>
        /// <param name="s">Range Gaussian sigma</param>
        /// <param name="samples">Number of quantization levels for intensity</param>
        /// <returns>Filtered image</returns>
        private static void Bilateralfilter(float[] input, int r, float s = 0.1f, int samples = 32)
        {
            int length = input.Length;
            float[] levels = new float[samples];

            for (int i = 0; i < samples; i++)
            {
                levels[i] = i / (float)(samples - 1);
            }

            float[] output = new float[length];
            float[] norm = new float[length];
            float[][] outputLocals = new float[samples][];
            float[][] normLocals = new float[samples][];

            float rangeSigma2 = 2 * s * s;

            Parallel.For(0, samples, l =>
            {
                float level = levels[l];
                float[] rangeWeight = new float[length];

                for (int i = 0; i < length; i++)
                {
                    float diff = input[i] - level;
                    rangeWeight[i] = Maths.Exp(-diff * diff / rangeSigma2);
                }

                float[] blurred = Matrix.Mean(input, rangeWeight, r);

                float[] outputLocal = new float[length];
                float[] normLocal = new float[length];

                for (int i = 0; i < length; i++)
                {
                    float w = rangeWeight[i];
                    outputLocal[i] = blurred[i] * w;
                    normLocal[i] = w;
                }

                outputLocals[l] = outputLocal;
                normLocals[l] = normLocal;
            });

            for (int l = 0; l < samples; l++)
            {
                var o = outputLocals[l];
                var n = normLocals[l];

                for (int i = 0; i < length; i++)
                {
                    output[i] += o[i];
                    norm[i] += n[i];
                }
            }

            for (int i = 0; i < length; i++)
            {
                input[i] = output[i] / (norm[i] + 1e-8f);
            }
        }
        /// <summary>
        /// Applies fast bilateral filter.
        /// </summary>
        /// <param name="input">Input signal</param>
        /// <param name="r">Radius</param>
        /// <param name="s">Range Gaussian sigma</param>
        /// <param name="samples">Number of quantization levels for intensity</param>
        /// <returns>Filtered image</returns>
        private static void Bilateralfilter(Complex32[] input, int r, float s = 0.1f, int samples = 32)
        {
            int length = input.Length;
            Complex32[] levels = new Complex32[samples];

            for (int i = 0; i < samples; i++)
            {
                levels[i] = i / (float)(samples - 1);
            }

            Complex32[] output = new Complex32[length];
            Complex32[] norm = new Complex32[length];
            Complex32[][] outputLocals = new Complex32[samples][];
            Complex32[][] normLocals = new Complex32[samples][];

            float rangeSigma2 = 2 * s * s;

            Parallel.For(0, samples, l =>
            {
                Complex32 level = levels[l];
                Complex32[] rangeWeight = new Complex32[length];

                for (int i = 0; i < length; i++)
                {
                    Complex32 diff = input[i] - level;
                    rangeWeight[i] = Maths.Exp(-diff * diff / rangeSigma2);
                }

                Complex32[] blurred = Matrix.Mean(input, rangeWeight, r);

                Complex32[] outputLocal = new Complex32[length];
                Complex32[] normLocal = new Complex32[length];

                for (int i = 0; i < length; i++)
                {
                    Complex32 w = rangeWeight[i];
                    outputLocal[i] = blurred[i] * w;
                    normLocal[i] = w;
                }

                outputLocals[l] = outputLocal;
                normLocals[l] = normLocal;
            });

            for (int l = 0; l < samples; l++)
            {
                var o = outputLocals[l];
                var n = normLocals[l];

                for (int i = 0; i < length; i++)
                {
                    output[i] += o[i];
                    norm[i] += n[i];
                }
            }

            for (int i = 0; i < length; i++)
            {
                input[i] = output[i] / (norm[i] + 1e-8f);
            }
        }
        #endregion
    }
}
