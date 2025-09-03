using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the bilateral filter.
    /// <remarks>
    /// Fast implementation of the bilateral filter (real-time).
    /// More information can be found on the website:
    /// https://www.researchgate.net/publication/220184523_Real-time_edge-aware_image_processing_with_the_bilateral_grid
    /// </remarks>
    /// </summary>
    [Serializable]
    public class BilateralGridFilter : IFilter
    {
        #region Private data
        private float sigmaSpatial;
        private float sigmaRange;
        private float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bilateral grid filter.
        /// </summary>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public BilateralGridFilter(float sigmaSpatial = 4f, float sigmaRange = 0.1f, float factor = -1.0f)
        {
            this.SigmaSpatial = sigmaSpatial;
            this.SigmaRange = sigmaRange;
            this.Factor = factor;
        }

        /// <summary>
        /// Gets or sets the spatial smoothing factor.
        /// </summary>
        public float SigmaSpatial
        {
            get
            {
                return this.sigmaSpatial;
            }
            set
            {
                this.sigmaSpatial = value;
            }
        }

        /// <summary>
        /// Gets or sets the range smoothing factor.
        /// </summary>
        public float SigmaRange
        {
            get
            {
                return this.sigmaRange;
            }
            set
            {
                this.sigmaRange = value;
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
                BilateralGridFilter.Bilateralgridfilter(copy, this.sigmaSpatial, this.sigmaRange);

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
                BilateralGridFilter.Bilateralgridfilter(copy, this.sigmaSpatial, this.sigmaRange);

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
                BilateralGridFilter.Bilateralgridfilter(copy, this.sigmaSpatial, this.sigmaRange);

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
                BilateralGridFilter.Bilateralgridfilter(copy, this.sigmaSpatial, this.sigmaRange);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Applies bilateral grid filter.
        /// </summary>
        /// <param name="input">Input</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        internal static void Bilateralgridfilter(float[] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            int length = input.Length;

            int gridZ = (int)(1.0f / sigmaRange) + 1;
            int gridX = (int)Math.Ceiling(length / sigmaSpatial) + 2;

            float[,] gridData = new float[gridZ, gridX];
            float[,] gridWeight = new float[gridZ, gridX];

            // SPLAT step
            for (int x = 0; x < length; x++)
            {
                float value = input[x];
                int z = (int)(value / sigmaRange);
                int gx = (int)(x / sigmaSpatial);

                gridData[z, gx] += value;
                gridWeight[z, gx] += 1f;
            }

            // BLUR step
            BilateralGridFilter.Boxblur2d(gridData, 1);
            BilateralGridFilter.Boxblur2d(gridWeight, 1);

            // SLICE step
            for (int x = 0; x < length; x++)
            {
                float value = input[x];
                float z = value / sigmaRange;
                float gx = x / sigmaSpatial;

                int z0 = (int)Math.Floor(z);
                int x0 = (int)Math.Floor(gx);

                float dz = z - z0;
                float dx = gx - x0;

                float v = 0f, w = 0f;

                for (int zz = 0; zz <= 1; zz++)
                {
                    for (int xx = 0; xx <= 1; xx++)
                    {
                        int zi = MathF.Range(z0 + zz, 0, gridZ - 1);
                        int xi = MathF.Range(x0 + xx, 0, gridX - 1);

                        float weight = (1 - Math.Abs(dz - zz)) * (1 - Math.Abs(dx - xx));

                        v += gridData[zi, xi] * weight;
                        w += gridWeight[zi, xi] * weight;
                    }
                }

                input[x] = (w > 1e-5f) ? v / w : value;
            }
        }
        /// <summary>
        /// Applies bilateral grid filter.
        /// </summary>
        /// <param name="input">Input</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        internal static void Bilateralgridfilter(float[,] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            int height = input.GetLength(0);
            int width = input.GetLength(1);

            int gridZ = (int)(1.0f / sigmaRange) + 1;
            int gridY = (int)Math.Ceiling(height / sigmaSpatial) + 2;
            int gridX = (int)Math.Ceiling(width / sigmaSpatial) + 2;

            float[,,] gridData = new float[gridZ, gridY, gridX];
            float[,,] gridWeight = new float[gridZ, gridY, gridX];

            // SPLAT step: accumulate data into the grid
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float value = input[y, x];
                    int z = (int)(value / sigmaRange);
                    int gy = (int)(y / sigmaSpatial);
                    int gx = (int)(x / sigmaSpatial);

                    gridData[z, gy, gx] += value;
                    gridWeight[z, gy, gx] += 1f;
                }
            }

            // BLUR step: apply 3D Gaussian blur (simplified using box filter)
            BilateralGridFilter.Boxblur3d(gridData, 1);
            BilateralGridFilter.Boxblur3d(gridWeight, 1);

            // SLICE step: reconstruct output image
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float value = input[y, x];
                    float z = value / sigmaRange;
                    float gy = y / sigmaSpatial;
                    float gx = x / sigmaSpatial;

                    // Trilinear interpolation
                    int z0 = (int)Math.Floor(z);
                    int y0 = (int)Math.Floor(gy);
                    int x0 = (int)Math.Floor(gx);

                    float dz = z - z0;
                    float dy = gy - y0;
                    float dx = gx - x0;

                    float v = 0f, w = 0f;

                    for (int zz = 0; zz <= 1; zz++)
                    {
                        for (int yy = 0; yy <= 1; yy++)
                        {
                            for (int xx = 0; xx <= 1; xx++)
                            {
                                int zi = MathF.Range(z0 + zz, 0, gridZ - 1);
                                int yi = MathF.Range(y0 + yy, 0, gridY - 1);
                                int xi = MathF.Range(x0 + xx, 0, gridX - 1);

                                float weight = (1 - Math.Abs(dz - zz)) *
                                               (1 - Math.Abs(dy - yy)) *
                                               (1 - Math.Abs(dx - xx));

                                v += gridData[zi, yi, xi] * weight;
                                w += gridWeight[zi, yi, xi] * weight;
                            }
                        }
                    }

                    input[y, x] = (w > 1e-5f) ? v / w : value;
                }
            }
        }
        /// <summary>
        /// Applies bilateral grid filter.
        /// </summary>
        /// <param name="input">Input</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        internal static void Bilateralgridfilter(Complex32[] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            int length = input.Length;

            int gridZ = (int)(1.0f / sigmaRange) + 1;
            int gridX = (int)Math.Ceiling(length / sigmaSpatial) + 2;

            Complex32[,] gridData = new Complex32[gridZ, gridX];
            Complex32[,] gridWeight = new Complex32[gridZ, gridX];

            // SPLAT step
            for (int x = 0; x < length; x++)
            {
                float value = input[x].Abs;
                int z = (int)(value / sigmaRange);
                int gx = (int)(x / sigmaSpatial);

                gridData[z, gx] += value;
                gridWeight[z, gx] += 1f;
            }

            // BLUR step
            BilateralGridFilter.Boxblur2d(gridData, 1);
            BilateralGridFilter.Boxblur2d(gridWeight, 1);

            // SLICE step
            for (int x = 0; x < length; x++)
            {
                float value = input[x].Abs;
                float z = value / sigmaRange;
                float gx = x / sigmaSpatial;

                int z0 = (int)Math.Floor(z);
                int x0 = (int)Math.Floor(gx);

                float dz = z - z0;
                float dx = gx - x0;

                Complex32 v = 0f, w = 0f;

                for (int zz = 0; zz <= 1; zz++)
                {
                    for (int xx = 0; xx <= 1; xx++)
                    {
                        int zi = MathF.Range(z0 + zz, 0, gridZ - 1);
                        int xi = MathF.Range(x0 + xx, 0, gridX - 1);

                        float weight = (1 - Math.Abs(dz - zz)) * (1 - Math.Abs(dx - xx));

                        v += gridData[zi, xi] * weight;
                        w += gridWeight[zi, xi] * weight;
                    }
                }

                input[x] = (w.Abs > 1e-5f) ? v / w : value;
            }
        }
        /// <summary>
        /// Applies bilateral grid filter.
        /// </summary>
        /// <param name="input">Input</param>
        /// <param name="sigmaSpatial">Spatial smoothing factor</param>
        /// <param name="sigmaRange">Range smoothing factor</param>
        internal static void Bilateralgridfilter(Complex32[,] input, float sigmaSpatial = 4f, float sigmaRange = 0.1f)
        {
            int height = input.GetLength(0);
            int width = input.GetLength(1);

            int gridZ = (int)(1.0f / sigmaRange) + 1;
            int gridY = (int)Math.Ceiling(height / sigmaSpatial) + 2;
            int gridX = (int)Math.Ceiling(width / sigmaSpatial) + 2;

            Complex32[,,] gridData = new Complex32[gridZ, gridY, gridX];
            Complex32[,,] gridWeight = new Complex32[gridZ, gridY, gridX];

            // SPLAT step: accumulate data into the grid
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float value = input[y, x].Abs;
                    int z = (int)(value / sigmaRange);
                    int gy = (int)(y / sigmaSpatial);
                    int gx = (int)(x / sigmaSpatial);

                    gridData[z, gy, gx] += value;
                    gridWeight[z, gy, gx] += 1f;
                }
            }

            // BLUR step: apply 3D Gaussian blur (simplified using box filter)
            BilateralGridFilter.Boxblur3d(gridData, 1);
            BilateralGridFilter.Boxblur3d(gridWeight, 1);

            // SLICE step: reconstruct output image
            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < width; x++)
                {
                    float value = input[y, x].Abs;
                    float z = value / sigmaRange;
                    float gy = y / sigmaSpatial;
                    float gx = x / sigmaSpatial;

                    // Trilinear interpolation
                    int z0 = (int)Math.Floor(z);
                    int y0 = (int)Math.Floor(gy);
                    int x0 = (int)Math.Floor(gx);

                    float dz = z - z0;
                    float dy = gy - y0;
                    float dx = gx - x0;

                    Complex32 v = 0f, w = 0f;

                    for (int zz = 0; zz <= 1; zz++)
                    {
                        for (int yy = 0; yy <= 1; yy++)
                        {
                            for (int xx = 0; xx <= 1; xx++)
                            {
                                int zi = MathF.Range(z0 + zz, 0, gridZ - 1);
                                int yi = MathF.Range(y0 + yy, 0, gridY - 1);
                                int xi = MathF.Range(x0 + xx, 0, gridX - 1);

                                float weight = (1 - Math.Abs(dz - zz)) *
                                               (1 - Math.Abs(dy - yy)) *
                                               (1 - Math.Abs(dx - xx));

                                v += gridData[zi, yi, xi] * weight;
                                w += gridWeight[zi, yi, xi] * weight;
                            }
                        }
                    }

                    input[y, x] = (w.Abs > 1e-5f) ? v / w : value;
                }
            }
        }

        /// <summary>
        /// Implements 2D box blur filter.
        /// </summary>
        /// <param name="grid">Grid</param>
        /// <param name="radius">Radius</param>
        private static void Boxblur2d(float[,] grid, int radius)
        {
            int depth = grid.GetLength(0);
            int width = grid.GetLength(1);

            float[,] temp = (float[,])grid.Clone();

            for (int z = 0; z < depth; z++)
            {
                for (int x = 0; x < width; x++)
                {
                    float sum = 0f;
                    int count = 0;

                    for (int dz = -radius; dz <= radius; dz++)
                    {
                        for (int dx = -radius; dx <= radius; dx++)
                        {
                            int zi = MathF.Range(z + dz, 0, depth - 1);
                            int xi = MathF.Range(x + dx, 0, width - 1);

                            sum += temp[zi, xi];
                            count++;
                        }
                    }

                    grid[z, x] = sum / count;
                }
            }
        }
        /// <summary>
        /// Implements 3D box blur filter.
        /// </summary>
        /// <param name="grid">Grid</param>
        /// <param name="radius">Radius</param>
        private static void Boxblur3d(float[,,] grid, int radius)
        {
            int depth = grid.GetLength(0);
            int height = grid.GetLength(1);
            int width = grid.GetLength(2);

            float[,,] temp = (float[,,])grid.Clone();

            for (int z = 0; z < depth; z++)
            {
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        float sum = 0f;
                        int count = 0;

                        for (int dz = -radius; dz <= radius; dz++)
                        {
                            for (int dy = -radius; dy <= radius; dy++)
                            {
                                for (int dx = -radius; dx <= radius; dx++)
                                {
                                    int zi = MathF.Range(z + dz, 0, depth - 1);
                                    int yi = MathF.Range(y + dy, 0, height - 1);
                                    int xi = MathF.Range(x + dx, 0, width - 1);

                                    sum += temp[zi, yi, xi];
                                    count++;
                                }
                            }
                        }

                        grid[z, y, x] = sum / count;
                    }
                }
            }
        }
        /// <summary>
        /// Implements 2D box blur filter.
        /// </summary>
        /// <param name="grid">Grid</param>
        /// <param name="radius">Radius</param>
        private static void Boxblur2d(Complex32[,] grid, int radius)
        {
            int depth = grid.GetLength(0);
            int width = grid.GetLength(1);

            Complex32[,] temp = (Complex32[,])grid.Clone();

            for (int z = 0; z < depth; z++)
            {
                for (int x = 0; x < width; x++)
                {
                    Complex32 sum = 0f;
                    int count = 0;

                    for (int dz = -radius; dz <= radius; dz++)
                    {
                        for (int dx = -radius; dx <= radius; dx++)
                        {
                            int zi = MathF.Range(z + dz, 0, depth - 1);
                            int xi = MathF.Range(x + dx, 0, width - 1);

                            sum += temp[zi, xi];
                            count++;
                        }
                    }

                    grid[z, x] = sum / count;
                }
            }
        }
        /// <summary>
        /// Implements 3D box blur filter.
        /// </summary>
        /// <param name="grid">Grid</param>
        /// <param name="radius">Radius</param>
        private static void Boxblur3d(Complex32[,,] grid, int radius)
        {
            int depth = grid.GetLength(0);
            int height = grid.GetLength(1);
            int width = grid.GetLength(2);

            Complex32[,,] temp = (Complex32[,,])grid.Clone();

            for (int z = 0; z < depth; z++)
            {
                for (int y = 0; y < height; y++)
                {
                    for (int x = 0; x < width; x++)
                    {
                        Complex32 sum = 0f;
                        int count = 0;

                        for (int dz = -radius; dz <= radius; dz++)
                        {
                            for (int dy = -radius; dy <= radius; dy++)
                            {
                                for (int dx = -radius; dx <= radius; dx++)
                                {
                                    int zi = MathF.Range(z + dz, 0, depth - 1);
                                    int yi = MathF.Range(y + dy, 0, height - 1);
                                    int xi = MathF.Range(x + dx, 0, width - 1);

                                    sum += temp[zi, yi, xi];
                                    count++;
                                }
                            }
                        }

                        grid[z, y, x] = sum / count;
                    }
                }
            }
        }
        #endregion

    }

}
