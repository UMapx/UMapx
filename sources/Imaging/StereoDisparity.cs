using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the stereo disparity filter for a pair of images.
    /// <remarks>
    /// More information can be found on the website:
    /// https://docs.opencv.org/master/d3/d14/tutorial_ximgproc_disparity_filtering.html
    /// </remarks>
    /// </summary>
    [Serializable]
    public class StereoDisparity
    {
        #region Class components
        /// <summary>
        /// Initializes the stereo disparity filter for a pair of images.
        /// </summary>
        /// <param name="disparity">Disparity</param>
        /// <param name="window">Window size</param>
        /// <param name="weight">Gradient weight</param>
        public StereoDisparity(int disparity = 50, int window = 8, float weight = 5)
        {
            Disparity = disparity;
            Window = window;
            Weight = weight;
        }
        /// <summary>
        /// Gets or sets disparity.
        /// </summary>
        public int Disparity { get; set; }
        /// <summary>
        /// Gets or sets window size.
        /// </summary>
        public int Window { get; set; }
        /// <summary>
        /// Gets or sets weight.
        /// </summary>
        public float Weight { get; set; }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public float[,] Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // images to matrices
            var left = BitmapMatrix.ToGrayscale(bmData);
            var right = BitmapMatrix.ToGrayscale(bmSrc);

            // apply filter
            var output = DisparityEstimator(left, right, Window, Disparity, Weight);

            // return result
            return output;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public float[,] Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            var output = Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);

            return output;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Estimates a disparity map from stereo image pairs.
        /// </summary>
        /// <param name="left">Left image matrix</param>
        /// <param name="right">Right image matrix</param>
        /// <param name="win">Window size</param>
        /// <param name="max_dis">Maximum disparity value</param>
        /// <param name="weight">Gradient weight</param>
        /// <param name="apply_median">Apply median filtering or not</param>
        /// <returns>Disparity matrix</returns>
        private float[,] DisparityEstimator(float[,] left, float[,] right, int win, int max_dis, float weight, bool apply_median = false)
        {
            int x = left.GetLength(1);
            int y = left.GetLength(0);
            int z = 3;

            var gradient = new float[] { 1, -1 };
            var left_x = MatrixF.Conv(left, gradient, Direction.Horizontal);
            var left_y = MatrixF.Conv(left, gradient, Direction.Vertical);

            var right_x = MatrixF.Conv(right, gradient, Direction.Horizontal);
            var right_y = MatrixF.Conv(right, gradient, Direction.Vertical);

            var im_l = new float[][,] { left, left_x, left_y };
            var im_r = new float[][,] { right, right_x, right_y };

            var even_win = MathsF.IsEven(win) ? win : win + 1;
            var disparity = DisparityEstimator(im_l, im_r, even_win, max_dis, weight, x, y, z);

            if (apply_median)
            {
                // NOT USED
                disparity = MatrixF.Morph(disparity, even_win, even_win, MorphologyMode.Median);
            }

            return disparity;
        }
        /// <summary>
        /// Aggregates matching costs to produce the final disparity map.
        /// </summary>
        /// <param name="im_l">Left image stack</param>
        /// <param name="im_r">Right image stack</param>
        /// <param name="win">Window size</param>
        /// <param name="max_dis">Maximum disparity value</param>
        /// <param name="weight">Gradient weight</param>
        /// <param name="dim_x">Image width</param>
        /// <param name="dim_y">Image height</param>
        /// <param name="dim_z">Number of channels</param>
        /// <returns>Disparity matrix</returns>
        private float[,] DisparityEstimator(float[][,] im_l, float[][,] im_r, int win, int max_dis, float weight, int dim_x, int dim_y, int dim_z)
        {
            var disparity = new float[dim_y, dim_x];
            var min2_dis = new float[dim_y, dim_x].Add(float.MaxValue);

            //Parallel.For(0, max_dis, i =>
            for (int i = 0; i < max_dis; i++)
            {
                var min3_dis = new float[dim_z][,];

                // decomposition
                for (int z = 0; z < dim_z; z++)
                {
                    min3_dis[z] = new float[dim_y, dim_x];

                    for (int y = 0; y < dim_y; y++)
                    {
                        for (int x = 0; x < dim_x; x++)
                        {
                            var j = x - i;

                            if (j < 0)
                            {
                                min3_dis[z][y, x] = MathsF.Pow(im_l[z][y, x] - 0.0000000001f, 2);
                            }
                            else
                            {
                                min3_dis[z][y, x] = MathsF.Pow(im_l[z][y, x] - im_r[z][y, j], 2);
                            }
                        }
                    }
                }

                // mean gradients
                var flat_dis = new float[dim_y, dim_x];

                for (int y = 0; y < dim_y; y++)
                {
                    for (int x = 0; x < dim_x; x++)
                    {
                        flat_dis[y, x] = min3_dis[0][y, x] + weight * (min3_dis[1][y, x] + min3_dis[2][y, x]);
                    }
                }

                // box filter
                if (win > 1)
                {
                    flat_dis = MatrixF.Mean(flat_dis, win, win);
                }

                // mapping
                for (int y = 0; y < dim_y; y++)
                {
                    for (int x = 0; x < dim_x; x++)
                    {
                        if (flat_dis[y, x] < min2_dis[y, x])
                        {
                            disparity[y, x] = i;
                            min2_dis[y, x] = flat_dis[y, x];
                        }
                    }
                }
            }

            // normalizing
            return disparity.Div(max_dis);
        }
        #endregion
    }
}
