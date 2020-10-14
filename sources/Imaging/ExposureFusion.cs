using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the exposure fusion filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://web.stanford.edu/class/cs231m/project-1/exposure-fusion.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ExposureFusion
    {
        #region Private data
        private LaplacianPyramidTransform lap;
        private GaussianPyramidTransform gap;
        private double sigma;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the exposure fusion filter.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="sigma">Sigma (0, 1)</param>
        public ExposureFusion(int levels, double sigma = 0.2)
        {
            this.lap = new LaplacianPyramidTransform(levels);
            this.gap = new GaussianPyramidTransform(levels);
            this.sigma = sigma;
        }
        /// <summary>
        /// Gets or sets number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.lap.Levels;
            }
            set
            {
                this.lap.Levels =
                    this.gap.Levels = value;
            }
        }
        /// <summary>
        /// Gets or sets the sigma value (0, 1).
        /// </summary>
        public double Sigma
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
        /// Apply filter.
        /// </summary>
        /// <param name="images">Bitmap array</param>
        /// <returns>Bitmap</returns>
        public Bitmap Apply(params Bitmap[] images)
        {
            // data
            int N = images.GetLength(0);
            double[][][,] data = new double[N][][,];
            int height = images[0].Height;
            int width = images[1].Width;
            double[][,] weights = new double[N][,];

            // to rgb array
            for (int i = 0; i < N; i++)
            {
                data[i] = BitmapConverter.ToRGB(images[i]);
            }

            // initialize weights
            for (int i = 0; i < N; i++)
                weights[i] = Matrice.One(height, width);

            // applying params
            weights = Mul(weights, Exp(data, this.sigma));

            // normalizing
            double[,] z = new double[height, width];

            for (int i = 0; i < N; i++)
                z = z.Add(weights[i]);

            for (int i = 0; i < N; i++)
                weights[i] = weights[i].Div(z);

            // pyramids
            double[][,] pyrW;
            double[][,] pyrIr;
            double[][,] pyrIg;
            double[][,] pyrIb;

            // outputs
            double[][,] zero = gap.Forward(new double[height, width]);
            double[][,] r = (double[][,])zero.Clone();
            double[][,] g = (double[][,])zero.Clone();
            double[][,] b = (double[][,])zero.Clone();
            int levels = r.GetLength(0);

            // do job
            for (int i = 0; i < N; i++)
            {
                pyrW = gap.Forward(weights[i]);
                zero = data[i];
                pyrIr = lap.Forward(zero[0]);
                pyrIg = lap.Forward(zero[1]);
                pyrIb = lap.Forward(zero[2]);

                for (int l = 0; l < levels; l++)
                {
                    z = pyrW[l];

                    r[l] = r[l].Add(z.Mul(pyrIr[l]));
                    g[l] = g[l].Add(z.Mul(pyrIg[l]));
                    b[l] = b[l].Add(z.Mul(pyrIb[l]));
                }
            }

            // reconstruction
            Bitmap bitmap = BitmapConverter.FromRGB(new double[][,] {
                lap.Backward(r),
                lap.Backward(g),
                lap.Backward(b) });

            return bitmap;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="images">BitmapData array</param>
        /// <returns>Bitmap</returns>
        public Bitmap Apply(params BitmapData[] images)
        {
            // data
            int N = images.GetLength(0);
            double[][][,] data = new double[N][][,];
            int height = images[0].Height;
            int width = images[1].Width;
            double[][,] weights = new double[N][,];

            // to rgb array
            for (int i = 0; i < N; i++)
            {
                data[i] = BitmapConverter.ToRGB(images[i]);
            }

            // initialize weights
            for (int i = 0; i < N; i++)
                weights[i] = Matrice.One(height, width);

            // applying params
            weights = Mul(weights, Exp(data, this.sigma));

            // normalizing
            double[,] z = new double[height, width];

            for (int i = 0; i < N; i++)
                z = z.Add(weights[i]);

            for (int i = 0; i < N; i++)
                weights[i] = weights[i].Div(z);

            // pyramids
            double[][,] pyrW;
            double[][,] pyrIr;
            double[][,] pyrIg;
            double[][,] pyrIb;

            // outputs
            double[][,] zero = gap.Forward(new double[height, width]);
            double[][,] r = (double[][,])zero.Clone();
            double[][,] g = (double[][,])zero.Clone();
            double[][,] b = (double[][,])zero.Clone();
            int levels = r.GetLength(0);

            // do job
            for (int i = 0; i < N; i++)
            {
                pyrW = gap.Forward(weights[i]);
                zero = data[i];
                pyrIr = lap.Forward(zero[0]);
                pyrIg = lap.Forward(zero[1]);
                pyrIb = lap.Forward(zero[2]);

                for (int l = 0; l < levels; l++)
                {
                    z = pyrW[l];

                    r[l] = r[l].Add(z.Mul(pyrIr[l]));
                    g[l] = g[l].Add(z.Mul(pyrIg[l]));
                    b[l] = b[l].Add(z.Mul(pyrIb[l]));
                }
            }

            // reconstruction
            Bitmap bitmap = BitmapConverter.FromRGB(new double[][,] {
                lap.Backward(r),
                lap.Backward(g),
                lap.Backward(b) });

            return bitmap;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Exponent filter.
        /// </summary>
        /// <param name="I">Input data</param>
        /// <param name="sigma">Sigma</param>
        /// <returns>Output data</returns>
        private static double[][,] Exp(double[][][,] I, double sigma)
        {
            // params
            int length = I.GetLength(0);
            double[][,] C = new double[length][,];
            double[][,] current;
            double[,] r, g, b;
            double[,] result;
            int height, width;
            int i, j, k;

            // do job
            for (i = 0; i < length; i++)
            {
                current = I[i];
                r = current[0];
                g = current[1];
                b = current[2];

                height = r.GetLength(0);
                width = r.GetLength(1);
                result = new double[height, width];

                for (j = 0; j < height; j++)
                {
                    for (k = 0; k < width; k++)
                    {
                        result[j, k] = Math.Exp(-0.5 * Math.Pow((r[j, k] - 0.5), 2) / Math.Pow(sigma, 2)) *
                                       Math.Exp(-0.5 * Math.Pow((g[j, k] - 0.5), 2) / Math.Pow(sigma, 2)) *
                                       Math.Exp(-0.5 * Math.Pow((b[j, k] - 0.5), 2) / Math.Pow(sigma, 2));
                    }
                }
                C[i] = result;
            }

            return C;
        }
        /// <summary>
        /// Implements matrix array multiplication.
        /// </summary>
        /// <param name="A">Matrix array</param>
        /// <param name="B">Matrix array</param>
        /// <returns>Matrix array</returns>
        private static double[][,] Mul(double[][,] A, double[][,] B)
        {
            int length = A.GetLength(0);
            double[][,] R = new double[length][,];

            for (int i = 0; i < length; i++)
                R[i] = A[i].Mul(B[i]);

            return R;
        }
        #endregion
    }
}
