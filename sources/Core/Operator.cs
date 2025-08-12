using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to implement convolutional operators. 
    /// </summary>
    public static class Operator
    {
        #region Fixed
        /// <summary>
        /// Implements the construction of the Roberts operator [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Roberts()
        {
            return new float[2, 2] { { 1, 0 }, { 0, -1 } };
        }
        /// <summary>
        /// Implements the construction of the Prewitt operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Prewitt()
        {
            return new float[3, 3] { { -1, -1, -1 }, { 0, 0, 0 }, { 1, 1, 1 } };
        }
        /// <summary>
        /// Implements the construction of the Sobel operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Sobel()
        {
            return new float[3, 3] { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };
        }
        /// <summary>
        /// Implements the construction of the Scharr operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Scharr()
        {
            return new float[3, 3] { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } };
        }
        /// <summary>
        /// Implements the construction of the Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Laplacian()
        {
            return new float[3, 3] { { 0, 1, 0 }, { 1, -4, 1 }, { 0, 1, 0 } };
        }
        /// <summary>
        /// Implements the construction of the diagonal Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] LaplacianDiagonal()
        {
            return new float[3, 3] { { 1, 1, 1 }, { 1, -8, 1 }, { 1, 1, 1 } };
        }
        /// <summary>
        /// Implements the construction of the inverted Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] LaplacianInvert()
        {
            return new float[3, 3] { { -1, 0, -1 }, { 0, 4, 0 }, { -1, 0, -1 } };
        }
        #endregion

        #region Compass fixed
        /// <summary>
        /// Implements the construction of the Roberts operator [3 x 3]. [2 x 2].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static float[,] Roberts(Gradient direction)
        {
            float[,] H = new float[2, 2];

            if (direction == Gradient.North)
            {
                H[0, 0] = 0; H[0, 1] = 1;
                H[1, 0] = -1; H[1, 1] = 0;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = -1;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = 0; H[0, 1] = -1;
                H[1, 0] = 1; H[1, 1] = 0;
            }
            else if (direction == Gradient.East)
            {
                H[0, 0] = -1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = 1;
            }
            else
            {
                throw new NotSupportedException("Direction is not supported in this operator");
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Kirsh operator [3 x 3].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static float[,] Kirsh(Gradient direction)
        {
            float[,] H = new float[3, 3];

            if (direction == Gradient.North)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.NorthWest)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 5; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.SouthWest)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = -3;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.SouthEast)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.East)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = 5;
            }
            else if (direction == Gradient.NorthEast)
            {
                H[0, 0] = -3; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }

            return H;
        }
        #endregion

        #region Radius
        /// <summary>
        /// Implements the construction of the inverted Gausssian filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static float[,] LoGaussian(int m, int l, float sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            float sigma2 = sigma * sigma;
            float f0 = -1.0f / (Maths.Pi * sigma2 * sigma2);
            float f1 = 2.0f * sigma2;
            float kernel;
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    kernel = (x * x + y * y) / f1;
                    H[i, j] = f0 * (1.0f - kernel) * Maths.Exp(-kernel);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of the Gaussian blur filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static float[,] Gaussian(int m, int l, float sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    H[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of the "unsharp masking" filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static float[,] Unsharp(int m, int l, float sigma)
        {
            float[,] G = new float[m, l];
            int i, j, x, y;
            int r1 = m / 2;
            int r2 = l / 2;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    G[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }

            float[,] invG = new float[m, l];
            float max = G[0, 0];
            float summary = 0;
            float v, iv;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    v = G[i, j] / max;
                    iv = v > byte.MaxValue ? byte.MaxValue : (int)v;
                    invG[i, j] = -iv;
                    summary += iv;
                }
            }

            invG[m / 2, l / 2] = 2 * summary - invG[m / 2, l / 2];
            return invG;
        }
        /// <summary>
        /// Implements the construction of the high-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="boost">Boost</param>
        /// <returns>Matrix</returns>
        public static float[,] HighPass(int m, int l, float boost)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = -1;
                }
            }
            H[r1, r2] = boost;
            return H;
        }
        /// <summary>
        /// Implements the construction of the low-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] LowPass(int m, int l)
        {
            return Operator.HighPass(m, l, 1);
        }
        /// <summary>
        /// Implements the construction of the emboss filter.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Emboss(int n)
        {
            float[,] H = new float[n, n];
            int r = n - 1, r2 = r / 2;

            H[0, 0] = -2; H[0, r2] = -1; //         0;
            H[r2, 0] = -1; H[r2, r2] = 1; H[r2, r] = 1;
            H[r, r2] = 1; H[r, r] = 2; //         0;

            return H;
        }
        #endregion
    }
}
