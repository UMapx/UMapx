using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Gaussian pyramid transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class GaussianPyramidTransform : IPyramidTransform
    {
        #region Private data
        int radius;
        int levels;
        #endregion

        #region Pyramid components
        /// <summary>
        /// Initializes the Gaussian pyramid transform.
        /// </summary>
        public GaussianPyramidTransform()
        {
            this.Levels = int.MaxValue;
            this.Radius = 2;
        }
        /// <summary>
        /// Initializes the Gaussian pyramid transform.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="radius">Radius</param>
        public GaussianPyramidTransform(int levels, int radius = 2)
        {
            this.Levels = levels;
            this.Radius = radius;
        }
        /// <summary>
        /// Gets or sets number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.levels = value;
            }
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
        #endregion

        #region Apply voids
        // **************************************************
        //            Gaussian Pyramid Transform
        // **************************************************
        // ORIGINALS: Burt, P., and Adelson, E. H.
        // IEEE Transactions on Communication, COM-31:532-540 
        // (1983).
        // Designed by Valery Asiryan (c), 2015-2020
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public double[][,] Forward(double[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            double[][,] pyr = new double[nlev][,];
            double[,] dummy = (double[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public double[][] Forward(double[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            double[][] pyr = new double[nlev][];
            double[] dummy = (double[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public double[] Backward(double[][] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public Complex[][,] Forward(Complex[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            Complex[][,] pyr = new Complex[nlev][,];
            Complex[,] dummy = (Complex[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public Complex[][] Forward(Complex[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            Complex[][] pyr = new Complex[nlev][];
            Complex[] dummy = (Complex[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[][] pyramid)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Gaussian filter.
        /// </summary>
        /// <returns>Array</returns>
        public static double[] Filter
        {
            get
            {
                return new double[] { .0625, .25, .375, .25, .0625 };
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static double[,] upsample(double[,] u, int radius)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            double[,] v = new double[n, m];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                for (l = 0, j = 0; j < c; j++, l += 2)
                {
                    v[k + 1, l] = u[i, j];
                    v[k, l + 1] = u[i, j];
                    v[k, l] = u[i, j];
                    v[k + 1, l + 1] = u[i, j];
                }
            }

            return v.Mean(radius, radius);
        }
        /// <summary>
        ///  
        /// </summary>
        /// <param name="u">Array</param>
        /// <param name="radius">Radius</param>
        /// <returns>Array</returns>
        internal static double[] upsample(double[] u, int radius)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            double[] v = new double[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static double[,] downsample(double[,] u, int radius)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            double[,] v = new double[n, m];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                for (l = 0, j = 0; j < c; j += 2, l++)
                {
                    v[k, l] = u[i, j];
                }
            }

            return v.Mean(radius, radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static double[] downsample(double[] u, int radius)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            double[] v = new double[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static double[,] add(double[,] m, double[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static double[] add(double[] m, double[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            double[] v = new double[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static double[,] sub(double[,] m, double[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static double[] sub(double[] m, double[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            double[] v = new double[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] - n[i];
            }
            return v;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex[,] upsample(Complex[,] u, int radius)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            Complex[,] v = new Complex[n, m];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                for (l = 0, j = 0; j < c; j++, l += 2)
                {
                    v[k + 1, l] = u[i, j];
                    v[k, l + 1] = u[i, j];
                    v[k, l] = u[i, j];
                    v[k + 1, l + 1] = u[i, j];
                }
            }

            return v.Mean(radius, radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Array</param>
        /// <param name="radius">Radius</param>
        /// <returns>Array</returns>
        internal static Complex[] upsample(Complex[] u, int radius)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            Complex[] v = new Complex[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex[,] downsample(Complex[,] u, int radius)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            Complex[,] v = new Complex[n, m];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                for (l = 0, j = 0; j < c; j += 2, l++)
                {
                    v[k, l] = u[i, j];
                }
            }

            return v.Mean(radius, radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex[] downsample(Complex[] u, int radius)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            Complex[] v = new Complex[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static Complex[,] add(Complex[,] m, Complex[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex[,] H = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static Complex[] add(Complex[] m, Complex[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex[] v = new Complex[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static Complex[,] sub(Complex[,] m, Complex[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex[,] H = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static Complex[] sub(Complex[] m, Complex[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex[] v = new Complex[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] - n[i];
            }
            return v;
        }
        #endregion
    }
}
