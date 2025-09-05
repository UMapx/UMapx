using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Gaussian pyramid transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </remarks>
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
                    throw new ArgumentException("Invalid argument value");

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
        // **************************************************

        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public float[][,] Forward(float[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            float[][,] pyr = new float[nlev][,];
            float[,] dummy = (float[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = Downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public float[][] Forward(float[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            float[][] pyr = new float[nlev][];
            float[] dummy = (float[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = Downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public Complex32[][,] Forward(Complex32[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c))
                / Math.Log(2)), levels);

            Complex32[][,] pyr = new Complex32[nlev][,];
            Complex32[,] dummy = (Complex32[,])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = Downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Forward Gaussian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public Complex32[][] Forward(Complex32[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            Complex32[][] pyr = new Complex32[nlev][];
            Complex32[] dummy = (Complex32[])data.Clone();

            for (int i = 0; i < nlev; i++)
            {
                pyr[i] = dummy;
                dummy = Downsample(dummy, this.radius);
            }

            return pyr;
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] pyramid)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Gaussian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] pyramid)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Upsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static float[,] Upsample(float[,] u, int radius)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            float[,] v = new float[n, m];

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
        /// Upsample the input signal.
        /// </summary>
        /// <param name="u">Array</param>
        /// <param name="radius">Radius</param>
        /// <returns>Array</returns>
        internal static float[] Upsample(float[] u, int radius)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            float[] v = new float[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// Downsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static float[,] Downsample(float[,] u, int radius)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            float[,] v = new float[n, m];

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
        /// Downsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static float[] Downsample(float[] u, int radius)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            float[] v = new float[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// Add two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static float[,] Add(float[,] m, float[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            float[,] H = new float[ml, mr];
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
        /// Add two arrays.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static float[] Add(float[] m, float[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            float[] v = new float[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// Sub two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static float[,] Sub(float[,] m, float[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            float[,] H = new float[ml, mr];
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
        /// Sub two arrays.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static float[] Sub(float[] m, float[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            float[] v = new float[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] - n[i];
            }
            return v;
        }
        /// <summary>
        /// Upsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex32[,] Upsample(Complex32[,] u, int radius)
        {
            int r = u.GetLength(0), c = u.GetLength(1);
            int n = r * 2, m = c * 2;
            int i, j, k, l;
            Complex32[,] v = new Complex32[n, m];

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
        /// Upsample the input signal.
        /// </summary>
        /// <param name="u">Array</param>
        /// <param name="radius">Radius</param>
        /// <returns>Array</returns>
        internal static Complex32[] Upsample(Complex32[] u, int radius)
        {
            int r = u.GetLength(0);
            int n = r * 2;
            int i, k;
            Complex32[] v = new Complex32[n];

            for (k = 0, i = 0; i < r; i++, k += 2)
            {
                v[k] = u[i];
                v[k + 1] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// Downsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex32[,] Downsample(Complex32[,] u, int radius)
        {
            int r = u.GetLength(0);
            int c = u.GetLength(1);
            int n = (r + 1) / 2, m = (c + 1) / 2;
            int i, j, k, l;
            Complex32[,] v = new Complex32[n, m];

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
        /// Downsample the input signal.
        /// </summary>
        /// <param name="u">Matrix</param>
        /// <param name="radius">Radius</param>
        /// <returns>Matrix</returns>
        internal static Complex32[] Downsample(Complex32[] u, int radius)
        {
            int r = u.Length;
            int n = (r + 1) / 2;
            int i, k;
            Complex32[] v = new Complex32[n];

            for (k = 0, i = 0; i < r; i += 2, k++)
            {
                v[k] = u[i];
            }

            return v.Mean(radius);
        }
        /// <summary>
        /// Add two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static Complex32[,] Add(Complex32[,] m, Complex32[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex32[,] H = new Complex32[ml, mr];
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
        /// Add two arrays.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static Complex32[] Add(Complex32[] m, Complex32[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex32[] v = new Complex32[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                v[i] = m[i] + n[i];
            }
            return v;
        }
        /// <summary>
        /// Sub two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        internal static Complex32[,] Sub(Complex32[,] m, Complex32[,] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            int mr = (int)Math.Min(m.GetLength(1), n.GetLength(1));
            Complex32[,] H = new Complex32[ml, mr];
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
        /// Sub two arrays.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="n">Array</param>
        /// <returns>Array</returns>
        internal static Complex32[] Sub(Complex32[] m, Complex32[] n)
        {
            int ml = (int)Math.Min(m.GetLength(0), n.GetLength(0));
            Complex32[] v = new Complex32[ml];
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
