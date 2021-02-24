using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Laplacian pyramid transform.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LaplacianPyramidTransform : IPyramidTransform
    {
        #region Private data
        int radius;
        int levels;
        #endregion

        #region Pyramid components
        /// <summary>
        /// Initializes the Laplacian pyramid transform.
        /// </summary>
        public LaplacianPyramidTransform()
        {
            this.Levels = int.MaxValue;
            this.Radius = 2;
        }
        /// <summary>
        /// Initializes the Laplacian pyramid transform.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="radius">Radius</param>
        public LaplacianPyramidTransform(int levels, int radius = 2)
        {
            this.Radius = radius;
            this.Levels = levels;
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
        //            Laplacian Pyramid Transform
        // **************************************************
        // ORIGINALS: Burt, P., and Adelson, E. H.
        // IEEE Transactions on Communication, COM-31:532-540 
        // (1983).
        // Designed by Valery Asiryan (c), 2015-2020
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public float[][,] Forward(float[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c)) / Math.Log(2)), levels);
            float[][,] lapl = new float[nlev][,];
            float[,] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J, this.radius);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I, this.radius));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public float[][] Forward(float[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            float[][] lapl = new float[nlev][];
            float[] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J, this.radius);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I, this.radius));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Backward Laplacian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[][,] pyramid)
        {
            int nlev = pyramid.Length - 1;
            float[,] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I, this.radius));
            }

            return I;
        }
        /// <summary>
        /// Backward Laplacian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public float[] Backward(float[][] pyramid)
        {
            int nlev = pyramid.Length;
            float[] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I, this.radius));
            }

            return I;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <returns>Pyramid</returns>
        public Complex32[][,] Forward(Complex32[,] data)
        {
            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = (int)Math.Min((Math.Log(Math.Min(r, c)) / Math.Log(2)), levels);
            Complex32[][,] lapl = new Complex32[nlev][,];
            Complex32[,] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J, this.radius);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I, this.radius));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Pyramid</returns>
        public Complex32[][] Forward(Complex32[] data)
        {
            int r = data.Length;
            int nlev = (int)Math.Min((Math.Log(r) / Math.Log(2)), levels);

            Complex32[][] lapl = new Complex32[nlev][];
            Complex32[] I, J = data;

            for (int i = 0; i < nlev - 1; i++)
            {
                I = GaussianPyramidTransform.downsample(J, this.radius);
                lapl[i] = GaussianPyramidTransform.sub(J, GaussianPyramidTransform.upsample(I, this.radius));
                J = I;
            }

            lapl[nlev - 1] = J;
            return lapl;
        }
        /// <summary>
        /// Backward Laplacian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[][,] pyramid)
        {
            int nlev = pyramid.Length - 1;
            Complex32[,] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I, this.radius));
            }

            return I;
        }
        /// <summary>
        /// Backward Laplacian pyramid transform.
        /// </summary>
        /// <param name="pyramid">Pyramid</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[][] pyramid)
        {
            int nlev = pyramid.Length;
            Complex32[] I = pyramid[nlev];

            for (int i = nlev - 1; i >= 0; i--)
            {
                I = GaussianPyramidTransform.add(pyramid[i], GaussianPyramidTransform.upsample(I, this.radius));
            }

            return I;
        }
        #endregion

        #region Gaussian pyramid to Laplacian pyramid
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Gaussian pyramid</param>
        /// <returns>Pyramid</returns>
        public float[][,] Forward(float[][,] data)
        {
            int nlev = data.Length;
            float[][,] lapl = new float[nlev][,];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i], this.radius));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Gaussian pyramid</param>
        /// <returns>Pyramid</returns>
        public float[][] Forward(float[][] data)
        {
            int nlev = data.Length;
            float[][] lapl = new float[nlev][];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i], this.radius));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Gaussian pyramid</param>
        /// <returns>Pyramid</returns>
        public Complex32[][,] Forward(Complex32[][,] data)
        {
            int nlev = data.Length;
            Complex32[][,] lapl = new Complex32[nlev][,];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i], this.radius));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        /// <summary>
        /// Forward Laplacian pyramid transform.
        /// </summary>
        /// <param name="data">Gaussian pyramid</param>
        /// <returns>Pyramid</returns>
        public Complex32[][] Forward(Complex32[][] data)
        {
            int nlev = data.Length;
            Complex32[][] lapl = new Complex32[nlev][];

            for (int i = 1; i < nlev; i++)
            {
                lapl[i - 1] = GaussianPyramidTransform.sub(data[i - 1], GaussianPyramidTransform.upsample(data[i], this.radius));
            }

            lapl[nlev - 1] = data[nlev - 1];
            return lapl;
        }
        #endregion
    }
}
