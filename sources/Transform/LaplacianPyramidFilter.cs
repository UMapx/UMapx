using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Laplace pyramid filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cs.toronto.edu/~jepson/csc320/notes/pyramids.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class LaplacianPyramidFilter : IFilter
    {
        #region Private data
        private LaplacianPyramidTransform lap;
        private double factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Laplace pyramid filter.
        /// </summary>
        /// <param name="lap">Laplacian pyramid</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public LaplacianPyramidFilter(LaplacianPyramidTransform lap, double factor = -1.0)
        {
            this.lap = lap;
            this.factor = factor;
        }
        /// <summary>
        /// Gets or sets the Laplacian pyramid.
        /// </summary>
        public LaplacianPyramidTransform LaplacianPyramid
        {
            get
            {
                return lap;
            }
            set
            {
                lap = value;
            }
        }
        /// <summary>
        /// Gets or sets the factor [-1, 1].
        /// </summary>
        public double Factor
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

        #region Apply voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(double[,] data)
        {
            // forward pyramid transform
            double[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            double[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex[,] data)
        {
            // forward pyramid transform
            Complex[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            Complex[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(double[] data)
        {
            // forward pyramid transform
            double[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            double[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex[] data)
        {
            // forward pyramid transform
            Complex[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = Matrice.Mul(pA[i], 1.0 + this.factor);
            }

            // backward pyramid transform
            Complex[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }

            return;
        }
        #endregion
    }
}
