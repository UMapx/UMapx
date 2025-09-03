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
        private float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Laplace pyramid filter.
        /// </summary>
        /// <param name="lap">Laplacian pyramid</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public LaplacianPyramidFilter(LaplacianPyramidTransform lap, float factor = -1.0f)
        {
            this.LaplacianPyramid = lap;
            this.Factor = factor;
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

        #region Apply voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            // forward pyramid transform
            float[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = MatrixF.Mul(pA[i], 1.0f + this.factor / (1.0f + i));
            }

            // backward pyramid transform
            float[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(ComplexF[,] data)
        {
            // forward pyramid transform
            ComplexF[][,] pA = lap.Forward(data);

            int r = data.GetLength(0), c = data.GetLength(1);
            int nlev = pA.Length - 1, i, j;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = MatrixF.Mul(pA[i], 1.0f + this.factor / (1.0f + i));
            }

            // backward pyramid transform
            ComplexF[,] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                {
                    data[i, j] = dummy[i, j];
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            // forward pyramid transform
            float[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = MatrixF.Mul(pA[i], 1.0f + this.factor / (1.0f + i));
            }

            // backward pyramid transform
            float[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(ComplexF[] data)
        {
            // forward pyramid transform
            ComplexF[][] pA = lap.Forward(data);

            int r = data.GetLength(0);
            int nlev = pA.Length - 1, i;

            for (i = 0; i < nlev; i++)
            {
                pA[i] = MatrixF.Mul(pA[i], 1.0f + this.factor / (1.0f + i));
            }

            // backward pyramid transform
            ComplexF[] dummy = lap.Backward(pA);

            for (i = 0; i < r; i++)
            {
                data[i] = dummy[i];
            }
        }
        #endregion
    }
}
