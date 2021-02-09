using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the multidimensional transform.
    /// </summary>
    [Serializable]
    public class MultidimensionalTransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the multidimensional transform.
        /// </summary>
        /// <param name="transform">ITransform</param>
        public MultidimensionalTransform(ITransform transform)
        {
            Transform = transform;
        }
        /// <summary>
        /// Gets or sets transform.
        /// </summary>
        public ITransform Transform { get; set; }
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public double[][] Forward(params double[][] A)
        {
            int count = A.Length;
            double[][] B = new double[count][];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public double[][,] Forward(params double[][,] A)
        {
            int count = A.Length;
            double[][,] B = new double[count][,];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public Complex[][] Forward(params Complex[][] A)
        {
            int count = A.Length;
            Complex[][] B = new Complex[count][];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public Complex[][,] Forward(params Complex[][,] A)
        {
            int count = A.Length;
            Complex[][,] B = new Complex[count][,];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public double[][] Backward(params double[][] B)
        {
            int count = B.Length;
            double[][] A = new double[count][];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="B">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public double[][,] Backward(params double[][,] B)
        {
            int count = B.Length;
            double[][,] A = new double[count][,];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public Complex[][] Backward(params Complex[][] B)
        {
            int count = B.Length;
            Complex[][] A = new Complex[count][];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        /// <summary>
        /// Forward multidimensional transform.
        /// </summary>
        /// <param name="B">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public Complex[][,] Backward(params Complex[][,] B)
        {
            int count = B.Length;
            Complex[][,] A = new Complex[count][,];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        #endregion
    }
}
