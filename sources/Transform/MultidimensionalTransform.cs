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
        public float[][] Forward(params float[][] A)
        {
            int count = A.Length;
            float[][] B = new float[count][];

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
        public float[][,] Forward(params float[][,] A)
        {
            int count = A.Length;
            float[][,] B = new float[count][,];

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
        public Complex32[][] Forward(params Complex32[][] A)
        {
            int count = A.Length;
            Complex32[][] B = new Complex32[count][];

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
        public Complex32[][,] Forward(params Complex32[][,] A)
        {
            int count = A.Length;
            Complex32[][,] B = new Complex32[count][,];

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
        public float[][] Backward(params float[][] B)
        {
            int count = B.Length;
            float[][] A = new float[count][];

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
        public float[][,] Backward(params float[][,] B)
        {
            int count = B.Length;
            float[][,] A = new float[count][,];

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
        public Complex32[][] Backward(params Complex32[][] B)
        {
            int count = B.Length;
            Complex32[][] A = new Complex32[count][];

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
        public Complex32[][,] Backward(params Complex32[][,] B)
        {
            int count = B.Length;
            Complex32[][,] A = new Complex32[count][,];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        #endregion
    }
}
