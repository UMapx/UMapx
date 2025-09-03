using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the multidimensional pyramid transform.
    /// </summary>
    [Serializable]
    public class MultidimensionalPyramidTransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the multidimensional pyramid transform.
        /// </summary>
        /// <param name="transform">IPyramidTransform</param>
        public MultidimensionalPyramidTransform(IPyramidTransform transform)
        {
            Transform = transform;
        }
        /// <summary>
        /// Gets or sets pyramid transform.
        /// </summary>
        public IPyramidTransform Transform { get; set; }
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public float[][][] Forward(params float[][] A)
        {
            int count = A.Length;
            float[][][] B = new float[count][][];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public float[][][,] Forward(params float[][,] A)
        {
            int count = A.Length;
            float[][][,] B = new float[count][][,];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <returns>Jagged array</returns>
        public ComplexF[][][] Forward(params ComplexF[][] A)
        {
            int count = A.Length;
            ComplexF[][][] B = new ComplexF[count][][];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="A">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public ComplexF[][][,] Forward(params ComplexF[][,] A)
        {
            int count = A.Length;
            ComplexF[][][,] B = new ComplexF[count][][,];

            for (int i = 0; i < count; i++)
            {
                B[i] = Transform.Forward(A[i]);
            }

            return B;
        }
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public float[][] Backward(params float[][][] B)
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
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="B">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public float[][,] Backward(params float[][][,] B)
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
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public ComplexF[][] Backward(params ComplexF[][][] B)
        {
            int count = B.Length;
            ComplexF[][] A = new ComplexF[count][];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        /// <summary>
        /// Forward multidimensional pyramid transform.
        /// </summary>
        /// <param name="B">Jagged matrix</param>
        /// <returns>Jagged matrix</returns>
        public ComplexF[][,] Backward(params ComplexF[][][,] B)
        {
            int count = B.Length;
            ComplexF[][,] A = new ComplexF[count][,];

            for (int i = 0; i < count; i++)
            {
                A[i] = Transform.Backward(B[i]);
            }

            return A;
        }
        #endregion
    }
}
