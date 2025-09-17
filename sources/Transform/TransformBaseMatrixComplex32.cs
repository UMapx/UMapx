using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines an adapter class for a matrix transform.
    /// </summary>
    public abstract class TransformBaseMatrixComplex32 : TransformBase, ITransform
    {
        #region Matrix methods
        /// <summary>
        /// Implements the construction of the transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="backward">Backward transform or not</param>
        /// <returns>Matrix</returns>
        protected abstract Complex32[,] Matrix(int n, bool backward = false);
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            Complex32[,] U = Matrix(N, false);
            Complex32[] B = Matrice.Dot(A, U);

            if (Normalized)
            {
                B = B.Div(NormalizationFactor(N, false));
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            Complex32[,] U = Matrix(N, true);
            Complex32[] A = Matrice.Dot(B, Matrice.Hermitian(U));

            if (Normalized)
            {
                A = A.Div(NormalizationFactor(N, true));
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex32[,] U = Matrix(N, false);
            Complex32[,] V = Matrix(M, false);
            Complex32[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Hermitian());
                B = Normalized ? B.Div(NormalizationFactor(N * M, false)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = Normalized ? B.Div(NormalizationFactor(N, false)) : B;
            }
            else
            {
                B = A.Dot(V.Hermitian());
                B = Normalized ? B.Div(NormalizationFactor(N, false)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex32[,] U = Matrix(N, true);
            Complex32[,] V = Matrix(M, true);
            Complex32[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Hermitian().Dot(B).Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(N * M, true)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Hermitian().Dot(B);
                A = Normalized ? A.Div(NormalizationFactor(N, true)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(N, true)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
