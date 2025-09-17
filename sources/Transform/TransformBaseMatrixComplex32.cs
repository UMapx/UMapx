using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines an adapter class for a matrix transform.
    /// </summary>
    public abstract class TransformBaseMatrixComplex32 : TransformBaseMatrix, ITransform
    {
        #region Matrix methods
        /// <summary>
        /// Implements the construction of the transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        protected abstract Complex32[,] TransformationMatrix(int n);
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public virtual Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            int Nu = ForwardMatrixSize(N);
            Complex32[,] U = TransformationMatrix(Nu);
            Complex32[] B = Matrice.Dot(A, U);

            if (Normalized)
            {
                B = B.Div(ForwardNormalizationFactor(Nu));
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public virtual Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            var Nu = BackwardMatrixSize(N);
            Complex32[,] U = TransformationMatrix(Nu);
            Complex32[] A = Matrice.Dot(B, Matrice.Hermitian(U));

            if (Normalized)
            {
                A = A.Div(BackwardNormalizationFactor(Nu));
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            var Nu = ForwardMatrixSize(N);
            var Mv = ForwardMatrixSize(M);

            Complex32[,] U = TransformationMatrix(Nu);
            Complex32[,] V = TransformationMatrix(Mv);
            Complex32[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Hermitian());
                B = Normalized ? B.Div(ForwardNormalizationFactor(Nu * Mv)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = Normalized ? B.Div(ForwardNormalizationFactor(Nu)) : B;
            }
            else
            {
                B = A.Dot(V.Hermitian());
                B = Normalized ? B.Div(ForwardNormalizationFactor(Mv)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            var Nu = ForwardMatrixSize(N);
            var Mv = ForwardMatrixSize(M);

            Complex32[,] U = TransformationMatrix(Nu);
            Complex32[,] V = TransformationMatrix(Mv);
            Complex32[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Hermitian().Dot(B).Dot(V);
                A = Normalized ? A.Div(BackwardNormalizationFactor(Nu * Mv)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Hermitian().Dot(B);
                A = Normalized ? A.Div(BackwardNormalizationFactor(Nu)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = Normalized ? A.Div(BackwardNormalizationFactor(Mv)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public virtual float[] Forward(float[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public virtual float[] Backward(float[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual float[,] Forward(float[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public virtual float[,] Backward(float[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
