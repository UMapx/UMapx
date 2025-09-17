using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines an adapter class for a matrix transform.
    /// </summary>
    public abstract class TransformBaseMatrixFloat : TransformBase, ITransform
    {
        #region Matrix methods
        /// <summary>
        /// Implements the construction of the transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        protected abstract float[,] TransformationMatrix(int n);
        #endregion

        #region Transform methods
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = TransformationMatrix(N);
            float[] B = Matrice.Dot(A, U);

            if (Normalized)
            {
                B = B.Div(NormalizationFactor(N));
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = TransformationMatrix(N);
            float[] A = Matrice.Dot(B, Matrice.Transpose(U));

            if (Normalized)
            {
                A = A.Div(NormalizationFactor(N));
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = TransformationMatrix(N);
            float[,] V = TransformationMatrix(M);
            float[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = Normalized ? B.Div(NormalizationFactor(N * M)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = Normalized ? B.Div(NormalizationFactor(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transpose());
                B = Normalized ? B.Div(NormalizationFactor(M)) : B;
            }

            return B;
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = TransformationMatrix(N);
            float[,] V = TransformationMatrix(M);
            float[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(N * M)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Transpose().Dot(B);
                A = Normalized ? A.Div(NormalizationFactor(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(M)) : A;
            }

            return A;
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = TransformationMatrix(N);
            Complex32[] B = Matrice.Dot(A, U);

            if (Normalized)
            {
                B = B.Div(NormalizationFactor(N));
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
            float[,] U = TransformationMatrix(N);
            Complex32[] A = Matrice.Dot(B, Matrice.Transpose(U));

            if (Normalized)
            {
                A = A.Div(NormalizationFactor(N));
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
            float[,] U = TransformationMatrix(N);
            float[,] V = TransformationMatrix(M);
            Complex32[,] B;

            if (Direction == Direction.Both)
            {
                B = U.Dot(A).Dot(V.Transpose());
                B = Normalized ? B.Div(NormalizationFactor(N * M)) : B;
            }
            else if (Direction == Direction.Vertical)
            {
                B = U.Dot(A);
                B = Normalized ? B.Div(NormalizationFactor(N)) : B;
            }
            else
            {
                B = A.Dot(V.Transpose());
                B = Normalized ? B.Div(NormalizationFactor(M)) : B;
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
            float[,] U = TransformationMatrix(N);
            float[,] V = TransformationMatrix(M);
            Complex32[,] A;

            if (Direction == Direction.Both)
            {
                A = U.Transpose().Dot(B).Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(N * M)) : A;
            }
            else if (Direction == Direction.Vertical)
            {
                A = U.Transpose().Dot(B);
                A = Normalized ? A.Div(NormalizationFactor(N)) : A;
            }
            else
            {
                A = B.Dot(V);
                A = Normalized ? A.Div(NormalizationFactor(M)) : A;
            }

            return A;
        }
        #endregion
    }
}
