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
        /// Returns normalization factor.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="backward">Backward transform or not</param>
        /// <returns>Value</returns>
        public virtual float NormalizationFactor(int n, bool backward = false)
        {
            return 1.0f;
        }
        /// <summary>
        /// Implements the construction of the transform matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <param name="backward">Backward transform or not</param>
        /// <returns>Matrix</returns>
        public abstract float[,] Matrix(int n, bool backward = false);
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
            float[,] U = Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = Matrix(N, true);
            return Matrice.Dot(B, Matrice.Transpose(U));
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = Matrix(N);
            float[,] V = Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transpose());
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = Matrix(N, true);
            float[,] V = Matrix(M, true);

            if (Direction == Direction.Both)
            {
                return U.Transpose().Dot(B).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Transpose().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = Matrix(N);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = Matrix(N, true);
            return Matrice.Dot(B, Matrice.Transpose(U));
        }
        /// <summary>
        /// Forward transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = Matrix(N);
            float[,] V = Matrix(M);

            if (Direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transpose());
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transpose());
        }
        /// <summary>
        /// Backward transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = Matrix(N, true);
            float[,] V = Matrix(M, true);

            if (Direction == Direction.Both)
            {
                return U.Transpose().Dot(B).Dot(V);
            }
            else if (Direction == Direction.Vertical)
            {
                return U.Transpose().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion
    }
}
