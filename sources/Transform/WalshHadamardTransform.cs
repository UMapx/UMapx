using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Walsh-Hadamard transform.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://kibia.ru/teachers/kreindelin/pdf/2.pdf
    /// </remarks>
    [Serializable]
    public class WalshHadamardTransform : TransformBaseMatrixFloat, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes the Walsh-Hadamard transform.
        /// </summary>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
        public WalshHadamardTransform(bool normalized = true, Direction direction = Direction.Vertical)
        {
            this.Normalized = normalized; 
            this.Direction = direction;
        }
        #endregion

        #region Walsh-Hadamard static components
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix.
        /// </summary>
        /// <param name="powOf2">Power of 2</param>
        /// <returns>Matrix</returns>
        public static float[,] Matrix(int powOf2)
        {
            if (powOf2 < 0)
                throw new ArgumentOutOfRangeException(nameof(powOf2), "Power must be non-negative");

            if (powOf2 == 0)
                return new float[1, 1] { { 1 } };

            int iterations = powOf2 - 1;
            float[,] hadamard = WalshHadamardTransform.Matrix();

            if (iterations > 0)
            {
                float[,] H = WalshHadamardTransform.Matrix();

                for (int i = 0; i < iterations; i++)
                {
                    H = Core.Matrice.Kronecker(H, hadamard);
                }
                return H;
            }

            return hadamard;
        }
        /// <summary>
        /// Implements the construction of the Walsh-Hadamard matrix [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static float[,] Matrix()
        {
            return new float[2, 2] { { 1, 1 }, { 1, -1 } };
        }
        #endregion

        #region Walsh-Hadamard Transform
        /// <inheritdoc/>
        protected override float ForwardNormalizationFactor(int n)
        {
            return Maths.Sqrt(n);
        }
        /// <inheritdoc/>
        protected override float BackwardNormalizationFactor(int n)
        {
            return Maths.Sqrt(n);
        }
        /// <inheritdoc/>
        protected override float[,] TransformationMatrix(int n)
        {
            if (!Maths.IsPower(n, 2))
                throw new ArgumentException("Dimension of the signal must be a power of 2");

            int N = (int)Maths.Log2(n);
            return Matrix(N);
        }
        #endregion
    }
}
