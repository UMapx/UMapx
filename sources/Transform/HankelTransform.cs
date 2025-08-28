using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the Hankel transform.
    /// <remarks>
    /// NOT RECOMMENDED.
    /// There is no fast O(N log N) algorithm for this transform.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Hankel_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HankelTransform : ITransform
    {
        #region Private data
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Param.
        /// </summary>
        private int a;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hankel transform.
        /// </summary>
        /// <param name="a">Param</param>
        /// <param name="direction">Processing direction</param>
        public HankelTransform(int a = 0, Direction direction = Direction.Vertical)
        {
            this.direction = direction;
            this.a = a;
        }
        /// <summary>
        /// Gets or sets the processing direction.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Gets or sets the param of transform.
        /// </summary>
        public int A
        {
            get
            { 
                return this.a; 
            }
            set
            { 
                this.a = value; 
            }
        }
        #endregion

        #region Hankel static components
        /// <summary>
        /// Implements the construction of the Hankel transform matrix.
        /// </summary>
        /// <param name="N">Size</param>
        /// <param name="a">Param</param>
        /// <returns>Matrix</returns>
        /// <exception cref="ArgumentException">Exception</exception>
        public static float[,] Matrix(int N, int a)
        {
            if (N <= 0 || a < 0) throw new ArgumentException("Arguments could not be negative");

            float[] j = new float[N + 1 + 1];

            for (int k = 1; k <= N + 1; k++)
                j[k] = BesselZeroJ(a, k);

            float[] Jp = new float[N + 1 + 1];

            for (int k = 1; k <= N + 1; k++)
                Jp[k] = Special.J(j[k], a + 1);

            float denom = j[N + 1];

            var T = new float[N, N];

            for (int m = 1; m <= N; m++)
            {
                for (int n = 1; n <= N; n++)
                {
                    float arg = j[m] * j[n] / denom;
                    float num = Special.J((float)arg, a);
                    float val = 2.0f / denom * num / (Jp[m] * Jp[n]);
                    T[m - 1, n - 1] = val;
                }
            }
            return T;
        }
        /// <summary>
        /// K-th positive zero of J_a (ν = a):
        /// Start from McMahon’s asymptotic and refine with Newton’s method.
        /// Do all computations in double for ~1e-12 relative accuracy of zeros.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="k">Value</param>
        /// <returns>Value</returns>
        private static float BesselZeroJ(int a, int k)
        {
            float nu = a;
            float x = (k + 0.5f * nu - 0.25f) * Maths.Pi;
            float eps = 1e-16f;
            int iterations = 120;

            for (int it = 0; it < iterations; it++)
            {
                float Ja = Special.J(x, a);
                float Ja1 = Special.J(x, a + 1);
                float Ja_1 = Special.J(x, a - 1);

                float Jprime = 0.5f * (Ja_1 - Ja1);
                float dx = Ja / Jprime;

                x -= dx;

                if (Math.Abs(dx) <= eps * Math.Abs(x)) break;
            }

            return x;
        }
        #endregion

        #region Hankel Transform
        /// <summary>
        /// Forward Hankel transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public float[] Forward(float[] A)
        {
            int N = A.Length;
            float[,] U = HankelTransform.Matrix(N, a);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward Hankel transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public float[] Backward(float[] B)
        {
            int N = B.Length;
            float[,] U = HankelTransform.Matrix(N, a);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward Hankel transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Forward(float[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = HankelTransform.Matrix(N, a);
            float[,] V = HankelTransform.Matrix(M, a);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Backward Hankel transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public float[,] Backward(float[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = HankelTransform.Matrix(N, a);
            float[,] V = HankelTransform.Matrix(M, a);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        /// <summary>
        /// Forward Hankel transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            float[,] U = HankelTransform.Matrix(N, a);
            return Matrice.Dot(A, U);
        }
        /// <summary>
        /// Backward Hankel transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length;
            float[,] U = HankelTransform.Matrix(N, a);
            return Matrice.Dot(B, U.Transponate());
        }
        /// <summary>
        /// Forward Hankel transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            float[,] U = HankelTransform.Matrix(N, a);
            float[,] V = HankelTransform.Matrix(M, a);

            if (direction == Direction.Both)
            {
                return U.Dot(A).Dot(V.Transponate());
            }
            else if (direction == Direction.Vertical)
            {
                return U.Dot(A);
            }
            return A.Dot(V.Transponate());
        }
        /// <summary>
        /// Backward Hankel transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex32[,] Backward(Complex32[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            float[,] U = HankelTransform.Matrix(N, a);
            float[,] V = HankelTransform.Matrix(M, a);

            if (direction == Direction.Both)
            {
                return U.Transponate().Dot(B).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                return U.Transponate().Dot(B);
            }
            return B.Dot(V);
        }
        #endregion
    }
}
