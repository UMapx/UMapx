using System;
using System.Drawing;
using System.Threading.Tasks;

namespace UMapx.Core
{
    /// <summary>
    /// Used to implement standard algebraic operations on matrices and vectors.
    /// </summary>
    public static class Matrice
    {
        // Matrix voids

        #region Matrix booleans
        /// <summary>
        /// Checks the equality of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsEquals(this float[,] m, float[,] n)
        {
            int r = m.GetLength(0);
            int c = m.GetLength(1);

            if (r != n.GetLength(0) || c != n.GetLength(1))
                return false;

            int i, j;

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                    if (m[i, j] != n[i, j])
                        return false;
            }

            return true;
        }
        /// <summary>
        /// Checks the equality of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsEquals(this Complex32[,] m, Complex32[,] n)
        {
            int r = m.GetLength(0);
            int c = m.GetLength(1);

            if (r != n.GetLength(0) || c != n.GetLength(1))
                return false;

            int i, j;

            for (i = 0; i < r; i++)
            {
                for (j = 0; j < c; j++)
                    if (m[i, j] != n[i, j])
                        return false;
            }

            return true;
        }
        /// <summary>
        /// Checks if the matrix is a vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsVector(this float[,] m)
        {
            if (m.GetLength(0) == 1 || m.GetLength(1) == 1)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is square.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSquare(this float[,] m)
        {
            if (m.GetLength(0) == m.GetLength(1))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is positive.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsPositive(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    if (m[i, j] < 0)
                        return false;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if the matrix is symmetric.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSymmetric(this float[,] m)
        {
            if (Matrice.IsSquare(m))
            {
                // ?A = A'
                if (Matrice.IsEquals(m, m.Transponate()))
                {
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is skew-symmetric.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSkewSymmetric(this float[,] m)
        {
            if (Matrice.IsSquare(m))
            {
                // ?A' = -A:
                if (Matrice.IsEquals(m.Transponate(), m.ToNegate()))
                {
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is diagonal.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsDiagonal(this float[,] m)
        {
            int i, j;
            int ml = m.GetLength(0), mr = m.GetLength(1);

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    if (i != j)
                    {
                        if (m[i, j] != 0)
                        {
                            return false;
                        }
                    }
                    else continue;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if the matrix is a vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsVector(this Complex32[,] m)
        {
            if (m.GetLength(0) == 1 || m.GetLength(1) == 1)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is square.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSquare(this Complex32[,] m)
        {
            if (m.GetLength(0) == m.GetLength(1))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is symmetric (Hermitian).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSymmetric(this Complex32[,] m)
        {
            if (Matrice.IsSquare(m))
            {
                // ?A = A'
                if (Matrice.IsEquals(m, Matrice.Hermitian(m)))
                {
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is skew-symmetric (anti-Hermitian).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsSkewSymmetric(this Complex32[,] m)
        {
            if (Matrice.IsSquare(m))
            {
                // ?A' = -A
                if (Matrice.IsEquals(m.Hermitian(), m.ToNegate()))
                {
                    return true;
                }
            }
            return false;
        }
        /// <summary>
        /// Checks if the matrix is diagonal.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool IsDiagonal(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    if (i != j)
                    {
                        if (m[i, j] != 0)
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
        #endregion

        #region Matrix transform
        /// <summary>
        /// Implements the matrix inversion operation.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Invert(this float[,] m)
        {
            if (m.GetLength(0) != m.GetLength(1))
            {
                return m;
            }

            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Invert(Jagged.ToJagged(m)));
        }
        /// <summary>
        /// Implements the transpose of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Transponate(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r1, r0];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[j, i] = m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the matrix inversion operation.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Invert(this Complex32[,] m)
        {
            if (m.GetLength(0) != m.GetLength(1))
            {
                return m;
            }

            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Invert(Jagged.ToJagged(m)));
        }
        /// <summary>
        /// Implements the transpose of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Transponate(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r1, r0];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[j, i] = m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the complex conjugate matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conjugate(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Conjugate;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the Hermitian-conjugation operation of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Hermitian(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r1, r0];
            int i, j, x, y;

            for (i = 0, x = 0; (i < r0) && (x < r0); i++, x++)
            {
                for (j = 0, y = 0; (j < r1) && (y < r1); j++, y++)
                {
                    H[y, x] = m[i, j].Conjugate;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a Gram (Hermitian) matrix.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns></returns>
        public static Complex32[,] Gram(this Complex32[,] A)
        {
            int n = A.GetLength(0), m = A.GetLength(1);
            var H = new Complex32[m, m];

            for (int i = 0; i < m; i++)
            {
                for (int j = i; j < m; j++)
                {
                    float re = 0f, im = 0f;
                    for (int k = 0; k < n; k++)
                    {
                        var ai = A[k, i];  // a = x + i y
                        var aj = A[k, j];  // b = u + i v
                        // conj(a)*b = (x - i y)(u + i v) = (x u + y v) + i (x v - y u)
                        re += ai.Real * aj.Real + ai.Imag * aj.Imag;
                        im += ai.Real * aj.Imag - ai.Imag * aj.Real;
                    }
                    H[i, j] = new Complex32(re, im);
                    if (i != j) H[j, i] = new Complex32(re, -im);
                }
            }
            return H;
        }
        /// <summary>
        /// Returns a Gram matrix.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Gram(this float[,] A)
        {
            int n = A.GetLength(0), m = A.GetLength(1);
            var G = new float[m, m];

            for (int i = 0; i < m; i++)
            {
                for (int j = i; j < m; j++)
                {
                    double sum = 0.0; // accumulate in double for better numeric stability
                    for (int k = 0; k < n; k++)
                        sum += (double)A[k, i] * A[k, j];

                    float v = (float)sum;
                    G[i, j] = v;
                    if (i != j) G[j, i] = v; // mirror to keep symmetry
                }
            }
            return G;
        }
        #endregion

        #region Matrix properties
        /// <summary>
        /// Returns the trace value of a square matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Value</returns>
        public static float Trace(this float[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new ArgumentException("The matrix must be square");

            int d = m.GetLength(0);
            int i;
            float kernel = 0;

            for (i = 0; i < d; i++)
            {
                kernel += m[i, i];
            }
            return kernel;
        }
        /// <summary>
        /// Returns the value of the matrix determinant.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Value</returns>
        public static float Det(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new ArgumentException("The matrix must be square");

            unsafe
            {
                // copy array
                float[,] n = (float[,])m.Clone();

                fixed (float* pm = &n[0, 0])
                    return LinealgOptions.MatrixOperation.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Returns the P-norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Matrix</returns>
        public static float Norm(this float[,] m, float p)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float norm = 0;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    norm += (float)Math.Pow(Math.Abs(m[i, j]), p);
                }
            }

            return (float)Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float Norm(this float[,] m)
        {
            return Matrice.Norm(m, 2);
        }
        /// <summary>
        /// Selects the integer part of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="digits">Digits</param>
        /// <param name="mode">Midpoint rounding</param>
        /// <returns>Matrix</returns>
        public static float[,] Round(this float[,] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] H = new float[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = (float)Math.Round(m[i, j], digits, mode);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a square permutation matrix.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <returns>Square matrix</returns>
        public static float[,] Permutation(this float[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new ArgumentException("The matrix must be square");

            int i, j, r, n = m.GetLength(0);
            float[] temp; float diagonal;
            float[][] perm = Jagged.ToJagged(Matrice.Eye(n, n));

            for (i = 0; i < n; i++)
            {
                diagonal = m[i, i]; r = i;

                for (j = i; j < n; j++)
                {
                    if (m[j, i] > diagonal)
                    {
                        diagonal = m[j, i]; r = j;
                    }
                }

                if (i != r)
                {
                    temp = perm[i];
                    perm[i] = perm[r];
                    perm[r] = temp;
                }
            }
            return Jagged.FromJagged(perm);
        }
        /// <summary>
        /// Returns the trace value of a square matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Complex number</returns>
        public static Complex32 Trace(this Complex32[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new ArgumentException("The matrix must be square");

            int d = m.GetLength(0);
            int i;
            Complex32 kernel = 0;

            for (i = 0; i < d; i++)
            {
                kernel += m[i, i];
            }
            return kernel;
        }
        /// <summary>
        /// Returns the value of the matrix determinant.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Value</returns>
        public static Complex32 Det(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new ArgumentException("The matrix must be square");

            unsafe
            {
                // copy array
                Complex32[,] n = (Complex32[,])m.Clone();

                fixed (Complex32* pm = &n[0, 0])
                    return LinealgOptions.MatrixOperation.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Returns the P-norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Matrix</returns>
        public static float Norm(this Complex32[,] m, float p)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float norm = 0;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    norm += (float)Math.Pow(m[i, j].Abs, p);
                }
            }

            return (float)Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float Norm(this Complex32[,] m)
        {
            return Matrice.Norm(m, 2);
        }
        /// <summary>
        /// Selects the integer part of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="digits">Digits</param>
        /// <param name="mode">Midpoint rounding</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Round(this Complex32[,] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            Complex32 c;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    c = m[i, j];
                    H[i, j] = new Complex32((float)Math.Round(c.Real, digits, mode), (float)Math.Round(c.Imag, digits, mode));
                }
            }

            return H;
        }
        #endregion

        #region Matrix kronecker product
        /// <summary>
        /// Returns the Kronecker matrix product.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Kronecker(this float[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            float[,] H = new float[ml * nl, mr * nr];
            int i, j, k, l;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    for (k = 0; k < nl; k++)
                    {
                        for (l = 0; l < nr; l++)
                        {
                            H[i * nl + k, j * nr + l] = m[i, j] * n[k, l];
                        }
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the Kronecker matrix product.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Kronecker(this Complex32[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex32[,] H = new Complex32[ml * nl, mr * nr];
            int i, j, k, l;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    for (k = 0; k < nl; k++)
                    {
                        for (l = 0; l < nr; l++)
                        {
                            H[i * nl + k, j * nr + l] = m[i, j] * n[k, l];
                        }
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the Kronecker matrix product.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Kronecker(this Complex32[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex32[,] H = new Complex32[ml * nl, mr * nr];
            int i, j, k, l;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    for (k = 0; k < nl; k++)
                    {
                        for (l = 0; l < nr; l++)
                        {
                            H[i * nl + k, j * nr + l] = m[i, j] * n[k, l];
                        }
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the Kronecker matrix product.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Kronecker(this float[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex32[,] H = new Complex32[ml * nl, mr * nr];
            int i, j, k, l;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    for (k = 0; k < nl; k++)
                    {
                        for (l = 0; l < nr; l++)
                        {
                            H[i * nl + k, j * nr + l] = m[i, j] * n[k, l];
                        }
                    }
                }
            }

            return H;
        }
        #endregion

        #region Matrix add/sub
        /// <summary>
        /// Returns the sum of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Add(this float[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] H = new float[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Returns the sum of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this Complex32[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Returns the sum of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this Complex32[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Returns the sum of two matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this float[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] + n[i, j];
                }
            }
            return H;
        }

        /// <summary>
        /// Returns the sum of a matrix and a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Add(this float[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] + a;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a matrix and a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this Complex32[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] + a;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a matrix and a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this Complex32[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] + a;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a matrix and a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(this float[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] + a;
                }
            }

            return H;
        }

        /// <summary>
        /// Returns the sum of a number and a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Add(float a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a + m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a number and a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(Complex32 a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a + m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a number and a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(Complex32 a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a + m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns the sum of a number and a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Add(float a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a + m[i, j];
                }
            }

            return H;
        }

        /// <summary>
        /// Subtracts one matrix from another.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Sub(this float[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] H = new float[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Subtracts one matrix from another.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this Complex32[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Subtracts one matrix from another.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this Complex32[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Subtracts one matrix from another.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this float[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = m[i, j] - n[i, j];
                }
            }
            return H;
        }

        /// <summary>
        /// Subtracts a number from the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Sub(this float[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] - a;
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a number from the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this Complex32[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] - a;
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a number from the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this Complex32[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] - a;
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a number from the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(this float[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] - a;
                }
            }

            return H;
        }

        /// <summary>
        /// Subtracts a matrix from a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Sub(float a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a - m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a matrix from a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(Complex32 a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a - m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a matrix from a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(Complex32 a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a - m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Subtracts a matrix from a number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Sub(float a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a - m[i, j];
                }
            }

            return H;
        }
        #endregion

        #region Matrix mul
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Vector</param>
        /// <returns>Matrix</returns>
        public static float[,] Mul(this float[,] m, float[] v)
        {
            int mr = m.GetLength(0), ml = m.GetLength(1);
            var H = new float[mr, ml];

            for (int j = 0; j < ml; j++)
            {
                var s = v[j];
                for (int i = 0; i < mr; i++)
                {
                    var a = m[i, j];
                    H[i, j] = a * s;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Vector</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, Complex32[] v)
        {
            int mr = m.GetLength(0), ml = m.GetLength(1);
            var H = new Complex32[mr, ml];

            for (int j = 0; j < ml; j++)
            {
                var s = v[j];
                for (int i = 0; i < mr; i++)
                {
                    var a = m[i, j];
                    H[i, j] = a * s;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Vector</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, float[] v)
        {
            int mr = m.GetLength(0), ml = m.GetLength(1);
            var H = new Complex32[mr, ml];

            for (int j = 0; j < ml; j++)
            {
                var s = v[j];
                for (int i = 0; i < mr; i++)
                {
                    var a = m[i, j];
                    H[i, j] = a * s;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Vector</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this float[,] m, Complex32[] v)
        {
            int mr = m.GetLength(0), ml = m.GetLength(1);
            var H = new Complex32[mr, ml];

            for (int j = 0; j < ml; j++)
            {
                var s = v[j];
                for (int i = 0; i < mr; i++)
                {
                    var a = m[i, j];
                    H[i, j] = a * s;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Mul(this float[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            float[,] H = new float[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] * n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] * n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] * n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements matrix multiplication.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this float[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] * n[i, j];
                }
            }
            return H;
        }

        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Mul(this float[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this Complex32[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(this float[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }

        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Mul(float a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a * m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(Complex32 a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(Complex32 a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        /// <summary>
        /// Multiplies all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Mul(float a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] * a;
                }
            }

            return H;
        }
        #endregion

        #region Matrix div
        /// <summary>
        /// Divides a matrix by a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Div(this float[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            float[,] H = new float[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] / n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Divides a matrix by a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this Complex32[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] / n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Divides a matrix by a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this Complex32[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] / n[i, j];
                }
            }
            return H;
        }
        /// <summary>
        /// Divides a matrix by a matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this float[,] m, Complex32[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = m[i, j] / n[i, j];
                }
            }
            return H;
        }

        /// <summary>
        /// Divides all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Div(this float[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] / a;
                }
            }

            return H;
        }
        /// <summary>
        /// Divides all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this Complex32[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] / a;
                }
            }

            return H;
        }
        /// <summary>
        /// Divides all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this Complex32[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] / a;
                }
            }

            return H;
        }
        /// <summary>
        /// Divides all matrix elements by number.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(this float[,] m, Complex32 a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j] / a;
                }
            }

            return H;
        }

        /// <summary>
        /// Divides number into matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Div(float a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a / m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Divides number into matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(Complex32 a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a / m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Divides number into matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(Complex32 a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a / m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Divides number into matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Div(float a, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = a / m[i, j];
                }
            }

            return H;
        }
        #endregion

        #region Matrix pow
        /// <summary>
        /// Raises all matrix elements to a power.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="pow">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Pow(this Complex32[,] m, float pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Pow(m[i, j], pow);
                }
            }

            return H;
        }
        /// <summary>
        /// Raises all matrix elements to a power.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="pow">Value</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Pow(this float[,] m, Complex32 pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Pow(m[i, j], pow);
                }
            }

            return H;
        }
        /// <summary>
        /// Raises all matrix elements to a power.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="pow">Value</param>
        /// <returns>Matrix</returns>
        public static float[,] Pow(this float[,] m, float pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = (float)Math.Pow(m[i, j], pow);
                }
            }

            return H;
        }

        /// <summary>
        /// Raises the number to the power of the matrix.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Pow(float a, float[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = (float)Math.Pow(a, m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Raises the number to the power of the matrix.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Pow(Complex32 a, float[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Pow(a, m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Raises the number to the power of the matrix.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Pow(float a, Complex32[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Pow(a, m[i, j]);
                }
            }

            return H;
        }
        #endregion

        #region Matrix conversions
        /// <summary>
        /// Negates all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToNegate(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = -m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Negates all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] ToNegate(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = -m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a complex matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] ToComplex(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] H = new Complex32[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a matrix whose values belong to the interval [0, 255].
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToByte(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Byte(m[i, j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Returns a matrix whose values belong to the interval [0, 1].
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToFloat(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Float(m[i, j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToAbs(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = (float)Maths.Abs(m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToAbs(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Abs;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes an angle for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToAngle(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Angle;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes the real part for all elements of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToReal(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Real;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes the imaginary part for all elements of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] ToImag(this Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] H = new float[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Imag;
                }
            }

            return H;
        }
        #endregion

        #region Matrix statistics
        /// <summary>
        /// Returns the vector of matrix sums.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Sum(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] kernel = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    kernel[i] += m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the vector of matrix sums.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Sum(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[] kernel = new Complex32[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    kernel[i] += m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the matrix product vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Mul(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] kernel = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    kernel[i] *= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the matrix product vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[] kernel = new Complex32[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    kernel[i] *= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the matrix divide vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Div(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] kernel = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    kernel[i] /= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the matrix divide vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[] kernel = new Complex32[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    kernel[i] /= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Returns the vector of the matrix mode.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Mode(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] h = new float[mr];
            float[] kernel = new float[ml];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    kernel[j] = m[j, i];
                }

                h[i] = Matrice.Mode(kernel);
            }
            return h;
        }
        /// <summary>
        /// Returns the vector of the matrix mode.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Mode(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[] h = new Complex32[mr];
            Complex32[] kernel = new Complex32[ml];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    kernel[j] = m[j, i];
                }

                h[i] = Matrice.Mode(kernel);
            }
            return h;
        }
        /// <summary>
        /// Sorts the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        public static float[,] Sort(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] h = new float[ml, mr];
            float[] kernel = new float[ml];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    kernel[j] = m[j, i];
                }

                Array.Sort(kernel);

                for (j = 0; j < ml; j++)
                {
                    h[j, i] = kernel[j];
                }
            }
            return h;
        }
        /// <summary>
        /// Returns the maximum matrix vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Max(this float[,] m)
        {
            return Max(m, out _);
        }
        /// <summary>
        /// Returns the maximum matrix vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="index">Index array</param>
        /// <returns>Array</returns>
        public static float[] Max(this float[,] m, out int[] index)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] v = new float[mr];
            index = new int[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                v[i] = float.MinValue;
            }

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (m[j, i] > v[i])
                    {
                        v[i] = m[j, i];
                        index[i] = j;
                    }
                }
            }
            return v;
        }
        /// <summary>
        /// Returns the minimum matrix vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Min(this float[,] m)
        {
            return Min(m, out _);
        }
        /// <summary>
        /// Returns the minimum matrix vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="index">Index array</param>
        /// <returns>Array</returns>
        public static float[] Min(this float[,] m, out int[] index)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] v = new float[mr];
            index = new int[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                v[i] = float.MaxValue;
            }

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (m[j, i] < v[i])
                    {
                        v[i] = m[j, i];
                        index[i] = j;
                    }
                }
            }
            return v;
        }
        /// <summary>
        /// Returns the matrix vector corresponding to the specified morphology mode.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="mode">Morphology mode</param>
        /// <returns>Array</returns>
        public static float[] Morph(this float[,] m, MorphologyMode mode = MorphologyMode.Median)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[] v = new float[ml];
            float[] u = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[j] = m[j, i];
                }

                Array.Sort(v);
                u[i] = v[LinealgOptions.MorphologyFilter.GetFilterRank(mode, v.Length)];
            }
            return u;
        }
        /// <summary>
        /// Returns the vector of means of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Mean(this float[,] m)
        {
            return Matrice.Div(Matrice.Sum(m), m.GetLength(0));
        }
        /// <summary>
        /// Returns the vector of means of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Mean(this Complex32[,] m)
        {
            return Matrice.Div(Matrice.Sum(m), m.GetLength(0));
        }
        /// <summary>
        /// Returns the vector of variances of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Var(this float[,] m)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            float[] u = Matrice.Mean(m);
            float[] v = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += (float)Math.Pow(m[j, i] - u[i], 2);
                }
            }
            return Matrice.Div(v, ml - 1);
        }
        /// <summary>
        /// Returns the vector of variances of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Var(this Complex32[,] m)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            Complex32[] u = Matrice.Mean(m);
            Complex32[] v = new Complex32[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += Maths.Pow(m[j, i] - u[i], 2);
                }
            }
            return Matrice.Div(v, ml - 1);
        }
        /// <summary>
        /// Returns the vector of variances of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Var(this float[,] m, float[,] n)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            float[] v = new float[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += (float)Math.Pow(m[j, i] - n[j, i], 2);
                }
            }
            return Matrice.Div(v, ml - 1);
        }
        /// <summary>
        /// Returns the vector of variances of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Var(this Complex32[,] m, Complex32[,] n)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            Complex32[] v = new Complex32[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += Maths.Pow(m[j, i] - n[j, i], 2);
                }
            }
            return Matrice.Div(v, ml - 1);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] StnDev(this float[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5f);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] StnDev(this Complex32[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5f);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static float[] StnDev(this float[,] m, float[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5f);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] StnDev(this Complex32[,] m, Complex32[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5f);
        }
        /// <summary>
        /// Returns the covariance matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Cov(this float[,] m)
        {
            float[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            float[,] H = (float[,])m.Clone();
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    H[i, j] -= v[j];
                }
            }

            return H.Transponate().Dot(H).Div(height - 1);
        }
        /// <summary>
        /// Returns the covariance matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Cov(this Complex32[,] m)
        {
            Complex32[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            Complex32[,] H = (Complex32[,])m.Clone();
            int i, j;

            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    H[i, j] -= v[j];
                }
            }

            return H.Hermitian().Dot(H).Div(height - 1);
        }
        /// <summary>
        /// Returns the entropy vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Entropy(this float[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            float[] v = new float[width];
            int i, j;

            for (i = 0; i < width; i++)
            {
                for (j = 0; j < height; j++)
                {
                    if (m[j, i] > 0)
                    {
                        v[i] += -m[j, i] * (float)Maths.Log2(m[j, i]);
                    }
                }
            }
            return v;
        }
        /// <summary>
        /// Returns the normalized matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Normalized(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] u = new float[ml, mr];
            var min = m.Min().Min();
            var max = m.Max().Max();
            var dm = max - min;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    u[i, j] = (m[i, j] - min) / dm;
                }
            }

            return u;
        }
        #endregion

        #region Matrix concatenation
        /// <summary>
        /// Implements matrix concatenation.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="B">Matrix</param>
        /// <param name="direction">Direction</param>
        /// <returns>Matrix</returns>
        public static float[,] Concat(this float[,] A, float[,] B, Direction direction = Direction.Horizontal)
        {
            int aRows = A.GetLength(0), aCols = A.GetLength(1);
            int bRows = B.GetLength(0), bCols = B.GetLength(1);

            switch (direction)
            {
                case Direction.Horizontal:
                    if (aRows != bRows)
                        throw new ArgumentException("For horizontal concat, row counts must match");
                    {
                        var R = new float[aRows, aCols + bCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, 0, aCols);
                        return R;
                    }

                case Direction.Vertical:
                    if (aCols != bCols)
                        throw new ArgumentException("For vertical concat, column counts must match");
                    {
                        var R = new float[aRows + bRows, aCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, aRows, 0);
                        return R;
                    }

                case Direction.Both:
                default:
                    throw new NotSupportedException("Both direction does not have an unambiguous meaning for matrix concatenation");
            }
        }
        /// <summary>
        /// Implements matrix concatenation.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="B">Matrix</param>
        /// <param name="direction">Direction</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Concat(this Complex32[,] A, Complex32[,] B, Direction direction)
        {
            int aRows = A.GetLength(0), aCols = A.GetLength(1);
            int bRows = B.GetLength(0), bCols = B.GetLength(1);

            switch (direction)
            {
                case Direction.Horizontal:
                    if (aRows != bRows)
                        throw new ArgumentException("For horizontal concat, row counts must match");
                    {
                        var R = new Complex32[aRows, aCols + bCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, 0, aCols);
                        return R;
                    }

                case Direction.Vertical:
                    if (aCols != bCols)
                        throw new ArgumentException("For vertical concat, column counts must match");
                    {
                        var R = new Complex32[aRows + bRows, aCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, aRows, 0);
                        return R;
                    }

                case Direction.Both:
                default:
                    throw new NotSupportedException("Both direction does not have an unambiguous meaning for matrix concatenation");
            }
        }
        /// <summary>
        /// Implements matrix concatenation.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="B">Matrix</param>
        /// <param name="direction">Direction</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Concat(this Complex32[,] A, float[,] B, Direction direction)
        {
            int aRows = A.GetLength(0), aCols = A.GetLength(1);
            int bRows = B.GetLength(0), bCols = B.GetLength(1);

            switch (direction)
            {
                case Direction.Horizontal:
                    if (aRows != bRows)
                        throw new ArgumentException("For horizontal concat, row counts must match");
                    {
                        var R = new Complex32[aRows, aCols + bCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, 0, aCols);
                        return R;
                    }

                case Direction.Vertical:
                    if (aCols != bCols)
                        throw new ArgumentException("For vertical concat, column counts must match");
                    {
                        var R = new Complex32[aRows + bRows, aCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, aRows, 0);
                        return R;
                    }

                case Direction.Both:
                default:
                    throw new NotSupportedException("Both direction does not have an unambiguous meaning for matrix concatenation");
            }
        }
        /// <summary>
        /// Implements matrix concatenation.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <param name="B">Matrix</param>
        /// <param name="direction">Direction</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Concat(this float[,] A, Complex32[,] B, Direction direction)
        {
            int aRows = A.GetLength(0), aCols = A.GetLength(1);
            int bRows = B.GetLength(0), bCols = B.GetLength(1);

            switch (direction)
            {
                case Direction.Horizontal:
                    if (aRows != bRows)
                        throw new ArgumentException("For horizontal concat, row counts must match");
                    {
                        var R = new Complex32[aRows, aCols + bCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, 0, aCols);
                        return R;
                    }

                case Direction.Vertical:
                    if (aCols != bCols)
                        throw new ArgumentException("For vertical concat, column counts must match");
                    {
                        var R = new Complex32[aRows + bRows, aCols];
                        LinealgOptions.MatrixOperation.Copy(A, R, 0, 0);
                        LinealgOptions.MatrixOperation.Copy(B, R, aRows, 0);
                        return R;
                    }

                case Direction.Both:
                default:
                    throw new NotSupportedException("Both direction does not have an unambiguous meaning for matrix concatenation");
            }
        }
        #endregion

        // Matrix special

        #region Matrix dot
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Dot(this float[,] m, float[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Dot(this Complex32[,] m, Complex32[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Dot(this Complex32[,] m, float[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Dot(this float[,] m, Complex32[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.MatrixOperation.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        #endregion

        #region Matrix convolutions
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static float[,] Conv(this float[,] m, float[,] n, bool normalize = true)
        {
            return LinealgOptions.ConvolutionFilter.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this Complex32[,] m, Complex32[,] n, bool normalize = true)
        {
            return LinealgOptions.ConvolutionFilter.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this Complex32[,] m, float[,] n, bool normalize = true)
        {
            return LinealgOptions.ConvolutionFilter.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this float[,] m, Complex32[,] n, bool normalize = true)
        {
            return LinealgOptions.ConvolutionFilter.Conv(m, n, normalize);
        }
        #endregion

        #region Matrix convolutions (separable)
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static float[,] Conv(this float[,] m, float[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvolutionFilter.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvolutionFilter.ConvVertical(
                LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this float[,] m, Complex32[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvolutionFilter.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvolutionFilter.ConvVertical(
                LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this Complex32[,] m, float[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvolutionFilter.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvolutionFilter.ConvVertical(
                LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Conv(this Complex32[,] m, Complex32[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvolutionFilter.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvolutionFilter.ConvVertical(
                LinealgOptions.ConvolutionFilter.ConvHorizontal(m, n, normalize), n, normalize);
        }
        #endregion

        #region Matrix morphology
        /// <summary>
        /// Returns the matrix result of morphological minimum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        /// <returns>Matrix</returns>
        public static float[,] Min(this float[,] m, int r0, int r1)
        {
            return Morph(m, r0, r1, MorphologyMode.Erosion);
        }
        /// <summary>
        /// Returns the matrix result of morphological maximum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        /// <returns>Matrix</returns>
        public static float[,] Max(this float[,] m, int r0, int r1)
        {
            return Morph(m, r0, r1, MorphologyMode.Dilatation);
        }
        /// <summary>
        /// Returns the matrix result of morphology.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        /// <param name="mode">Morphology mode</param>
        /// <returns>Matrix</returns>
        public static float[,] Morph(this float[,] m, int r0, int r1, MorphologyMode mode = MorphologyMode.Median)
        {
            return LinealgOptions.MorphologySortFilter.Apply(m, r0 / 2, r1 / 2, mode);
        }
        #endregion

        #region Matrix mean (separable)
        /// <summary>
        /// Returns the result matrix of local averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static float[,] Mean(this float[,] m, int r0, int r1)
        {
            return LinealgOptions.MeanFilter.MeanVertical(LinealgOptions.MeanFilter.MeanHorizontal(m, r1), r0);
        }
        /// <summary>
        /// Returns the result matrix of local averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static Complex32[,] Mean(this Complex32[,] m, int r0, int r1)
        {
            return LinealgOptions.MeanFilter.MeanVertical(LinealgOptions.MeanFilter.MeanHorizontal(m, r1), r0);
        }

        /// <summary>
        /// Returns the result matrix of local weighted averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="w">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static float[,] Mean(this float[,] m, float[,] w, int r0, int r1)
        {
            return LinealgOptions.MeanFilter.MeanVerticalWeighted(LinealgOptions.MeanFilter.MeanHorizontalWeighted(m, w, r1), w, r0);
        }
        /// <summary>
        /// Returns the result matrix of local weighted averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="w">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static Complex32[,] Mean(this Complex32[,] m, Complex32[,] w, int r0, int r1)
        {
            return LinealgOptions.MeanFilter.MeanVerticalWeighted(LinealgOptions.MeanFilter.MeanHorizontalWeighted(m, w, r1), w, r0);
        }
        #endregion

        // Vector voids

        #region Vector transform
        /// <summary>
        /// Implements the diagonal vector inversion operation.
        /// </summary>
        /// <param name="v">Vector</param>
        /// <returns>Vector</returns>
        public static float[] Invert(this float[] v)
        {
            var inv = new float[v.Length];

            for (int i = 0; i < v.Length; i++) 
                inv[i] = (v[i] > 0f) ? 1f / v[i] : 0f;

            return inv;
        }
        /// <summary>
        /// Implements the diagonal vector inversion operation.
        /// </summary>
        /// <param name="v">Vector</param>
        /// <returns>Vector</returns>
        public static Complex32[] Invert(this Complex32[] v)
        {
            var inv = new Complex32[v.Length];

            for (int i = 0; i < v.Length; i++)
            {
                float x = v[i].Real;
                float y = v[i].Imag;
                float den = x * x + y * y; // |z|^2

                inv[i] = (den != 0f)
                    ? new Complex32(x / den, -y / den) // conj(z) / |z|^2
                    : Complex32.Zero;
            }
            return inv;
        }
        #endregion

        #region Vector booleans
        /// <summary>
        /// Checks the equality of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsEquals(this float[] a, float[] b)
        {
            int n = a.Length;

            if (n != b.Length)
                return false;

            for (int i = 0; i < n; ++i)
            {
                if (a[i] != b[i])
                    return false;
            }
            return true;
        }
        /// <summary>
        /// Checks the equality of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsEquals(this Complex32[] a, Complex32[] b)
        {
            int n = a.Length;

            if (n != b.Length)
                return false;

            for (int i = 0; i < n; ++i)
            {
                if (a[i] != b[i])
                    return false;
            }
            return true;
        }
        /// <summary>
        /// Checks if vectors are collinear.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsPositive(this float[] v)
        {
            int N = v.Length;

            for (int i = 0; i < N; i++)
            {
                if (v[i] < 0)
                {
                    return false;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if vectors are collinear.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsCollinear(this float[] a, float[] b)
        {
            int N = a.Length, i, j;
            float k;

            for (i = 0; i < N; ++i)
            {
                if (a[i] == 0 &&
                    b[i] == 0) continue;

                k = a[i] / b[i];

                for (j = i; j < N; j++)
                {
                    if (a[j] != b[j] * k) return false;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if vectors are collinear.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsCollinear(this Complex32[] a, Complex32[] b)
        {
            int N = a.Length, i, j;
            Complex32 k;

            for (i = 0; i < N; ++i)
            {
                if (a[i] == 0 &&
                    b[i] == 0) continue;

                k = a[i] / b[i];

                for (j = i; j < N; j++)
                {
                    if (a[j] != b[j] * k) return false;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if vectors are collinear.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsCollinear(this Complex32[] a, float[] b)
        {
            int N = a.Length, i, j;
            Complex32 k;

            for (i = 0; i < N; ++i)
            {
                if (a[i] == 0 &&
                    b[i] == 0) continue;

                k = a[i] / b[i];

                for (j = i; j < N; j++)
                {
                    if (a[j] != b[j] * k) return false;
                }
            }
            return true;
        }
        /// <summary>
        /// Checks if vectors are collinear.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Boolean</returns>
        public static bool IsCollinear(this float[] a, Complex32[] b)
        {
            int N = a.Length, i, j;
            Complex32 k;

            for (i = 0; i < N; ++i)
            {
                if (a[i] == 0 &&
                    b[i] == 0) continue;

                k = a[i] / b[i];

                for (j = i; j < N; j++)
                {
                    if (a[j] != b[j] * k) return false;
                }
            }
            return true;
        }
        #endregion

        #region Vector properties
        /// <summary>
        /// Returns the P-norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Value</returns>
        public static float Norm(this float[] a, float p)
        {
            int length = a.Length, i;
            float norm = 0;

            for (i = 0; i < length; i++)
            {
                norm += (float)Math.Pow(Math.Abs(a[i]), p);
            }
            return (float)Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <returns>Value</returns>
        public static float Norm(this float[] a)
        {
            return Norm(a, 2);
        }
        /// <summary>
        /// Returns the P-norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Value</returns>
        public static float Norm(this Complex32[] a, float p)
        {
            int length = a.Length, i;
            float norm = 0;

            for (i = 0; i < length; i++)
            {
                norm += (float)Maths.Pow(Maths.Abs(a[i]), p);
            }
            return (float)Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <returns>Value</returns>
        public static float Norm(this Complex32[] a)
        {
            return Norm(a, 2);
        }
        /// <summary>
        /// Selects the integer part of the matrix.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="digits">Digits</param>
        /// <param name="mode">Midpoint rounding</param>
        /// <returns>Array</returns>
        public static float[] Round(this float[] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0);
            float[] H = new float[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                H[i] = (float)Math.Round(H[i], digits, mode);
            }

            return H;
        }
        /// <summary>
        /// Selects the integer part of the matrix.
        /// </summary>
        /// <param name="m">Array</param>
        /// <param name="digits">Digits</param>
        /// <param name="mode">Midpoint rounding</param>
        /// <returns>Array</returns>
        public static Complex32[] Round(this Complex32[] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0);
            Complex32[] H = new Complex32[ml];
            Complex32 c;
            int i;

            for (i = 0; i < ml; i++)
            {
                c = m[i];
                H[i] = new Complex32((float)Math.Round(c.Real, digits, mode), (float)Math.Round(c.Imag, digits, mode));
            }

            return H;
        }
        #endregion

        #region Vector angle, projection, cosine
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static float Angle(this float[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Angle(this Complex32[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Angle(this float[] a, Complex32[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Angle(this Complex32[] a, Complex32[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }

        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static float Proj(this float[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Proj(this Complex32[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Proj(this float[] a, Complex32[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Proj(this Complex32[] a, Complex32[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }

        /// <summary>
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
        public static float[] GramProj(this float[] e, float[] a)
        {
            int length = e.Length;
            float[] proj = new float[length];
            int i;
            float ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * a[i];
                ee += e[i] * e[i];
            }

            float div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }
        /// <summary>
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] GramProj(this Complex32[] e, Complex32[] a)
        {
            int length = e.Length;
            Complex32[] proj = new Complex32[length];
            int i;
            Complex32 ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * a[i].Conjugate;
                ee += e[i] * e[i].Conjugate;
            }

            Complex32 div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }
        /// <summary>
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] GramProj(this float[] e, Complex32[] a)
        {
            int length = e.Length;
            Complex32[] proj = new Complex32[length];
            int i;
            Complex32 ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * a[i].Conjugate;
                ee += e[i] * e[i];
            }

            Complex32 div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }
        /// <summary>
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] GramProj(this Complex32[] e, float[] a)
        {
            int length = e.Length;
            Complex32[] proj = new Complex32[length];
            int i;
            Complex32 ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * a[i];
                ee += e[i] * e[i].Conjugate;
            }

            Complex32 div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }

        /// <summary>
        /// Returns the direction cosines of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] Cosines(this float[] v)
        {
            int length = v.Length, i;
            float abs = Matrice.Norm(v);
            float[] cos = new float[length];

            for (i = 0; i < length; i++)
            {
                cos[i] = v[i] / abs;
            }
            return cos;
        }
        /// <summary>
        /// Returns the direction cosines of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Cosines(this Complex32[] v)
        {
            int length = v.Length, i;
            Complex32 abs = Matrice.Norm(v);
            Complex32[] cos = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                cos[i] = v[i] / abs;
            }
            return cos;
        }
        #endregion

        #region Vector add/sub
        /// <summary>
        /// Returns the sum of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Add(this float[] a, float[] b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this Complex32[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this Complex32[] a, float[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this float[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b[i];
            }
            return c;
        }

        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static float[] Add(this float[] a, float b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b;
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this Complex32[] a, Complex32 b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b;
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this Complex32[] a, float b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b;
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(this float[] a, Complex32 b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] + b;
            }
            return c;
        }

        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static float[] Add(float b, float[] a)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b + a[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(Complex32 b, Complex32[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b + a[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(float b, Complex32[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b + a[i];
            }
            return c;
        }
        /// <summary>
        /// Returns the sum of a vector and a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Add(Complex32 b, float[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b + a[i];
            }
            return c;
        }

        /// <summary>
        /// Subtracts one vector from another.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this Complex32[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts one vector from another.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Sub(this float[] a, float[] b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts one vector from another.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this Complex32[] a, float[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts one vector from another.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this float[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b[i];
            }
            return c;
        }

        /// <summary>
        /// Subtracts a number from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static float[] Sub(this float[] a, float b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b;
            }
            return c;
        }
        /// <summary>
        /// Subtracts a number from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this Complex32[] a, Complex32 b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b;
            }
            return c;
        }
        /// <summary>
        /// Subtracts a number from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this Complex32[] a, float b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b;
            }
            return c;
        }
        /// <summary>
        /// Subtracts a number from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(this float[] a, Complex32 b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] - b;
            }
            return c;
        }

        /// <summary>
        /// Subtracts a vector from a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static float[] Sub(float b, float[] a)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b - a[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts a vector from a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(Complex32 b, Complex32[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b - a[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts a vector from a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(Complex32 b, float[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b - a[i];
            }
            return c;
        }
        /// <summary>
        /// Subtracts a vector from a number.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Sub(float b, Complex32[] a)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = b - a[i];
            }
            return c;
        }
        #endregion

        #region Vector mul
        /// <summary>
        /// Implements element-wise product of vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Mul(this float[] a, float[] b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] * b[i];
            }
            return c;
        }
        /// <summary>
        /// Implements element-wise product of vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this Complex32[] a, float[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] * b[i];
            }
            return c;
        }
        /// <summary>
        /// Implements element-wise product of vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this float[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] * b[i];
            }
            return c;
        }
        /// <summary>
        /// Implements element-wise product of vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this Complex32[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] * b[i];
            }
            return c;
        }

        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static float[] Mul(this float[] v, float a)
        {
            int length = v.Length, i;
            float[] H = new float[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] * a;
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this float[] v, Complex32 a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] * a;
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this Complex32[] v, float a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] * a;
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(this Complex32[] v, Complex32 a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] * a;
            }
            return H;
        }

        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static float[] Mul(float a, float[] v)
        {
            int length = v.Length, i;
            float[] H = new float[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a * v[i];
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(Complex32 a, Complex32[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a * v[i];
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(Complex32 a, float[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a * v[i];
            }
            return H;
        }
        /// <summary>
        /// Implements the multiplication of the vector by number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Mul(float a, Complex32[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a * v[i];
            }
            return H;
        }
        #endregion

        #region Vector div
        /// <summary>
        /// Divides a vector by a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Div(this float[] a, float[] b)
        {
            int length = a.Length, i;
            float[] c = new float[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] / b[i];
            }
            return c;
        }
        /// <summary>
        /// Divides a vector by a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this Complex32[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] / b[i];
            }
            return c;
        }
        /// <summary>
        /// Divides a vector by a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this Complex32[] a, float[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] / b[i];
            }
            return c;
        }
        /// <summary>
        /// Divides a vector by a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this float[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32[] c = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                c[i] = a[i] / b[i];
            }
            return c;
        }

        /// <summary>
        /// Divides a vector by a number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static float[] Div(this float[] v, float a)
        {
            int length = v.Length, i;
            float[] H = new float[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] / a;
            }
            return H;
        }
        /// <summary>
        /// Divides a vector by a number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this float[] v, Complex32 a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] / a;
            }
            return H;
        }
        /// <summary>
        /// Divides a vector by a number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this Complex32[] v, float a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] / a;
            }
            return H;
        }
        /// <summary>
        /// Divides a vector by a number.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(this Complex32[] v, Complex32 a)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = v[i] / a;
            }
            return H;
        }

        /// <summary>
        /// Divides a number by a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static float[] Div(float a, float[] v)
        {
            int length = v.Length, i;
            float[] H = new float[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a / v[i];
            }
            return H;
        }
        /// <summary>
        /// Divides a number by a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(Complex32 a, Complex32[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a / v[i];
            }
            return H;
        }
        /// <summary>
        /// Divides a number by a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(float a, Complex32[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a / v[i];
            }
            return H;
        }
        /// <summary>
        /// Divides a number by a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Value</param>
        /// <returns>Array</returns>
        public static Complex32[] Div(Complex32 a, float[] v)
        {
            int length = v.Length, i;
            Complex32[] H = new Complex32[length];

            for (i = 0; i < length; i++)
            {
                H[i] = a / v[i];
            }
            return H;
        }
        #endregion

        #region Vector pow
        /// <summary>
        /// Raises the elements of a vector to a power.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="power">Power</param>
        /// <returns>Array</returns>
        public static float[] Pow(this float[] v, float power)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = (float)Math.Pow(v[i], power);
            }
            return H;
        }
        /// <summary>
        /// Raises the elements of a vector to a power.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="power">Power</param>
        /// <returns>Array</returns>
        public static Complex32[] Pow(this Complex32[] v, float power)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Maths.Pow(v[i], power);
            }
            return H;
        }
        /// <summary>
        /// Raises the elements of a vector to a power.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="power">Power</param>
        /// <returns>Array</returns>
        public static Complex32[] Pow(this float[] v, Complex32 power)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Maths.Pow(v[i], power);
            }
            return H;
        }

        /// <summary>
        /// Raises the number to the power of the vector.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] Pow(float a, float[] v)
        {
            int n = v.GetLength(0);
            float[] H = new float[n];
            int i;

            for (i = 0; i < n; i++)
            {
                H[i] = (float)Math.Pow(a, v[i]);
            }

            return H;
        }
        /// <summary>
        /// Raises the number to the power of the vector.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Pow(Complex32 a, float[] v)
        {
            int n = v.GetLength(0);
            Complex32[] H = new Complex32[n];
            int i;

            for (i = 0; i < n; i++)
            {
                H[i] = Maths.Pow(a, v[i]);
            }

            return H;
        }
        /// <summary>
        /// Raises the number to the power of the vector.
        /// </summary>
        /// <param name="a">Value</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Pow(float a, Complex32[] v)
        {
            int n = v.GetLength(0);
            Complex32[] H = new Complex32[n];
            int i;

            for (i = 0; i < n; i++)
            {
                H[i] = Maths.Pow(a, v[i]);
            }

            return H;
        }
        #endregion

        #region Vector conversions
        /// <summary>
        /// Returns a vector whose values belong to the interval [0, 1].
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToFloat(this float[] v)
        {
            int length = v.Length;
            float max = Matrice.Max(v);
            float min = Matrice.Min(v);
            float range = max - min;
            float[] u = new float[length];

            for (int i = 0; i < length; i++)
            {
                u[i] = (v[i] - min) / range;
            }
            return u;
        }
        /// <summary>
        /// Returns a vector whose values belong to the interval [0, 255].
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToByte(this float[] v)
        {
            int length = v.Length;
            float[] u = new float[length];

            for (int i = 0; i < length; i++)
            {
                u[i] = Maths.Byte(v[i]);
            }
            return u;
        }
        /// <summary>
        /// Returns the module of the elements of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToAbs(this float[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Math.Abs(v[i]);
            }
            return H;
        }
        /// <summary>
        /// Negates all elements of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToNegate(this float[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = -v[i];
            }
            return H;
        }
        /// <summary>
        /// Negates all elements of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] ToNegate(this Complex32[] v)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = -v[i];
            }
            return H;
        }
        /// <summary>
        /// Returns a complex vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] ToComplex(this float[] v)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i];
            }
            return H;
        }
        /// <summary>
        /// Returns the module of elements of a complex vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToAbs(this Complex32[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Abs;
            }
            return H;
        }
        /// <summary>
        /// Returns the angle of the elements of a complex vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToAngle(this Complex32[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Angle;
            }
            return H;
        }
        /// <summary>
        /// Returns the real part of the elements of a complex vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToReal(this Complex32[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Real;
            }
            return H;
        }
        /// <summary>
        /// Returns the imaginary part of the elements of a complex vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] ToImag(this Complex32[] v)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Imag;
            }
            return H;
        }
        /// <summary>
        /// Returns a complex conjugate vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] ToConjugate(this Complex32[] v)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Conjugate;
            }
            return H;
        }
        #endregion

        #region Vector statistics
        /// <summary>
        /// Returns the total value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Sum(this float[] v)
        {
            float total = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total += v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the total value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Sum(this Complex32[] v)
        {
            Complex32 total = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total += v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the total product of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Mul(this float[] v)
        {
            float total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total *= v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the total product of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Mul(this Complex32[] v)
        {
            Complex32 total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total *= v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the common quotient of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Div(this float[] v)
        {
            float total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total /= v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the common quotient of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Div(this Complex32[] v)
        {
            Complex32 total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                //if (Maths.IsSingular(v[i])) continue;
                total /= v[i];
            }
            return total;
        }
        /// <summary>
        /// Returns the average value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Mean(this float[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Returns the average value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Mean(this Complex32[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Returns the variance value.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Var(this float[] v)
        {
            int length = v.Length;
            float mean = Matrice.Mean(v);
            float sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += (float)Math.Pow(v[i] - mean, 2);
            }

            return sum / (length - 1);
        }
        /// <summary>
        /// Returns the variance value.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Var(this Complex32[] v)
        {
            int length = v.Length;
            Complex32 mean = Matrice.Mean(v);
            Complex32 sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(v[i] - mean, 2);
            }

            return sum / (length - 1);
        }
        /// <summary>
        /// Returns the variance value.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Value</returns>
        public static float Var(this float[] x, float[] y)
        {
            int length = x.Length;
            float sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += (float)Math.Pow(x[i] - y[i], 2);
            }

            return sum / (length - 1);
        }
        /// <summary>
        /// Returns the variance value.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Var(this Complex32[] x, Complex32[] y)
        {
            int length = x.Length;
            Complex32 sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += Maths.Pow(x[i] - y[i], 2);
            }

            return sum / (length - 1);
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float StnDev(this float[] v)
        {
            return (float)Math.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 StnDev(this Complex32[] v)
        {
            return Maths.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Value</returns>
        public static float StnDev(this float[] x, float[] y)
        {
            return (float)Math.Sqrt(Matrice.Var(x, y));
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 StnDev(this Complex32[] x, Complex32[] y)
        {
            return Maths.Sqrt(Matrice.Var(x, y));
        }
        /// <summary>
        /// Returns the value of the vector mode.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Mode(this float[] v)
        {
            int count = 0;
            int length = v.Length;
            float frequent = 0;
            int i, j, k;

            for (i = 0; i < length; i++)
            {
                k = 1;
                for (j = i + 1; j < length; j++)
                {
                    if (v[i] == v[j]) k++;
                }
                if (k > count)
                {
                    count = k;
                    frequent = v[i];
                }
            }
            return frequent;
        }
        /// <summary>
        /// Returns the value of the vector mode.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Mode(this Complex32[] v)
        {
            int count = 0;
            int length = v.Length;
            Complex32 frequent = 0;
            int i, j, k;

            for (i = 0; i < length; i++)
            {
                k = 1;
                for (j = i + 1; j < length; j++)
                {
                    if (v[i] == v[j]) k++;
                }
                if (k > count)
                {
                    count = k;
                    frequent = v[i];
                }
            }
            return frequent;
        }
        /// <summary>
        /// Sorts the vector.
        /// </summary>
        /// <param name="v">Array</param>
        public static float[] Sort(this float[] v)
        {
            float[] w = (float[])v.Clone();
            Array.Sort(w);
            return w;
        }
        /// <summary>
        /// Gets the value of the minimum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Min(this float[] v)
        {
            return Min(v, out _);
        }
        /// <summary>
        /// Gets the value of the minimum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="index">Max index</param>
        /// <returns>Value</returns>
        public static float Min(this float[] v, out int index)
        {
            int length = v.Length;
            float minimum = float.MaxValue;
            float c;
            index = 0;

            for (int i = 0; i < length; i++)
            {
                c = v[i];

                if (c < minimum)
                {
                    minimum = c;
                    index = i;
                }
            }

            return minimum;
        }
        /// <summary>
        /// Gets the value of the maximum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Max(this float[] v)
        {
            return Max(v, out _);
        }
        /// <summary>
        /// Gets the value of the maximum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="index">Max index</param>
        /// <returns>Value</returns>
        public static float Max(this float[] v, out int index)
        {
            int length = v.Length;
            float maximum = float.MinValue;
            float c;
            index = 0;

            for (int i = 0; i < length; i++)
            {
                c = v[i];

                if (c > maximum)
                {
                    maximum = c;
                    index = i;
                }
            }
            return maximum;
        }
        /// <summary>
        /// Gets the value of the vector element corresponding to the morphology mode.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="mode">Morphology mode</param>
        /// <returns>Value</returns>
        public static float Morph(this float[] v, MorphologyMode mode = MorphologyMode.Median)
        {
            float[] u = (float[])v.Clone();
            Array.Sort(u);
            return u[LinealgOptions.MorphologyFilter.GetFilterRank(mode, u.Length)];
        }
        /// <summary>
        /// Returns the covariance value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Cov(this float[] v)
        {
            int xlength = v.Length;
            float xv = Matrice.Mean(v);
            float total = 0;
            int i;

            for (i = 0; i < xlength; i++)
            {
                total += v[i] * v[i] - xv * xv;
            }
            return total / (xlength - 1);
        }
        /// <summary>
        /// Returns the covariance value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex32 Cov(this Complex32[] v)
        {
            int xlength = v.Length;
            Complex32 xv = Matrice.Mean(v);
            Complex32 total = 0;
            int i;

            for (i = 0; i < xlength; i++)
            {
                total += v[i] * v[i] - xv * xv;
            }
            return total / (xlength - 1);
        }
        /// <summary>
        /// Returns the entropy of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Value</returns>
        public static float Entropy(this float[] v)
        {
            float H = 0;
            int length = v.Length, i;

            for (i = 0; i < length; i++)
            {
                if (v[i] > 0)
                {
                    H += -v[i] * (float)Maths.Log2(v[i]);
                }
            }
            return H;
        }
        /// <summary>
        /// Returns the normalized matrix.
        /// </summary>
        /// <param name="m">Vector</param>
        /// <returns>Vector</returns>
        public static float[] Normalized(this float[] m)
        {
            int ml = m.GetLength(0);
            float[] u = new float[ml];
            var min = m.Min();
            var max = m.Max();
            var dm = max - min;

            for (int i = 0; i < ml; i++)
            {
                u[i] = (m[i] - min) / dm;
            }

            return u;
        }
        #endregion

        #region Vector concatenation
        /// <summary>
        /// Implements vector concatenation.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Concat(this float[] a, float[] b)
        {
            int na = a.Length, nb = b.Length, i;
            float[] v = new float[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        /// <summary>
        /// Implements vector concatenation.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Concat(this Complex32[] a, Complex32[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex32[] v = new Complex32[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        /// <summary>
        /// Implements vector concatenation.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Concat(this Complex32[] a, float[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex32[] v = new Complex32[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        /// <summary>
        /// Implements vector concatenation.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Concat(this float[] a, Complex32[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex32[] v = new Complex32[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        #endregion

        #region Vector as diagonal matrix
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: A * diag(v).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static float[,] Dot(this float[,] m, float[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] temp = new float[r0, r1];
            float alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: A * diag(v).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this Complex32[,] m, Complex32[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: A * diag(v).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this Complex32[,] m, float[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: A * diag(v).
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this float[,] m, Complex32[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }

        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: diag(v) * A.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static float[,] Dot(this float[] v, float[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[,] temp = new float[r0, r1];
            float alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: diag(v) * A.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this Complex32[] v, Complex32[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: diag(v) * A.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this Complex32[] v, float[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        /// <summary>
        /// Implements the scalar product of a matrix by a vector of the form: diag(v) * A.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="v">Array</param>
        /// <param name="inverse">Use inverse to diagonal matrix or not</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dot(this float[] v, Complex32[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[,] temp = new Complex32[r0, r1];
            Complex32 alpha;
            int i, j;

            if (!inverse)
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];
                    for (i = 0; i < r0; i++)
                    {
                        temp[i, j] = m[i, j] * alpha;
                    }
                }
            }
            else
            {
                for (j = 0; j < r1; j++)
                {
                    alpha = v[j];

                    if (alpha != 0)
                    {
                        for (i = 0; i < r0; i++)
                        {
                            temp[i, j] = m[i, j] / alpha;
                        }
                    }
                }
            }

            return temp;
        }
        #endregion

        // Vector special

        #region Vector dot
        /// <summary>
        /// Implements a scalar product of vectors of the form: a * b'.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static float Dot(this float[] a, float[] b)
        {
            int length = a.Length, i;
            float sum = 0;

            for (i = 0; i < length; i++)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a * b'.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Dot(this Complex32[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32 sum = 0;

            for (i = 0; i < length; i++)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a * b'.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Dot(this Complex32[] a, float[] b)
        {
            int length = a.Length, i;
            Complex32 sum = 0;

            for (i = 0; i < length; i++)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a * b'.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Value</returns>
        public static Complex32 Dot(this float[] a, Complex32[] b)
        {
            int length = a.Length, i;
            Complex32 sum = 0;

            for (i = 0; i < length; i++)
            {
                sum += a[i] * b[i];
            }
            return sum;
        }
        /// <summary>
        /// Implements scalar multiplication of a vector by a matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Dot(this float[] v, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            float[] temp = new float[r0];

            Parallel.For(0, r0, i =>
            {
                int j;
                for (j = 0; j < r1; j++)
                {
                    temp[i] += v[j] * m[i, j];
                }
            }
            );

            return temp;
        }
        /// <summary>
        /// Implements scalar multiplication of a vector by a matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Dot(this float[] v, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[] temp = new Complex32[r0];

            Parallel.For(0, r0, i =>
            {
                int j;
                for (j = 0; j < r1; j++)
                {
                    temp[i] += v[j] * m[i, j];
                }
            }
            );

            return temp;
        }
        /// <summary>
        /// Implements scalar multiplication of a vector by a matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Dot(this Complex32[] v, Complex32[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[] temp = new Complex32[r0];

            Parallel.For(0, r0, i =>
            {
                int j;
                for (j = 0; j < r1; j++)
                {
                    temp[i] += v[j] * m[i, j];
                }
            }
            );

            return temp;
        }
        /// <summary>
        /// Implements scalar multiplication of a vector by a matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Dot(this Complex32[] v, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex32[] temp = new Complex32[r0];

            Parallel.For(0, r0, i =>
            {
                int j;
                for (j = 0; j < r1; j++)
                {
                    temp[i] += v[j] * m[i, j];
                }
            }
            );

            return temp;
        }
        #endregion

        #region Vector/matrix multiply
        /// <summary>
        /// Implements a scalar product of vectors of the form: a' * b, 
        /// where ' is the transpose sign.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[,] Dotp(this float[] a, float[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            float[,] H = new float[l0, l1];
            float c;
            int i, j;

            for (j = 0; j < l0; j++)
            {
                c = a[j];
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += c * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a' * b, 
        /// where ' is the transpose sign.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dotp(this Complex32[] a, Complex32[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex32[,] H = new Complex32[l0, l1];
            Complex32 c;
            int i, j;

            for (j = 0; j < l0; j++)
            {
                c = a[j];
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += c * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a' * b, 
        /// where ' is the transpose sign.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dotp(this Complex32[] a, float[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex32[,] H = new Complex32[l0, l1];
            Complex32 c;
            int i, j;

            for (j = 0; j < l0; j++)
            {
                c = a[j];
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += c * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements a scalar product of vectors of the form: a' * b, 
        /// where ' is the transpose sign.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[,] Dotp(this float[] a, Complex32[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex32[,] H = new Complex32[l0, l1];
            Complex32 c;
            int i, j;

            for (j = 0; j < l0; j++)
            {
                c = a[j];
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += c * b[i];
                }
            }
            return H;
        }
        #endregion

        #region Vector convolutions
        /// <summary>
        /// Implements discrete convolution of vectors.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="u">Array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Array</returns>
        public static float[] Conv(this float[] v, float[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            float sum, div;
            float[] uv = new float[n];

            if (normalize)
            {
                // normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0; div = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                            div += u[j];
                        }
                    }

                    uv[i] = sum / div;
                }
            }
            else
            {
                // non-normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                        }
                    }

                    uv[i] = sum;
                }
            }

            return uv;
        }
        /// <summary>
        /// Implements discrete convolution of vectors.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="u">Array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Array</returns>
        public static Complex32[] Conv(this Complex32[] v, Complex32[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex32 sum, div;
            Complex32[] uv = new Complex32[n];

            if (normalize)
            {
                // normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0; div = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                            div += u[j];
                        }
                    }

                    uv[i] = sum / div;
                }
            }
            else
            {
                // non-normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                        }
                    }

                    uv[i] = sum;
                }
            }

            return uv;
        }
        /// <summary>
        /// Implements discrete convolution of vectors.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="u">Array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Array</returns>
        public static Complex32[] Conv(this Complex32[] v, float[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex32 sum, div;
            Complex32[] uv = new Complex32[n];

            if (normalize)
            {
                // normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0; div = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                            div += u[j];
                        }
                    }

                    uv[i] = sum / div;
                }
            }
            else
            {
                // non-normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                        }
                    }

                    uv[i] = sum;
                }
            }

            return uv;
        }
        /// <summary>
        /// Implements discrete convolution of vectors.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="u">Array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Array</returns>
        public static Complex32[] Conv(this float[] v, Complex32[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex32 sum, div;
            Complex32[] uv = new Complex32[n];

            if (normalize)
            {
                // normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0; div = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                            div += u[j];
                        }
                    }

                    uv[i] = sum / div;
                }
            }
            else
            {
                // non-normalize convolution:
                for (i = 0; i < n; i++)
                {
                    sum = 0;
                    k = i - r2;

                    for (j = 0; j < r; j++)
                    {
                        c = k + j;

                        if (c >= 0 && c < n)
                        {
                            sum += v[c] * u[j];
                        }
                    }

                    uv[i] = sum;
                }
            }

            return uv;
        }
        #endregion

        #region Vector morphology
        /// <summary>
        /// Returns the vector result of morphology minimum.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Min(this float[] v, int r)
        {
            return Morph(v, r, MorphologyMode.Erosion);
        }
        /// <summary>
        /// Returns the vector result of morphology maximum.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Max(this float[] v, int r)
        {
            return Morph(v, r, MorphologyMode.Dilatation);
        }
        /// <summary>
        /// Returns the vector result of morphology.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        /// <param name="mode">Morphology mode</param>
        /// <returns>Array</returns>
        public static float[] Morph(this float[] v, int r, MorphologyMode mode = MorphologyMode.Median)
        {
            return LinealgOptions.MorphologySortFilter.Apply(v, r / 2, mode);
        }
        #endregion

        #region Vector mean
        /// <summary>
        /// Returns the result vector of local averaging.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Mean(this float[] v, int r)
        {
            return LinealgOptions.MeanFilter.Mean(v, r);
        }
        /// <summary>
        /// Returns the result vector of local averaging.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static Complex32[] Mean(this Complex32[] v, int r)
        {
            return LinealgOptions.MeanFilter.Mean(v, r);
        }
        /// <summary>
        /// Returns the result vector of local weighted averaging.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="w">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Mean(this float[] v, float[] w, int r)
        {
            return LinealgOptions.MeanFilter.MeanWeighted(v, w, r);
        }
        /// <summary>
        /// Returns the result vector of local weighted averaging.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="w">Array</param>
        /// <param name="r">Radius</param>
        public static Complex32[] Mean(this Complex32[] v, Complex32[] w, int r)
        {
            return LinealgOptions.MeanFilter.MeanWeighted(v, w, r);
        }
        #endregion

        // Imaging voids

        #region Rotate voids
        /// <summary>
        /// Rotates matrix by rotation value.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="rotation">Rotation</param>
        /// <returns>Matrix</returns>
        public static float[,] Rotate(this float[,] matrix, RotationMode rotation)
        {
            return rotation switch
            {
                RotationMode.R0 => matrix,
                RotationMode.R90 => Rotate90(matrix),
                RotationMode.R180 => Rotate180(matrix),
                RotationMode.R270 => Rotate270(matrix),
                _ => matrix,
            };
        }

        #region Private rotation
        /// <summary>
        /// Rotates the matrix by 90 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static float[,] Rotate90(float[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            float[,] H = new float[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 180 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static float[,] Rotate180(float[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            float[,] H = new float[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 270 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static float[,] Rotate270(float[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            float[,] H = new float[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[i, w - j - 1];
                }
            }

            return H;
        }
        #endregion

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static float[,] Rotate(this float[,] matrix, float angle, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            return Rotate(matrix, angle, 0, interpolationMode);
        }
        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static float[,] Rotate(this float[,] matrix, float angle, float value, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return RotateBicubic(matrix, angle, value);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return RotateBilinear(matrix, angle, value);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return RotateNearestNeighbor(matrix, angle, value);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private rotate

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static float[,] RotateNearestNeighbor(this float[,] matrix, float angle, float value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy;
            // source pixel's coordinates
            int ox, oy;
            // output
            float[,] H = new float[newHeight, newWidth];

            // check pixel format
            // ARGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;
                for (int x = 0; x < newWidth; x++)
                {
                    // coordinate of the nearest point
                    ox = (int)(angleCos * cx + angleSin * cy + oldXradius);
                    oy = (int)(-angleSin * cx + angleCos * cy + oldYradius);

                    // validate source pixel's coordinates
                    if ((ox < 0) || (oy < 0) || (ox >= width) || (oy >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        // fill destination image with pixel from source image
                        H[y, x] = matrix[oy, ox];
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static float[,] RotateBilinear(this float[,] matrix, float angle, float value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy;
            // coordinates of source points
            double ox, oy, tx, ty, dx1, dy1, dx2, dy2;
            int ox1, oy1, ox2, oy2;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // output
            float[,] H = new float[newHeight, newWidth];

            // RGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                // do some pre-calculations of source points' coordinates
                // (calculate the part which depends on y-loop, but does not
                // depend on x-loop)
                tx = angleSin * cy + oldXradius;
                ty = angleCos * cy + oldYradius;

                cx = -newXradius;
                for (int x = 0; x < newWidth; x++)
                {
                    // coordinates of source point
                    ox = tx + angleCos * cx;
                    oy = ty - angleSin * cx;

                    // top-left coordinate
                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        // bottom-right coordinate
                        ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                        oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;

                        if ((dx1 = ox - (float)ox1) < 0)
                            dx1 = 0;
                        dx2 = 1.0f - dx1;

                        if ((dy1 = oy - (float)oy1) < 0)
                            dy1 = 0;
                        dy2 = 1.0f - dy1;

                        // get four points
                        var p1 = matrix[oy1, ox1];
                        var p2 = matrix[oy1, ox2];
                        var p3 = matrix[oy2, ox1];
                        var p4 = matrix[oy2, ox2];

                        // interpolate using 4 points
                        H[y, x] = (float)(
                            dy2 * (dx2 * p1 + dx1 * p2) +
                            dy1 * (dx2 * p3 + dx1 * p4));
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static float[,] RotateBicubic(this float[,] matrix, float angle, float value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            float oldXradius = (float)(width - 1) / 2;
            float oldYradius = (float)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            float newXradius = (float)(newWidth - 1) / 2;
            float newYradius = (float)(newHeight - 1) / 2;

            // angle's sine and cosine
            float angleRad = -angle * Maths.Pi / 180.0f;
            float angleCos = Maths.Cos(angleRad);
            float angleSin = Maths.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            float cx, cy;
            // coordinates of source points and cooefficiens
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            // destination pixel values
            float g;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // output
            float[,] H = new float[newHeight, newWidth];

            // grayscale
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;

                for (int x = 0; x < newWidth; x++)
                {
                    // coordinates of source point
                    ox = angleCos * cx + angleSin * cy + oldXradius;
                    oy = -angleSin * cx + angleCos * cy + oldYradius;

                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        dx = ox - ox1;
                        dy = oy - oy1;

                        // initial pixel value
                        g = 0;

                        for (int n = -1; n < 3; n++)
                        {
                            // get Y coefficient
                            k1 = Kernel.Bicubic((float)(dy - n));

                            oy2 = oy1 + n;
                            if (oy2 < 0)
                                oy2 = 0;
                            if (oy2 > ymax)
                                oy2 = ymax;

                            for (int m = -1; m < 3; m++)
                            {
                                // get X coefficient
                                k2 = k1 * Kernel.Bicubic((float)(m - dx));

                                ox2 = ox1 + m;
                                if (ox2 < 0)
                                    ox2 = 0;
                                if (ox2 > xmax)
                                    ox2 = xmax;

                                g += k2 * matrix[oy2, ox2];
                            }
                        }
                        H[y, x] = g;
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        #endregion

        /// <summary>
        /// Rotates matrix by rotation value.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="rotation">Rotation</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Rotate(this Complex32[,] matrix, RotationMode rotation)
        {
            return rotation switch
            {
                RotationMode.R0 => matrix,
                RotationMode.R90 => Rotate90(matrix),
                RotationMode.R180 => Rotate180(matrix),
                RotationMode.R270 => Rotate270(matrix),
                _ => matrix,
            };
        }

        #region Private rotation
        /// <summary>
        /// Rotates the matrix by 90 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] Rotate90(Complex32[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            Complex32[,] H = new Complex32[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 180 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] Rotate180(Complex32[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            Complex32[,] H = new Complex32[h, w];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[i, j] = input[h - i - 1, w - j - 1];
                }
            }

            return H;
        }
        /// <summary>
        /// Rotates the depth by 270 degrees.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] Rotate270(Complex32[,] input)
        {
            int h = input.GetLength(0);
            int w = input.GetLength(1);

            Complex32[,] H = new Complex32[w, h];

            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    H[j, i] = input[i, w - j - 1];
                }
            }

            return H;
        }
        #endregion

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Rotate(this Complex32[,] matrix, float angle, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            return Rotate(matrix, angle, 0, interpolationMode);
        }
        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Rotate(this Complex32[,] matrix, float angle, Complex32 value, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return RotateBicubic(matrix, angle, value);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return RotateBilinear(matrix, angle, value);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return RotateNearestNeighbor(matrix, angle, value);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private rotate

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] RotateNearestNeighbor(this Complex32[,] matrix, float angle, Complex32 value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy;
            // source pixel's coordinates
            int ox, oy;
            // output
            Complex32[,] H = new Complex32[newHeight, newWidth];

            // check pixel format
            // ARGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;
                for (int x = 0; x < newWidth; x++)
                {
                    // coordinate of the nearest point
                    ox = (int)(angleCos * cx + angleSin * cy + oldXradius);
                    oy = (int)(-angleSin * cx + angleCos * cy + oldYradius);

                    // validate source pixel's coordinates
                    if ((ox < 0) || (oy < 0) || (ox >= width) || (oy >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        // fill destination image with pixel from source image
                        H[y, x] = matrix[oy, ox];
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] RotateBilinear(this Complex32[,] matrix, float angle, Complex32 value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            double oldXradius = (double)(width - 1) / 2;
            double oldYradius = (double)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            double newXradius = (double)(newWidth - 1) / 2;
            double newYradius = (double)(newHeight - 1) / 2;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy;
            // coordinates of source points
            double ox, oy, tx, ty, dx1, dy1, dx2, dy2;
            int ox1, oy1, ox2, oy2;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // output
            Complex32[,] H = new Complex32[newHeight, newWidth];

            // RGB
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                // do some pre-calculations of source points' coordinates
                // (calculate the part which depends on y-loop, but does not
                // depend on x-loop)
                tx = angleSin * cy + oldXradius;
                ty = angleCos * cy + oldYradius;

                cx = -newXradius;
                for (int x = 0; x < newWidth; x++)
                {
                    // coordinates of source point
                    ox = tx + angleCos * cx;
                    oy = ty - angleSin * cx;

                    // top-left coordinate
                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        // bottom-right coordinate
                        ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                        oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;

                        if ((dx1 = ox - (float)ox1) < 0)
                            dx1 = 0;
                        dx2 = 1.0f - dx1;

                        if ((dy1 = oy - (float)oy1) < 0)
                            dy1 = 0;
                        dy2 = 1.0f - dy1;

                        // get four points
                        var p1 = matrix[oy1, ox1];
                        var p2 = matrix[oy1, ox2];
                        var p3 = matrix[oy2, ox1];
                        var p4 = matrix[oy2, ox2];

                        // interpolate using 4 points
                        H[y, x] = (Complex32)(
                            dy2 * (dx2 * p1 + dx1 * p2) +
                            dy1 * (dx2 * p3 + dx1 * p4));
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        /// <summary>
        /// Rotates matrix by angle.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="angle">Angle</param>
        /// <param name="value">Value</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] RotateBicubic(this Complex32[,] matrix, float angle, Complex32 value)
        {
            // get source image size
            int width = matrix.GetLength(1);
            int height = matrix.GetLength(0);
            float oldXradius = (float)(width - 1) / 2;
            float oldYradius = (float)(height - 1) / 2;

            // get destination image size
            int newWidth = width;
            int newHeight = height;
            float newXradius = (float)(newWidth - 1) / 2;
            float newYradius = (float)(newHeight - 1) / 2;

            // angle's sine and cosine
            float angleRad = -angle * Maths.Pi / 180.0f;
            float angleCos = Maths.Cos(angleRad);
            float angleSin = Maths.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            float cx, cy;
            // coordinates of source points and cooefficiens
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            // destination pixel values
            Complex32 g;
            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;
            // output
            Complex32[,] H = new Complex32[newHeight, newWidth];

            // grayscale
            cy = -newYradius;
            for (int y = 0; y < newHeight; y++)
            {
                cx = -newXradius;

                for (int x = 0; x < newWidth; x++)
                {
                    // coordinates of source point
                    ox = angleCos * cx + angleSin * cy + oldXradius;
                    oy = -angleSin * cx + angleCos * cy + oldYradius;

                    ox1 = (int)ox;
                    oy1 = (int)oy;

                    // validate source pixel's coordinates
                    if ((ox1 < 0) || (oy1 < 0) || (ox1 >= width) || (oy1 >= height))
                    {
                        // fill destination image with filler
                        H[y, x] = value;
                    }
                    else
                    {
                        dx = ox - ox1;
                        dy = oy - oy1;

                        // initial pixel value
                        g = 0;

                        for (int n = -1; n < 3; n++)
                        {
                            // get Y coefficient
                            k1 = Kernel.Bicubic((float)(dy - n));

                            oy2 = oy1 + n;
                            if (oy2 < 0)
                                oy2 = 0;
                            if (oy2 > ymax)
                                oy2 = ymax;

                            for (int m = -1; m < 3; m++)
                            {
                                // get X coefficient
                                k2 = k1 * Kernel.Bicubic((float)(m - dx));

                                ox2 = ox1 + m;
                                if (ox2 < 0)
                                    ox2 = 0;
                                if (ox2 > xmax)
                                    ox2 = xmax;

                                g += k2 * matrix[oy2, ox2];
                            }
                        }
                        H[y, x] = g;
                    }
                    cx++;
                }
                cy++;
            }

            return H;
        }

        #endregion
        #endregion

        #region Resize voids
        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static float[,] Resize(this float[,] input, int h, int w, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return ResizeBicubic(input, h, w);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return ResizeBilinear(input, h, w);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return ResizeNearestNeighbor(input, h, w);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private methods

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static float[,] ResizeBicubic(this float[,] input, int h, int w)
        {
            // get source size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // coordinates of source points and cooefficiens
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            float g;

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;

            // output
            float[,] H = new float[h, w];

            // grayscale
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                for (int x = 0; x < w; x++)
                {
                    // X coordinates
                    ox = x * xFactor - 0.5f;
                    ox1 = (int)ox;
                    dx = ox - ox1;

                    // initial pixel value
                    g = 0;

                    for (int n = -1; n < 3; n++)
                    {
                        // get Y coefficient
                        k1 = Kernel.Bicubic((float)(dy - n));

                        oy2 = oy1 + n;
                        if (oy2 < 0)
                            oy2 = 0;
                        if (oy2 > ymax)
                            oy2 = ymax;

                        for (int m = -1; m < 3; m++)
                        {
                            // get X coefficient
                            k2 = k1 * Kernel.Bicubic((float)(m - dx));

                            ox2 = ox1 + m;
                            if (ox2 < 0)
                                ox2 = 0;
                            if (ox2 > xmax)
                                ox2 = xmax;

                            g += k2 * input[oy2, ox2];
                        }
                    }

                    H[y, x] = g;
                }
            }

            return H;
        }

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static float[,] ResizeBilinear(this float[,] input, int h, int w)
        {
            // get source image size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;

            // output
            float[,] H = new float[h, w];

            // for each line
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                double oy = (double)y * yFactor;
                int oy1 = (int)oy;
                int oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;
                double dy1 = oy - (double)oy1;
                double dy2 = 1.0 - dy1;

                // for each pixel
                for (int x = 0; x < w; x++)
                {
                    // X coordinates
                    double ox = (double)x * xFactor;
                    int ox1 = (int)ox;
                    int ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                    double dx1 = ox - (double)ox1;
                    double dx2 = 1.0 - dx1;

                    // get four points
                    var p1 = input[oy1, ox1];
                    var p2 = input[oy1, ox2];
                    var p3 = input[oy2, ox1];
                    var p4 = input[oy2, ox2];

                    // interpolate using 4 points
                    H[y, x] = (float)(
                            dy2 * (dx2 * (p1) + dx1 * (p2)) +
                            dy1 * (dx2 * (p3) + dx1 * (p4)));
                }
            }

            return H;
        }

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static float[,] ResizeNearestNeighbor(this float[,] input, int h, int w)
        {
            // get source image size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // output
            float[,] H = new float[h, w];

            // for each line
            for (int y = 0; y < h; y++)
            {
                var oy = (int)(y * yFactor);

                // for each pixel
                for (int x = 0; x < w; x++)
                {
                    var ox = (int)(x * xFactor);

                    H[y, x] = input[oy, ox];
                }
            }

            return H;
        }

        #endregion

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Resize(this Complex32[,] input, int h, int w, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return ResizeBicubic(input, h, w);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return ResizeBilinear(input, h, w);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return ResizeNearestNeighbor(input, h, w);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private methods

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] ResizeBicubic(this Complex32[,] input, int h, int w)
        {
            // get source size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // coordinates of source points and cooefficiens
            float ox, oy, dx, dy, k1, k2;
            int ox1, oy1, ox2, oy2;
            Complex32 g;

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;

            // output
            Complex32[,] H = new Complex32[h, w];

            // grayscale
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                for (int x = 0; x < w; x++)
                {
                    // X coordinates
                    ox = x * xFactor - 0.5f;
                    ox1 = (int)ox;
                    dx = ox - ox1;

                    // initial pixel value
                    g = 0;

                    for (int n = -1; n < 3; n++)
                    {
                        // get Y coefficient
                        k1 = Kernel.Bicubic((float)(dy - n));

                        oy2 = oy1 + n;
                        if (oy2 < 0)
                            oy2 = 0;
                        if (oy2 > ymax)
                            oy2 = ymax;

                        for (int m = -1; m < 3; m++)
                        {
                            // get X coefficient
                            k2 = k1 * Kernel.Bicubic((float)(m - dx));

                            ox2 = ox1 + m;
                            if (ox2 < 0)
                                ox2 = 0;
                            if (ox2 > xmax)
                                ox2 = xmax;

                            g += k2 * input[oy2, ox2];
                        }
                    }

                    H[y, x] = g;
                }
            }

            return H;
        }

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] ResizeBilinear(this Complex32[,] input, int h, int w)
        {
            // get source image size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // width and height decreased by 1
            int ymax = height - 1;
            int xmax = width - 1;

            // output
            Complex32[,] H = new Complex32[h, w];

            // for each line
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                double oy = (double)y * yFactor;
                int oy1 = (int)oy;
                int oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;
                double dy1 = oy - (double)oy1;
                double dy2 = 1.0 - dy1;

                // for each pixel
                for (int x = 0; x < w; x++)
                {
                    // X coordinates
                    double ox = (double)x * xFactor;
                    int ox1 = (int)ox;
                    int ox2 = (ox1 == xmax) ? ox1 : ox1 + 1;
                    double dx1 = ox - (double)ox1;
                    double dx2 = 1.0 - dx1;

                    // get four points
                    var p1 = input[oy1, ox1];
                    var p2 = input[oy1, ox2];
                    var p3 = input[oy2, ox1];
                    var p4 = input[oy2, ox2];

                    // interpolate using 4 points
                    H[y, x] = (Complex32)(
                        dy2 * (dx2 * (p1) + dx1 * (p2)) +
                        dy1 * (dx2 * (p3) + dx1 * (p4)));
                }
            }

            return H;
        }

        /// <summary>
        /// Returns resized matrix.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] ResizeNearestNeighbor(this Complex32[,] input, int h, int w)
        {
            // get source image size
            int width = input.GetLength(1);
            int height = input.GetLength(0);

            float xFactor = (float)width / w;
            float yFactor = (float)height / h;

            // output
            Complex32[,] H = new Complex32[h, w];

            // for each line
            for (int y = 0; y < h; y++)
            {
                var oy = (int)(y * yFactor);

                // for each pixel
                for (int x = 0; x < w; x++)
                {
                    var ox = (int)(x * xFactor);

                    H[y, x] = input[oy, ox];
                }
            }

            return H;
        }

        #endregion

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Array</returns>
        public static float[] Resize(this float[] input, int h, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return ResizeBicubic(input, h);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return ResizeBilinear(input, h);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return ResizeNearestNeighbor(input, h);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private methods

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static float[] ResizeBicubic(this float[] input, int h)
        {
            // get source size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // coordinates of source points and cooefficiens
            float oy, dy, k1;
            int oy1, oy2;
            float g;

            // width and height decreased by 1
            int ymax = height - 1;

            // output
            float[] H = new float[h];

            // grayscale
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                // initial pixel value
                g = 0;

                for (int n = -1; n < 3; n++)
                {
                    // get Y coefficient
                    k1 = Kernel.Bicubic((float)(dy - n));

                    oy2 = oy1 + n;
                    if (oy2 < 0)
                        oy2 = 0;
                    if (oy2 > ymax)
                        oy2 = ymax;

                    g += k1 * input[oy2];
                }

                H[y] = g;
            }

            return H;
        }

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static float[] ResizeBilinear(this float[] input, int h)
        {
            // get source image size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // width and height decreased by 1
            int ymax = height - 1;

            // output
            float[] H = new float[h];

            // for each line
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                double oy = (double)y * yFactor;
                int oy1 = (int)oy;
                int oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;
                double dy1 = oy - (double)oy1;
                double dy2 = 1.0 - dy1;

                // get 2 points
                var p1 = input[oy1];
                var p3 = input[oy2];

                // interpolate using 2 points
                H[y] = (float)(
                        dy2 * p1 +
                        dy1 * p3);
            }

            return H;
        }

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static float[] ResizeNearestNeighbor(this float[] input, int h)
        {
            // get source image size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // output
            float[] H = new float[h];

            // for each line
            for (int y = 0; y < h; y++)
            {
                var oy = (int)(y * yFactor);

                // for each pixel
                H[y] = input[oy];
            }

            return H;
        }

        #endregion

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Array</returns>
        public static Complex32[] Resize(this Complex32[] input, int h, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            if (interpolationMode == InterpolationMode.Bicubic)
            {
                return ResizeBicubic(input, h);
            }
            else if (interpolationMode == InterpolationMode.Bilinear)
            {
                return ResizeBilinear(input, h);
            }
            else if (interpolationMode == InterpolationMode.NearestNeighbor)
            {
                return ResizeNearestNeighbor(input, h);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        #region Private methods

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static Complex32[] ResizeBicubic(this Complex32[] input, int h)
        {
            // get source size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // coordinates of source points and cooefficiens
            float oy, dy, k1;
            int oy1, oy2;
            Complex32 g;

            // width and height decreased by 1
            int ymax = height - 1;

            // output
            Complex32[] H = new Complex32[h];

            // grayscale
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                oy = y * yFactor - 0.5f;
                oy1 = (int)oy;
                dy = oy - oy1;

                // initial pixel value
                g = 0;

                for (int n = -1; n < 3; n++)
                {
                    // get Y coefficient
                    k1 = Kernel.Bicubic((float)(dy - n));

                    oy2 = oy1 + n;
                    if (oy2 < 0)
                        oy2 = 0;
                    if (oy2 > ymax)
                        oy2 = ymax;

                    g += k1 * input[oy2];
                }

                H[y] = g;
            }

            return H;
        }

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static Complex32[] ResizeBilinear(this Complex32[] input, int h)
        {
            // get source image size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // width and height decreased by 1
            int ymax = height - 1;

            // output
            Complex32[] H = new Complex32[h];

            // for each line
            for (int y = 0; y < h; y++)
            {
                // Y coordinates
                double oy = (double)y * yFactor;
                int oy1 = (int)oy;
                int oy2 = (oy1 == ymax) ? oy1 : oy1 + 1;
                double dy1 = oy - (double)oy1;
                double dy2 = 1.0 - dy1;

                // get 2 points
                var p1 = input[oy1];
                var p3 = input[oy2];

                // interpolate using 2 points
                H[y] = (Complex32)(
                        dy2 * p1 +
                        dy1 * p3);
            }

            return H;
        }

        /// <summary>
        /// Returns resized vector.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="h">Length</param>
        /// <returns>Array</returns>
        private static Complex32[] ResizeNearestNeighbor(this Complex32[] input, int h)
        {
            // get source image size
            int height = input.GetLength(0);

            float yFactor = (float)height / h;

            // output
            Complex32[] H = new Complex32[h];

            // for each line
            for (int y = 0; y < h; y++)
            {
                var oy = (int)(y * yFactor);

                // for each pixel
                H[y] = input[oy];
            }

            return H;
        }

        #endregion

        // preserving proportions

        /// <summary>
        /// Resize method with preserving proportions.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="value">Background value</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static float[,] ResizePreserved(this float[,] input, int h, int w, float value, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            int width = input.GetLength(1);
            int height = input.GetLength(0);
            int max = Math.Max(width, height);
            var rect = new Rectangle((max - width) / 2, (max - height) / 2, width, height);
            var temp = new float[max, max].Add(value);

            for (int y = 0; y < rect.Height; y++)
            {
                for (int x = 0; x < rect.Width; x++)
                {
                    temp[y + rect.Y, x + rect.X] = input[y, x];
                }
            }

            return temp.Resize(h, w, interpolationMode);
        }

        /// <summary>
        /// Resize method with preserving proportions.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static float[,] ResizePreserved(this float[,] input, int h, int w, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            int width = w;
            int height = h;
            int max = Math.Max(width, height);
            var rect = new Rectangle((max - width) / 2, (max - height) / 2, width, height);
            var resized = input.Resize(max, max, interpolationMode);
            var temp = new float[rect.Height, rect.Width];

            for (int y = 0; y < rect.Height; y++)
            {
                for (int x = 0; x < rect.Width; x++)
                {
                    temp[y, x] = resized[y + rect.Y, x + rect.X];
                }
            }

            return temp;
        }

        /// <summary>
        /// Resize method with preserving proportions.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="value">Background value</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] ResizePreserved(this Complex32[,] input, int h, int w, Complex32 value, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            int width = input.GetLength(1);
            int height = input.GetLength(0);
            int max = Math.Max(width, height);
            var rect = new Rectangle((max - width) / 2, (max - height) / 2, width, height);
            var temp = new Complex32[max, max].Add(value);

            for (int y = 0; y < rect.Height; y++)
            {
                for (int x = 0; x < rect.Width; x++)
                {
                    temp[y + rect.Y, x + rect.X] = input[y, x];
                }
            }

            return temp.Resize(h, w, interpolationMode);
        }

        /// <summary>
        /// Resize method with preserving proportions.
        /// </summary>
        /// <param name="input">Matrix</param>
        /// <param name="h">Height</param>
        /// <param name="w">Width</param>
        /// <param name="interpolationMode">Interpolation mode</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] ResizePreserved(this Complex32[,] input, int h, int w, InterpolationMode interpolationMode = InterpolationMode.Bicubic)
        {
            int width = w;
            int height = h;
            int max = Math.Max(width, height);
            var rect = new Rectangle((max - width) / 2, (max - height) / 2, width, height);
            var resized = input.Resize(max, max, interpolationMode);
            var temp = new Complex32[rect.Height, rect.Width];

            for (int y = 0; y < rect.Height; y++)
            {
                for (int x = 0; x < rect.Width; x++)
                {
                    temp[y, x] = resized[y + rect.Y, x + rect.X];
                }
            }

            return temp;
        }

        #endregion

        #region Shift voids
        /// <summary>
        /// Implements a shift of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="m">The number of positions to which a shift in height occurs</param>
        /// <param name="l">The number of positions by which the shift occurs in width</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Shift(this Complex32[,] a, int m, int l)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            Complex32[,] temp = new Complex32[l0, l1];
            int i, j;

            for (i = 0; i < l0; i++)
            {
                for (j = 0; j < l1; j++)
                {
                    temp[i, j] = a[Maths.Mod(i - m, l1), Maths.Mod(j - l, l0)];
                }
            }
            return temp;
        }
        /// <summary>
        /// Implements a shift of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="m">The number of positions to which a shift in height occurs</param>
        /// <param name="l">The number of positions by which the shift occurs in width</param>
        /// <returns>Matrix</returns>
        public static float[,] Shift(this float[,] a, int m, int l)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            float[,] temp = new float[l0, l1];
            int i, j;

            for (i = 0; i < l0; i++)
            {
                for (j = 0; j < l1; j++)
                {
                    temp[i, j] = a[Maths.Mod(i - m, l1), Maths.Mod(j - l, l0)];
                }
            }
            return temp;
        }
        /// <summary>
        /// Implements a shift of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="l">Number of positions to shift</param>
        /// <returns>Array</returns>
        public static Complex32[] Shift(this Complex32[] v, int l)
        {
            int N = v.Length;
            Complex32[] temp = new Complex32[N];

            for (int i = 0; i < N; i++)
            {
                temp[i] = v[Maths.Mod(i - l, N)];
            }

            return temp;
        }
        /// <summary>
        /// Implements a shift of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="l">Number of positions to shift</param>
        /// <returns>Array</returns>
        public static float[] Shift(this float[] v, int l)
        {
            int N = v.Length;
            float[] temp = new float[N];

            for (int i = 0; i < N; i++)
            {
                temp[i] = v[Maths.Mod(i - l, N)];
            }

            return temp;
        }
        #endregion

        #region Flip voids
        /// <summary>
        /// Flips matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <returns>Matrix</returns>
        public static float[,] Flip(this float[,] m, Direction direction)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[,] H = new float[ml, mr];
            int i, j;

            // horizontal flipping:
            if (direction == Direction.Horizontal)
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[i, mr - j - 1];
                    }
                }
            }
            // vertical flipping:
            else if (direction == Direction.Vertical)
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[ml - i - 1, j];
                    }
                }
            }
            // both flipping:
            else
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[ml - i - 1, mr - j - 1];
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Flips matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Flip(this Complex32[,] m, Direction direction)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[,] H = new Complex32[ml, mr];
            int i, j;

            // horizontal flipping:
            if (direction == Direction.Horizontal)
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[i, mr - j - 1];
                    }
                }
            }
            // vertical flipping:
            else if (direction == Direction.Vertical)
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[ml - i - 1, j];
                    }
                }
            }
            // both flipping:
            else
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = 0; j < mr; j++)
                    {
                        H[i, j] = m[ml - i - 1, mr - j - 1];
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Flips vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] Flip(this float[] v)
        {
            int mr = v.Length;
            float[] H = new float[mr];

            for (int j = 0; j < mr; j++)
            {
                H[j] = v[mr - j - 1];
            }

            return H;
        }
        /// <summary>
        /// Flips vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Flip(this Complex32[] v)
        {
            int mr = v.Length;
            Complex32[] H = new Complex32[mr];

            for (int j = 0; j < mr; j++)
            {
                H[j] = v[mr - j - 1];
            }

            return H;
        }
        #endregion

        #region Crop voids
        /// <summary>
        /// Returns the specified part of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="start">Starting position</param>
        /// <param name="length">Vector length</param>
        /// <param name="clamp">Clamp crop or not</param>
        /// <returns>Array</returns>
        public static float[] Crop(this float[] a, int start, int length, bool clamp = true)
        {
            // vector param
            int _length = a.GetLength(0);

            // range processing
            int y = clamp ? Maths.Range(start, 0, a.Length) : start;
            int h = clamp ? Maths.Range(length, 0, a.Length - start) : length;

            float[] v = new float[h];

            for (int i = 0; i < h; i++)
            {
                var iy = i + y;

                if (iy < 0 || iy >= _length) continue;

                v[i] = a[iy];
            }

            return v;
        }
        /// <summary>
        /// Returns the specified part of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="start">Starting position</param>
        /// <param name="length">Vector length</param>
        /// <param name="clamp">Clamp crop or not</param>
        /// <returns>Array</returns>
        public static Complex32[] Crop(this Complex32[] a, int start, int length, bool clamp = true)
        {
            // vector param
            int _length = a.GetLength(0);

            // range processing
            int y = clamp ? Maths.Range(start, 0, a.Length) : start;
            int h = clamp ? Maths.Range(length, 0, a.Length - start) : length;

            Complex32[] v = new Complex32[h];

            for (int i = 0; i < h; i++)
            {
                var iy = i + y;

                if (iy < 0 || iy >= _length) continue;

                v[i] = a[iy];
            }

            return v;
        }
        /// <summary>
        /// Crops the matrix to the specified size.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="y">Starting position in height</param>
        /// <param name="x">Starting position in width</param>
        /// <param name="height">Height</param>
        /// <param name="width">Width</param>
        /// <param name="clamp">Clamp crop or not</param>
        /// <returns>Matrix</returns>
        public static float[,] Crop(this float[,] m, int y, int x, int height, int width, bool clamp = true)
        {
            // image params
            int _width = m.GetLength(1);
            int _height = m.GetLength(0);

            // range processing
            int xx = clamp ? Maths.Range(x, 0, m.GetLength(1)) : x;
            int yy = clamp ? Maths.Range(y, 0, m.GetLength(0)) : y;
            int ww = clamp ? Maths.Range(width, 0, m.GetLength(1) - xx) : width;
            int hh = clamp ? Maths.Range(height, 0, m.GetLength(0) - yy) : height;

            float[,] array = new float[hh, ww];

            for (int i = 0; i < hh; i++)
            {
                var iyy = i + yy;

                if (iyy < 0 || iyy >= _height) continue;

                for (int j = 0; j < ww; j++)
                {
                    var jxx = j + xx;

                    if (jxx < 0 || jxx >= _width) continue;

                    array[i, j] = m[iyy, jxx];
                }
            }

            return array;
        }
        /// <summary>
        /// Crops the matrix to the specified size.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="y">Starting position in height</param>
        /// <param name="x">Starting position in width</param>
        /// <param name="height">Height</param>
        /// <param name="width">Width</param>
        /// <param name="clamp">Clamp crop or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Crop(this Complex32[,] m, int y, int x, int height, int width, bool clamp = true)
        {
            // image params
            int _width = m.GetLength(1);
            int _height = m.GetLength(0);

            // range processing
            int xx = clamp ? Maths.Range(x, 0, m.GetLength(1)) : x;
            int yy = clamp ? Maths.Range(y, 0, m.GetLength(0)) : y;
            int ww = clamp ? Maths.Range(width, 0, m.GetLength(1) - xx) : width;
            int hh = clamp ? Maths.Range(height, 0, m.GetLength(0) - yy) : height;

            Complex32[,] array = new Complex32[hh, ww];

            for (int i = 0; i < hh; i++)
            {
                var iyy = i + yy;

                if (iyy < 0 || iyy >= _height) continue;

                for (int j = 0; j < ww; j++)
                {
                    var jxx = j + xx;

                    if (jxx < 0 || jxx >= _width) continue;

                    array[i, j] = m[iyy, jxx];
                }
            }

            return array;
        }
        #endregion

        #region Merge voids
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Vector</returns>
        public static float[] Merge(this float[] a, float[] b)
        {
            return Merge(a, b, 0, b.GetLength(0));
        }
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <param name="start">Start position</param>
        /// <param name="length">Length</param>
        /// <returns>Vector</returns>
        public static float[] Merge(this float[] a, float[] b, int start, int length)
        {
            float[] c = Resize(b, length);
            float[] d = (float[])a.Clone();
            int h = Math.Min(length, a.GetLength(0) - start);

            for (int i = start; i < h; i++)
            {
                d[i] = c[i - start];
            }

            return d;
        }
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Vector</returns>
        public static Complex32[] Merge(this Complex32[] a, float[] b)
        {
            return Merge(a, b, 0, b.GetLength(0));
        }
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <param name="start">Start position</param>
        /// <param name="length">Length</param>
        /// <returns>Vector</returns>
        public static Complex32[] Merge(this Complex32[] a, float[] b, int start, int length)
        {
            float[] c = Resize(b, length);
            Complex32[] d = (Complex32[])a.Clone();
            int h = Math.Min(length, a.GetLength(0) - start);

            for (int i = start; i < h; i++)
            {
                d[i] = c[i - start];
            }

            return d;
        }
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <returns>Vector</returns>
        public static Complex32[] Merge(this Complex32[] a, Complex32[] b)
        {
            return Merge(a, b, 0, b.GetLength(0));
        }
        /// <summary>
        /// Merges two vectors.
        /// </summary>
        /// <param name="a">Vector</param>
        /// <param name="b">Vector</param>
        /// <param name="start">Start position</param>
        /// <param name="length">Length</param>
        /// <returns>Vector</returns>
        public static Complex32[] Merge(this Complex32[] a, Complex32[] b, int start, int length)
        {
            Complex32[] c = Resize(b, length);
            Complex32[] d = (Complex32[])a.Clone();
            int h = Math.Min(length, a.GetLength(0) - start);

            for (int i = start; i < h; i++)
            {
                d[i] = c[i - start];
            }

            return d;
        }

        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        public static float[,] Merge(this float[,] a, float[,] b)
        {
            return Merge(a, b, 0, 0, b.GetLength(0), b.GetLength(1));
        }
        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        /// <param name="y">Y</param>
        /// <param name="x">X</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <returns>Matrix</returns>
        public static float[,] Merge(this float[,] a, float[,] b, int y, int x, int height, int width)
        {
            float[,] c = Resize(b, height, width);
            float[,] d = (float[,])a.Clone();
            int h = Math.Min(height, a.GetLength(0) - y);
            int w = Math.Min(width, a.GetLength(1) - x);

            for (int i = y; i < h; i++)
            {
                for (int j = x; j < w; j++)
                {
                    d[i, j] = c[i - y, j - x];
                }
            }

            return d;
        }
        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        public static Complex32[,] Merge(this Complex32[,] a, float[,] b)
        {
            return Merge(a, b, 0, 0, b.GetLength(0), b.GetLength(1));
        }
        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        /// <param name="y">Y</param>
        /// <param name="x">X</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Merge(this Complex32[,] a, float[,] b, int y, int x, int height, int width)
        {
            float[,] c = Resize(b, height, width);
            Complex32[,] d = (Complex32[,])a.Clone();
            int h = Math.Min(height, a.GetLength(0) - y);
            int w = Math.Min(width, a.GetLength(1) - x);

            for (int i = y; i < h; i++)
            {
                for (int j = x; j < w; j++)
                {
                    d[i, j] = c[i - y, j - x];
                }
            }

            return d;
        }
        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        public static Complex32[,] Merge(this Complex32[,] a, Complex32[,] b)
        {
            return Merge(a, b, 0, 0, b.GetLength(0), b.GetLength(1));
        }
        /// <summary>
        /// Merges two matrices.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="b">Matrix</param>
        /// <param name="y">Y</param>
        /// <param name="x">X</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Merge(this Complex32[,] a, Complex32[,] b, int y, int x, int height, int width)
        {
            Complex32[,] c = Resize(b, height, width);
            Complex32[,] d = (Complex32[,])a.Clone();
            int h = Math.Min(height, a.GetLength(0) - y);
            int w = Math.Min(width, a.GetLength(1) - x);

            for (int i = y; i < h; i++)
            {
                for (int j = x; j < w; j++)
                {
                    d[i, j] = c[i - y, j - x];
                }
            }

            return d;
        }
        #endregion

        // MATLAB voids

        #region Abs/Angle
        /// <summary>
        /// Returns vector module.
        /// </summary>
        /// <param name="vector">Vector</param>
        /// <param name="squared">Squared or not</param>
        /// <returns>Value</returns>
        public static float Abs(this float[] vector, bool squared = false)
        {
            int length = vector.Length;
            float v = 0;

            for (int i = 0; i < length; i++)
            {
                v += vector[i] * vector[i];
            }

            if (squared)
                return v;

            return (float)Math.Sqrt(v);
        }
        /// <summary>
        /// Returns matrix module.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="squared">Squared or not</param>
        /// <returns>Vector</returns>
        public static float[] Abs(this float[,] matrix, bool squared = false)
        {
            int r = matrix.GetLength(0);
            int c = matrix.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = matrix[i, j];
                }

                v[i] = t.Abs(squared);
            }

            return v;
        }
        /// <summary>
        /// Returns vector module.
        /// </summary>
        /// <param name="vector">Vector</param>
        /// <param name="squared">Squared or not</param>
        /// <returns>Value</returns>
        public static Complex32 Abs(this Complex32[] vector, bool squared = false)
        {
            int length = vector.Length;
            Complex32 v = 0;

            for (int i = 0; i < length; i++)
            {
                v += vector[i] * vector[i];
            }

            if (squared)
                return v;

            return Maths.Sqrt(v);
        }
        /// <summary>
        /// Returns matrix module.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <param name="squared">Squared or not</param>
        /// <returns>Vector</returns>
        public static Complex32[] Abs(this Complex32[,] matrix, bool squared = false)
        {
            int r = matrix.GetLength(0);
            int c = matrix.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = matrix[i, j];
                }

                v[i] = t.Abs(squared);
            }

            return v;
        }

        /// <summary>
        /// Returns vector angle.
        /// </summary>
        /// <param name="vector">Vector</param>
        /// <returns>Value</returns>
        public static float Angle(this float[] vector)
        {
            return Angle(vector, vector);
        }
        /// <summary>
        /// Returns matrix module.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <returns>Vector</returns>
        public static float[] Angle(this float[,] matrix)
        {
            int r = matrix.GetLength(0);
            int c = matrix.GetLength(1);
            float[] v = new float[r];

            for (int i = 0; i < r; i++)
            {
                float[] t = new float[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = matrix[i, j];
                }

                v[i] = t.Angle();
            }

            return v;
        }
        /// <summary>
        /// Returns vector angle.
        /// </summary>
        /// <param name="vector">Vector</param>
        /// <returns>Value</returns>
        public static Complex32 Angle(this Complex32[] vector)
        {
            return Angle(vector, vector);
        }
        /// <summary>
        /// Returns matrix module.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        /// <returns>Vector</returns>
        public static Complex32[] Angle(this Complex32[,] matrix)
        {
            int r = matrix.GetLength(0);
            int c = matrix.GetLength(1);
            Complex32[] v = new Complex32[r];

            for (int i = 0; i < r; i++)
            {
                Complex32[] t = new Complex32[c];

                for (int j = 0; j < c; j++)
                {
                    t[j] = matrix[i, j];
                }

                v[i] = t.Angle();
            }

            return v;
        }

        #endregion

        #region Get/set rows and columns
        /// <summary>
        /// Returns the matrix column vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r">Column number</param>
        /// <returns>Array</returns>
        public static float[] GetCol(this float[,] m, int r)
        {
            int w = m.GetLength(0);
            float[] U = new float[w];
            for (int i = 0; i < w; i++)
            {
                U[i] = m[i, r];
            }
            return U;
        }
        /// <summary>
        /// Specifies the matrix column vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Array</param>
        /// <param name="r">Column number</param>
        /// <returns>Matrix</returns>
        public static float[,] SetCol(this float[,] m, float[] n, int r)
        {
            int w = m.GetLength(0);
            float[,] U = (float[,])m.Clone();
            for (int i = 0; i < w; i++)
            {
                U[i, r] = n[i];
            }
            return U;
        }
        /// <summary>
        /// Returns the row vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r">Row number</param>
        /// <returns>Array</returns>
        public static float[] GetRow(this float[,] m, int r)
        {
            int w = m.GetLength(1);
            float[] U = new float[w];
            for (int i = 0; i < w; i++)
            {
                U[i] = m[r, i];
            }
            return U;
        }
        /// <summary>
        /// Specifies the row vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Array</param>
        /// <param name="r">Row number</param>
        /// <returns>Matrix</returns>
        public static float[,] SetRow(this float[,] m, float[] n, int r)
        {
            int w = m.GetLength(1);
            float[,] U = (float[,])m.Clone();
            for (int i = 0; i < w; i++)
            {
                U[r, i] = n[i];
            }
            return U;
        }
        /// <summary>
        /// Returns the matrix column vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r">Column number</param>
        /// <returns>Array</returns>
        public static Complex32[] GetCol(this Complex32[,] m, int r)
        {
            int w = m.GetLength(0);
            Complex32[] U = new Complex32[w];
            for (int i = 0; i < w; i++)
            {
                U[i] = m[i, r];
            }
            return U;
        }
        /// <summary>
        /// Specifies the matrix column vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Array</param>
        /// <param name="r">Column number</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] SetCol(this Complex32[,] m, Complex32[] n, int r)
        {
            int w = m.GetLength(0);
            Complex32[,] U = (Complex32[,])m.Clone();
            for (int i = 0; i < w; i++)
            {
                U[i, r] = n[i];
            }
            return U;
        }
        /// <summary>
        /// Returns the row vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r">Row number</param>
        /// <returns>Array</returns>
        public static Complex32[] GetRow(this Complex32[,] m, int r)
        {
            int w = m.GetLength(1);
            Complex32[] U = new Complex32[w];
            for (int i = 0; i < w; i++)
            {
                U[i] = m[r, i];
            }
            return U;
        }
        /// <summary>
        /// Specifies the row vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Array</param>
        /// <param name="r">Row number</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] SetRow(this Complex32[,] m, Complex32[] n, int r)
        {
            int w = m.GetLength(1);
            Complex32[,] U = (Complex32[,])m.Clone();
            for (int i = 0; i < w; i++)
            {
                U[r, i] = n[i];
            }
            return U;
        }
        #endregion

        #region Diff voids
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Order</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        public static float[,] Diff(this float[,] a, int n, Direction direction, bool reverse = false)
        {
            // start
            int i, r, c;

            // direction of processing
            if (direction == Direction.Horizontal)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    c = a.GetLength(1);

                    if (c == 0)
                        break;

                    a = DiffHorizontal(a, reverse);
                }
            }
            else if (direction == Direction.Vertical)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);

                    if (r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                }
            }
            else
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);
                    c = a.GetLength(1);

                    if (c == 0 || r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                    a = DiffHorizontal(a, reverse);
                }
            }

            return a;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static float[,] DiffVertical(float[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new float[0, m];

            // new array
            float[,] y = new float[r, m];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i - 1, k] - a[i, k];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i, k] - a[i - 1, k];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static float[,] DiffHorizontal(float[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new float[m, 0];

            // new array
            float[,] y = new float[m, c];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i - 1] - a[k, i];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i] - a[k, i - 1];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Order</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Diff(this Complex32[,] a, int n, Direction direction, bool reverse = false)
        {
            // start
            int i, r, c;

            // direction of processing
            if (direction == Direction.Horizontal)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    c = a.GetLength(1);

                    if (c == 0)
                        break;

                    a = DiffHorizontal(a, reverse);
                }
            }
            else if (direction == Direction.Vertical)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);

                    if (r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                }
            }
            else
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);
                    c = a.GetLength(1);

                    if (c == 0 || r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                    a = DiffHorizontal(a, reverse);
                }
            }

            return a;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] DiffVertical(Complex32[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new Complex32[0, m];

            // new array
            Complex32[,] y = new Complex32[r, m];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i - 1, k] - a[i, k];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i, k] - a[i - 1, k];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static Complex32[,] DiffHorizontal(Complex32[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new Complex32[m, 0];

            // new array
            Complex32[,] y = new Complex32[m, c];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i - 1] - a[k, i];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i] - a[k, i - 1];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="n">Order</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Array</returns>
        public static float[] Diff(this float[] v, int n, bool reverse = false)
        {
            // start
            float[] z;
            float[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new float[0];

                y = new float[length];

                if (reverse)
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i - 1] - z[i];
                }
                else
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i] - z[i - 1];
                }
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="n">Order</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Array</returns>
        public static Complex32[] Diff(this Complex32[] v, int n, bool reverse = false)
        {
            // start
            Complex32[] z;
            Complex32[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new Complex32[0];

                y = new Complex32[length];

                if (reverse)
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i - 1] - z[i];
                }
                else
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i] - z[i - 1];
                }
            }

            return y;
        }
        #endregion

        #region Reshape voids
        /// <summary>
        /// Returns a matrix formed from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="height">Height</param>
        /// <returns>Matrix</returns>
        public static float[,] Reshape(this float[] a, int height)
        {
            int n = a.Length;
            int width = (int)Math.Ceiling(n / (float)height);
            float[,] H = new float[height, width];
            int i, j, k;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    k = j * height + i;

                    if (k < n)
                    {
                        H[i, j] = a[k];
                    }
                    else break;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a matrix formed from a vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="height">Height</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Reshape(this Complex32[] a, int height)
        {
            int n = a.Length;
            int width = (int)Math.Ceiling(n / (float)height);
            Complex32[,] H = new Complex32[height, width];
            int i, j, k;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    k = j * height + i;

                    if (k < n)
                    {
                        H[i, j] = a[k];
                    }
                    else break;
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a vector formed from a matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Reshape(this float[,] a, int length)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            float[] v = new float[length];

            int i, j, k;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    k = j * height + i;

                    if (k < length)
                    {
                        v[j * height + i] = a[i, j];
                    }
                    else break;
                }
            }

            return v;
        }
        /// <summary>
        /// Returns a vector formed from a matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static Complex32[] Reshape(this Complex32[,] a, int length)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            Complex32[] v = new Complex32[length];

            int i, j, k;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    k = j * height + i;

                    if (k < length)
                    {
                        v[j * height + i] = a[i, j];
                    }
                    else break;
                }
            }

            return v;
        }
        #endregion

        #region Diagonalize
        /// <summary>
        /// Implements the reduction of a vector to a diagonal matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Diag(this float[] v)
        {
            int n = v.Length, i;
            float[,] diag = new float[n, n];

            for (i = 0; i < n; i++)
            {
                diag[i, i] = v[i];
            }
            return diag;
        }
        /// <summary>
        /// Implements the reduction of a vector to a diagonal matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Diag(this Complex32[] v)
        {
            int n = v.Length, i;
            Complex32[,] diag = new Complex32[n, n];

            for (i = 0; i < n; i++)
            {
                diag[i, i] = v[i];
            }
            return diag;
        }
        /// <summary>
        /// Returns a vector whose elements lie on the diagonal of the matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <returns>Array</returns>
        public static float[] Diag(this float[,] a)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            float[] v = new float[height];
            int i;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                if (i < width)
                {
                    v[i] = a[i, i];
                }
            }

            return v;
        }
        /// <summary>
        /// Returns a vector whose elements lie on the diagonal of the matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Diag(this Complex32[,] a)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            Complex32[] v = new Complex32[height];
            int i;

            // reshaping:
            for (i = 0; i < height; i++)
            {
                if (i < width)
                {
                    v[i] = a[i, i];
                }
            }

            return v;
        }
        #endregion

        #region Swap voids
        /// <summary>
        /// Implements a permutation of the vectors of the matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="i">First row or column</param>
        /// <param name="j">Second row or column</param>
        /// <param name="direction">Processing direction</param>
        public static void Swap(this float[,] a, int i, int j, Direction direction = Direction.Horizontal)
        {
            // properties:
            float temp;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z;

            // by rows:
            if (direction == Direction.Horizontal)
            {
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
            }
            // by columns:
            else if (direction == Direction.Vertical)
            {
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
            // by rows and columns
            else
            {
                // rows:
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
                // columns:
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
        }
        /// <summary>
        /// Implements a permutation of the vectors of the matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="i">First row or column</param>
        /// <param name="j">Second row or column</param>
        /// <param name="direction">Processing direction</param>
        public static void Swap(this Complex32[,] a, int i, int j, Direction direction = Direction.Horizontal)
        {
            // properties:
            Complex32 temp;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z;

            // by rows:
            if (direction == Direction.Horizontal)
            {
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
            }
            // by columns:
            else if (direction == Direction.Vertical)
            {
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
            // by rows and columns
            else
            {
                // rows:
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
                // columns:
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
        }
        /// <summary>
        /// Implements a permutation of the elements of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="i">First element position</param>
        /// <param name="j">Second element position</param>
        public static void Swap(this float[] v, int i, int j)
        {
            // get elements:
            float e1 = v[i], e2 = v[j];
            // swapping vector elements:
            v[j] = e1; v[i] = e2;
        }
        /// <summary>
        /// Implements a permutation of the elements of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="i">First element position</param>
        /// <param name="j">Second element position</param>
        public static void Swap(this Complex32[] v, int i, int j)
        {
            // get elements:
            Complex32 e1 = v[i], e2 = v[j];
            // swapping vector elements:
            v[j] = e1; v[i] = e2;
        }
        #endregion

        #region Remove voids
        /// <summary>
        /// Implements the removal of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="i">First row or column</param>
        /// <param name="length">Length</param>
        /// <param name="direction">Processing direction</param>
        /// <returns>Matrix</returns>
        public static float[,] Remove(this float[,] a, int i, int length, Direction direction = Direction.Horizontal)
        {
            // properties:
            float[,] H;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z, y;

            // by columns:
            if (direction == Direction.Vertical)
            {
                H = new float[col, length];

                for (int j = 0; j < col; j++)
                {
                    for (z = 0, y = i; z < length; z++, y++)
                    {
                        H[j, z] = a[j, y];
                    }
                }
            }
            // by rows:
            else if (direction == Direction.Horizontal)
            {
                H = new float[length, row];

                for (int j = 0; j < row; j++)
                {
                    for (z = 0, y = i; z < length; z++, y++)
                    {
                        H[z, j] = a[y, j];
                    }
                }
            }
            // by columns and rows:
            else
            {
                H = new float[length, length];
                int w, t;

                for (z = 0, y = i; z < length; z++, y++)
                {
                    for (w = 0, t = i; w < length; w++, t++)
                    {
                        H[z, w] = a[y, t];
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the removal of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="i">First row or column</param>
        /// <param name="length">Length</param>
        /// <param name="direction">Processing direction</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Remove(this Complex32[,] a, int i, int length, Direction direction = Direction.Horizontal)
        {
            // properties:
            Complex32[,] H;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z, y;

            // by columns:
            if (direction == Direction.Vertical)
            {
                H = new Complex32[col, length];

                for (int j = 0; j < col; j++)
                {
                    for (z = 0, y = i; z < length; z++, y++)
                    {
                        H[j, z] = a[j, y];
                    }
                }
            }
            // by rows:
            else if (direction == Direction.Horizontal)
            {
                H = new Complex32[length, row];

                for (int j = 0; j < row; j++)
                {
                    for (z = 0, y = i; z < length; z++, y++)
                    {
                        H[z, j] = a[y, j];
                    }
                }
            }
            // by columns and rows:
            else
            {
                H = new Complex32[length, length];
                int w, t;

                for (z = 0, y = i; z < length; z++, y++)
                {
                    for (w = 0, t = i; w < length; w++, t++)
                    {
                        H[z, w] = a[y, t];
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the removal of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="i">Number of element</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Remove(this float[] v, int i, int length)
        {
            float[] w = new float[length];

            for (int z = 0, y = i; z < length; z++, y++)
            {
                w[z] = v[y];
            }
            return w;
        }
        /// <summary>
        /// Implements the removal of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="i">Number of element</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static Complex32[] Remove(this Complex32[] v, int i, int length)
        {
            Complex32[] w = new Complex32[length];

            for (int z = 0, y = i; z < length; z++, y++)
            {
                w[z] = v[y];
            }
            return w;
        }
        #endregion

        #region Minor voids
        /// <summary>
        /// Implements the operation of taking the minor of the matrix.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <param name="n">Row and column number</param>
        /// <returns>Square matrix</returns>
        public static float[,] Minor(this float[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new ArgumentException("The matrix must be square");
            if (n >= height || n < 0) throw new ArgumentException("Row and column number specified is invalid");

            // new matrix:
            float[,] H = new float[height - 1, width - 1];
            int i, j, x = 0, y = 0;

            // dub:
            for (i = 0; i < height; i++)
            {
                if (i == n)
                    continue;
                else
                {
                    for (j = 0; j < width; j++)
                    {
                        if (j == n)
                            continue;

                        else { H[y, x] = m[i, j]; x++; }
                    }
                    x = 0;
                    y++;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the operation of taking the minor of the matrix.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <param name="n">Row and column number</param>
        /// <returns>Square matrix</returns>
        public static Complex32[,] Minor(this Complex32[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new ArgumentException("The matrix must be square");
            if (n >= height || n < 0) throw new ArgumentException("Row and column number specified is invalid");

            // new matrix:
            Complex32[,] H = new Complex32[height - 1, width - 1];
            int i, j, x = 0, y = 0;

            // dub:
            for (i = 0; i < height; i++)
            {
                if (i == n)
                    continue;
                else
                {
                    for (j = 0; j < width; j++)
                    {
                        if (j == n)
                            continue;

                        else { H[y, x] = m[i, j]; x++; }
                    }
                    x = 0;
                    y++;
                }
            }
            return H;
        }
        #endregion

        #region Extend voids
        /// <summary>
        /// Extends the vector to the specified length.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Extend(this float[] v, int length)
        {
            int r0 = v.GetLength(0);
            int rr = (length - r0) / 2;
            int dr = length - rr;
            float[] b = new float[length];
            int i;

            for (i = 0; i < rr; i++)
                b[i] = v[rr - i];

            for (i = 0; i < r0; i++)
                b[i + rr] = v[i];

            for (i = 0; i <= rr; i++)
                b[i + dr - 1] = v[r0 - i - 1];

            return b;
        }
        /// <summary>
        /// Extends the matrix to the specified size.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="height">Height</param>
        /// <param name="width">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] Extend(this float[,] m, int height, int width)
        {
            int r = m.GetLength(0);
            int c = m.GetLength(1);

            if (height > r)
                m = Matrice.ExtendVertical(m, height);
            if (width > c)
                m = Matrice.ExtendHorizontal(m, width);

            return m;
        }
        /// <summary>
        /// extend vertical.
        /// </summary>
        /// <param name="m"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private static float[,] ExtendVertical(float[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int rr = (length - r0) / 2;
            float[,] B = new float[length, c0];
            int i, j;

            // do job
            for (i = 0; i < length; i++)
                for (j = 0; j < c0; j++)
                    B[i, j] = m[Maths.Mod(rr - i - 1, r0), Maths.Mod(j, c0)];


            for (i = 0; i < r0; i++)
                for (j = 0; j < c0; j++)
                    B[rr + i, j] = m[i, j];

            return B;
        }
        /// <summary>
        /// extend horizontal.
        /// </summary>
        /// <param name="m"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private static float[,] ExtendHorizontal(float[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int cc = (length - c0) / 2;
            float[,] B = new float[r0, length];
            int i, j;

            // do job
            for (i = 0; i < r0; i++)
                for (j = 0; j < length; j++)
                    B[i, j] = m[Maths.Mod(i, r0), Maths.Mod(cc - j - 1, c0)];


            for (i = 0; i < r0; i++)
                for (j = 0; j < c0; j++)
                    B[i, cc + j] = m[i, j];

            return B;
        }

        /// <summary>
        /// Extends the vector to the specified length.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static Complex32[] Extend(this Complex32[] v, int length)
        {
            int r0 = v.GetLength(0);
            int rr = (length - r0) / 2;
            int dr = length - rr;
            Complex32[] b = new Complex32[length];
            int i;

            for (i = 0; i < rr; i++)
                b[i] = v[rr - i];

            for (i = 0; i < r0; i++)
                b[i + rr] = v[i];

            for (i = 0; i <= rr; i++)
                b[i + dr - 1] = v[r0 - i - 1];

            return b;
        }
        /// <summary>
        /// Extends the matrix to the specified size.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="height">Height</param>
        /// <param name="width">Width</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Extend(this Complex32[,] m, int height, int width)
        {
            int r = m.GetLength(0);
            int c = m.GetLength(1);

            if (height > r)
                m = Matrice.ExtendVertical(m, height);
            if (width > c)
                m = Matrice.ExtendHorizontal(m, width);

            return m;
        }
        /// <summary>
        /// extend vertical.
        /// </summary>
        /// <param name="m"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private static Complex32[,] ExtendVertical(Complex32[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int rr = (length - r0) / 2;
            Complex32[,] B = new Complex32[length, c0];
            int i, j;

            // do job
            for (i = 0; i < length; i++)
                for (j = 0; j < c0; j++)
                    B[i, j] = m[Maths.Mod(rr - i - 1, r0), Maths.Mod(j, c0)];


            for (i = 0; i < r0; i++)
                for (j = 0; j < c0; j++)
                    B[rr + i, j] = m[i, j];

            return B;
        }
        /// <summary>
        /// extend horizontal.
        /// </summary>
        /// <param name="m"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private static Complex32[,] ExtendHorizontal(Complex32[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int cc = (length - c0) / 2;
            Complex32[,] B = new Complex32[r0, length];
            int i, j;

            // do job
            for (i = 0; i < r0; i++)
                for (j = 0; j < length; j++)
                    B[i, j] = m[Maths.Mod(i, r0), Maths.Mod(cc - j - 1, c0)];


            for (i = 0; i < r0; i++)
                for (j = 0; j < c0; j++)
                    B[i, cc + j] = m[i, j];

            return B;
        }
        #endregion

        // Extra voids

        #region Compute methods
        /// <summary>
        /// Returns an array of function values.
        /// </summary>
        /// <param name="min">Minimum</param>
        /// <param name="max">Maximum</param>
        /// <param name="step">Step</param>
        /// <returns>Array</returns>
        public static float[] Compute(float min, float max, float step)
        {
            // ******************************
            // MATLAB vector computing void
            // designed by Valery Asiryan
            // ******************************

            // shifts and variables:
            float dy = max - min + step, i;
            int dx = (int)Maths.Round(dy / step), j;

            // C# has a significant bug, which you can check with:
            // min = 0.5, max = 1, step = 0.001,
            // maxz = max, j = 63.

            // output vector and eps:
            float[] x = new float[dx];
            float eps = max / 1e8f, maxz;

            // compute:
            if (min > max)
            {
                // limit value:
                maxz = max - eps;

                // for arrays like [6,5...-5,-6]:
                for (j = 0, i = min; i >= maxz; i += step)
                {
                    if (j < dx)
                    {
                        x[j] = i; j++;
                    }
                    else break;
                }
            }
            else
            {
                // limit value:
                maxz = max + eps;

                // for arrays like [-6,-5...5,6]:
                for (j = 0, i = min; i <= maxz; i += step)
                {
                    if (j < dx)
                    {
                        x[j] = i; j++;
                    }
                    else break;
                }
            }
            return x;
        }
        /// <summary>
        /// Returns an array of function values.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Array</returns>
        public static float[] Compute(this float[] v, IFloat function)
        {
            int length = v.Length;
            float[] H = new float[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = function(v[i]);
            }
            return H;
        }
        /// <summary>
        /// Returns an array of function values.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Array</returns>
        public static Complex32[] Compute(this Complex32[] v, IComplex function)
        {
            int length = v.Length;
            Complex32[] H = new Complex32[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = function(v[i]);
            }
            return H;
        }
        /// <summary>
        /// Returns a matrix of function values.
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static float[,] Compute(this float[] x, float[] y, IFloatMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            float[,] z = new float[xlength, ylength];
            int i, j;

            for (i = 0; i < xlength; i++)
            {
                for (j = 0; j < ylength; j++)
                {
                    z[i, j] = function(x[i], y[i]);
                }
            }
            return z;
        }
        /// <summary>
        /// Returns a matrix of function values.
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Compute(this float[] x, Complex32[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex32[,] z = new Complex32[xlength, ylength];
            int i, j;

            for (i = 0; i < xlength; i++)
            {
                for (j = 0; j < ylength; j++)
                {
                    z[i, j] = function(x[i], y[i]);
                }
            }
            return z;
        }
        /// <summary>
        /// Returns a matrix of function values.
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Compute(this Complex32[] x, float[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex32[,] z = new Complex32[xlength, ylength];
            int i, j;

            for (i = 0; i < xlength; i++)
            {
                for (j = 0; j < ylength; j++)
                {
                    z[i, j] = function(x[i], y[i]);
                }
            }
            return z;
        }
        /// <summary>
        /// Returns a matrix of function values.
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Compute(this Complex32[] x, Complex32[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex32[,] z = new Complex32[xlength, ylength];
            int i, j;

            for (i = 0; i < xlength; i++)
            {
                for (j = 0; j < ylength; j++)
                {
                    z[i, j] = function(x[i], y[i]);
                }
            }
            return z;
        }
        /// <summary>
        /// Returns an array of function values.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static float[,] Compute(this float[,] m, IFloat function)
        {
            int i, j;
            int ml = m.GetLength(1), mr = m.GetLength(0);
            float[,] H = new float[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = function(m[i, j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Returns a matrix of function values.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="function">Continuous function delegate</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Compute(this Complex32[,] m, IComplex function)
        {
            int i, j;
            int ml = m.GetLength(1), mr = m.GetLength(0);
            Complex32[,] H = new Complex32[mr, ml];

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    H[i, j] = function(m[i, j]);
                }
            }
            return H;
        }
        #endregion

        #region Radius vector
        /// <summary>
        /// Implements the construction of a vector of ones.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static float[] One(int n)
        {
            float[] v = new float[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = 1;
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a vector of zeros.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static float[] Zero(int n)
        {
            return new float[n];
        }
        #endregion

        #region Matrix products
        /// <summary>
        /// Returns the Householder vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static float[] Householder(this float[] v)
        {
            int length = v.Length;

            // length checking:
            if (length > 0)
            {
                float norm = v.Norm();
                float[] u = new float[length];

                // householder vector:
                if (norm != 0)
                {
                    u[0] = v[0] / norm;
                    u[0] = u[0] + Maths.Sign(u[0]);
                    u[0] = u[0] / (float)Math.Sqrt(Math.Abs(u[0]));

                    for (int i = 1; i < length; i++)
                    {
                        u[i] = v[i] / norm;
                        u[i] = u[i] / u[0];
                    }
                }
                else
                {
                    u = v;
                    u[0] = (float)Maths.Sqrt2;
                }
                return u;
            }
            return v;
        }
        /// <summary>
        /// Implements the construction of the companion matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Companion(this float[] v)
        {
            int n = v.Length, i;
            float[,] H = new float[n, n];

            // last column:
            for (i = 0; i < n; i++)
            {
                H[0, i] = -v[i];
            }
            // eyes:
            for (i = 1; i < n; i++)
            {
                H[i, i - 1] = 1;
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Vandermond matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Vander(this float[] v)
        {
            int n = v.Length;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = (float)Math.Pow(v[i], j);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of an incomplete Hankel matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Hankeli(this float[] v)
        {
            int n = v.Length;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n - i; j++)
                {
                    H[j, i] = v[i + j];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Hankel matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Hankel(this float[] v)
        {
            int n = v.Length / 2;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[i + j];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Toeplitz matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Toeplitz(this float[] v)
        {
            int n = v.Length / 2;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Cauchy matrix.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Cauchy(this float[] x, float[] y)
        {
            int m = x.Length, l = y.Length;
            float[,] H = new float[m, l];
            float kern;
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    kern = x[i] - y[j];
                    H[i, j] = (kern != 0) ? 1.0f / kern : 0;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a circulant matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Circulant(this float[] v)
        {
            int n = v.Length;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a symmetric matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static float[,] Symmetric(this float[] v)
        {
            int n = v.Length;
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
                for (j = 0; j < i; j++)
                {
                    H[i, j] = H[j, i];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of the companion matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Companion(this Complex32[] v)
        {
            int n = v.Length, i;
            Complex32[,] H = new Complex32[n, n];

            // last column:
            for (i = 0; i < n; i++)
            {
                H[i, n - 1] = -v[i];
            }
            // eyes:
            for (i = 1; i < n; i++)
            {
                H[i, i - 1] = 1;
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Vandermond matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Vander(this Complex32[] v)
        {
            int n = v.Length;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Pow(v[i], j);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of an incomplete Hankel matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Hankeli(this Complex32[] v)
        {
            int n = v.Length;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n - i; j++)
                {
                    H[j, i] = v[i + j];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Hankel matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Hankel(this Complex32[] v)
        {
            int n = v.Length / 2;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[i + j];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Toeplitz matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Toeplitz(this Complex32[] v)
        {
            int n = v.Length / 2;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Cauchy matrix.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Cauchy(this Complex32[] x, Complex32[] y)
        {
            int m = x.Length, l = y.Length;
            Complex32[,] H = new Complex32[m, l];
            Complex32 kern;
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    kern = x[i] - y[j];
                    H[i, j] = (kern != 0) ? 1.0 / kern : 0;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a circulant matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Circulant(this Complex32[] v)
        {
            int n = v.Length;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a symmetric matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Symmetric(this Complex32[] v)
        {
            int n = v.Length;
            Complex32[,] H = new Complex32[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    H[i, j] = v[Maths.Mod(j - i, n)];
                }
                for (j = 0; j < i; j++)
                {
                    H[i, j] = H[j, i];
                }
            }
            return H;
        }
        #endregion

        #region Radius matrix
        /// <summary>
        /// Implements the construction of a zero matrix.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] Zero(int m, int l)
        {
            return new float[m, l];
        }
        /// <summary>
        /// Implements the construction of a eye matrix.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] Eye(int m, int l)
        {
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = i == j ? 1 : 0;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of ones.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] One(int m, int l)
        {
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = 1;
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of the exchange matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Exchange(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = (j == n - i - 1) ? 1 : 0;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Lehmer matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Lehmer(int n)
        {
            float[,] H = new float[n, n];
            int i, j;
            float x, y;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    x = (float)i; y = (float)j;
                    H[i - 1, j - 1] = (j >= i) ? x / y : y / x;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of the Redheffer matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Redheffer(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    if (i == j)
                    {
                        H[i - 1, j - 1] = 1;
                    }
                    if (i % j == 0)
                    {
                        H[i - 1, j - 1] = 1;
                    }
                    else if (i == 1 || j == 1)
                    {
                        H[i - 1, j - 1] = 1;
                    }
                    else
                    {
                        H[i - 1, j - 1] = 0;
                    }
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a Hilbert matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Hilbert(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    H[i - 1, j - 1] = 1.0f / (i + j - 1);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a cyclic matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Circulant(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Mod(j - i, n);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a symmetric matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] Symmetric(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    H[i, j] = Maths.Mod(j - i, n);
                }
                for (j = 0; j < i; j++)
                {
                    H[i, j] = H[j, i];
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of GCD.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static float[,] GCD(int n)
        {
            float[,] H = new float[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Maths.Gcd(i + 1, j + 1);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of the Stirling matrix of the first or second kind.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <param name="second">Second kind or not</param>
        /// <returns>Matrix</returns>
        public static float[,] Stirling(int n, bool second = false)
        {
            // Stirling's matrix 
            // of the second kind
            float[,] S = new float[n, n];
            int i, j;

            // S(n, 0) = 0
            for (i = 0; i < n; i++)
                S[i, 0] = 0;

            // S(n, n) = 1
            for (i = 0; i < n; i++)
                S[i, i] = 1;

            // Create Stirling's matrix?
            if (n > 0)
            {
                // Second kind or not?
                if (second)
                {
                    // S(n, m) = S(n - 1, m - 1) + m * S(n - 1, m):
                    for (i = 1; i < n; i++)
                    {
                        for (j = 1; j <= i - 1; j++)
                        {
                            S[i, j] = S[i - 1, j - 1] + j * S[i - 1, j];
                        }
                    }
                }
                else
                {
                    // S(n, m) = S(n - 1, m - 1) + (n - 1) * S(n - 1, m):
                    for (i = 1; i < n; i++)
                    {
                        for (j = 1; j <= i - 1; j++)
                        {
                            S[i, j] = S[i - 1, j - 1] + (i - 1) * S[i - 1, j];
                        }
                    }
                }
            }

            return S;
        }
        /// <summary>
        /// Implements the construction of a magic square.
        /// </summary>
        /// <param name="n">Size (odd number)</param>
        /// <returns>Matrix</returns>
        public static float[,] Magic(int n)
        {
            if (Maths.Mod(n, 2) != 1)
                throw new ArgumentException("Dimension of the matrix must be an odd number");

            float[,] m = new float[n, n];
            int i, j;

            // MATLAB construction
            // for magic square:
            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    m[i - 1, j - 1] = n * Maths.Mod(i + j - (n + 3) / 2, n) + Maths.Mod(i + 2 * (j - 1), n) + 1;
                }
            }
            return m;
        }
        #endregion

        #region Random matrices and vectors

        /// <summary>
        /// 
        /// </summary>
        private static readonly Random rnd = new Random();

        #region Rand and randc voids
        /// <summary>
        /// Implements the construction of a vector of random numbers, the values of which are distributed UMapxing to a uniform distribution.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static float[] Rand(int n)
        {
            float[] v = new float[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = (float)rnd.NextDouble();
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a vector of random numbers, the values of which are distributed UMapxing to a uniform distribution.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static Complex32[] Randc(int n)
        {
            Complex32[] v = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = new Complex32((float)rnd.NextDouble(), (float)rnd.NextDouble());
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a matrix of random numbers, the values of which are distributed UMapxing to a uniform distribution.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] Rand(int m, int l)
        {
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = (float)rnd.NextDouble();
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of random numbers, the values of which are distributed UMapxing to a uniform distribution.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Randc(int m, int l)
        {
            Complex32[,] H = new Complex32[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = new Complex32((float)rnd.NextDouble(), (float)rnd.NextDouble());
                }
            }

            return H;
        }
        #endregion

        #region Randi and randic voids
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static float[] Randi(int n)
        {
            return Randi(n, 1, n + 1);
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static Complex32[] Randic(int n)
        {
            return Randic(n, 1, n + 1);
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Array</returns>
        public static float[] Randi(int n, int a, int b)
        {
            float[] v = new float[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = rnd.Next(a, b);
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Array</returns>
        public static Complex32[] Randic(int n, int a, int b)
        {
            Complex32[] v = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = new Complex32(rnd.Next(a, b), rnd.Next(a, b));
            }

            return v;
        }

        /// <summary>
        /// Implements the construction of a matrix of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[,] Randi(int m, int l)
        {
            return Randi(m, l, 1, l + 1);
        }
        /// <summary>
        /// Implements the construction of a matrix of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Matrix</returns>
        public static float[,] Randi(int m, int l, int a, int b)
        {
            float[,] H = new float[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = rnd.Next(a, b);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Randic(int m, int l)
        {
            return Randic(m, l, 1, l + 1);
        }
        /// <summary>
        /// Implements the construction of a matrix of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Randic(int m, int l, int a, int b)
        {
            Complex32[,] H = new Complex32[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = new Complex32(rnd.Next(a, b), rnd.Next(a, b));
                }
            }

            return H;
        }
        #endregion
        #endregion

        #region Parse methods
        /// <summary>
        /// Parses the original string into a matrix of float numbers.
        /// <remarks>
        /// Example: "[1, 2, 3; 4, 5, 6; 7, 8, 9]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static float[,] Parse(this float[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
            int r = rows.Length, n = nums.Length, k;
            float[,] H = new float[r, n];
            int i, j;

            // first row
            for (j = 0; j < n; j++)
            {
                H[0, j] = float.Parse(nums[j]);
            }

            // other rows
            for (i = 1; i < r; i++)
            {
                nums = rows[i].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
                k = Math.Min(n, nums.Length);

                for (j = 0; j < k; j++)
                {
                    H[i, j] = float.Parse(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of float numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, out float[,] result)
        {
            float[,] zero = null;
            try
            {
                result = Matrice.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }
        /// <summary>
        /// Parses the original string into a matrix of complex numbers.
        /// </summary>
        /// <remarks>
        /// Example: "[1 + 2i, 2 + 4i; 3 + 6i, 4 + 8i]";
        /// </remarks>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static Complex32[,] Parse(this Complex32[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
            int r = rows.Length, n = nums.Length, k;
            Complex32[,] H = new Complex32[r, n];
            int i, j;

            // first row
            for (j = 0; j < n; j++)
            {
                H[0, j] = StringOptions.Compar(nums[j]);
            }

            // other rows
            for (i = 1; i < r; i++)
            {
                nums = rows[i].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
                k = Math.Min(n, nums.Length);

                for (j = 0; j < k; j++)
                {
                    H[i, j] = StringOptions.Compar(nums[j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of complex numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, out Complex32[,] result)
        {
            Complex32[,] zero = null;
            try
            {
                result = Matrice.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }

        /// <summary>
        /// Parses the original string into a vector of float numbers.
        /// <remarks>
        /// Example: "[1, 2, 3, 4]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static float[] Parse(this float[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
                string[] nums = rows[0].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
                int n = nums.Length, i;
                float[] H = new float[n];

                // collecting rows:
                for (i = 0; i < n; i++)
                {
                    H[i] = float.Parse(nums[i]);
                }
                return H;
            }
            else
            {
                throw new ArgumentException("The input string was in the wrong format");
            }
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of float numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, out float[] result)
        {
            float[] zero = null;
            try
            {
                result = Matrice.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }
        /// <summary>
        /// Parses the original string into a vector of complex numbers.
        /// <remarks>
        /// Example: "[1 + 2i, 2 + 0.3i, 3 + i, 4 - 11i]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static Complex32[] Parse(this Complex32[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
                string[] nums = rows[0].Split(new char[] { '|' }, StringSplitOptions.RemoveEmptyEntries);
                int n = nums.Length, i;
                Complex32[] H = new Complex32[n];

                // collecting rows:
                for (i = 0; i < n; i++)
                {
                    H[i] = StringOptions.Compar(nums[i]);
                }
                return H;
            }
            else
            {
                throw new ArgumentException("The input string was in the wrong format");
            }
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of complex numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, out Complex32[] result)
        {
            Complex32[] zero = null;
            try
            {
                result = Matrice.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }
        #endregion

        // Solvers

        #region Gauss-Jordan elimination
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Extended matrix</param>
        /// <returns>Array</returns>
        public static float[] Solve(this float[,] A)
        {
            int N = A.GetLength(0);
            float[,] Q = new float[N, N];
            float[] b = new float[N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Q[i, j] = A[i, j];
                }
                b[i] = A[i, N];
            }

            return Q.Solve(b);
        }
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Solve(this float[,] A, float[] b)
        {
            // Input data
            if (!IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            int M = A.GetLength(0);
            int N = b.GetLength(0);

            if (N != M)
                throw new ArgumentException("Vector length should be equal to the height of the matrix");

            float[][] a = Jagged.ToJagged(A);
            float[] q = (float[])b.Clone();
            float eps = 1e-16f;

            // method of Gauss 
            for (int p = 0; p < N; p++)
            {

                int max = p;
                for (int i = p + 1; i < N; i++)
                {
                    if (Math.Abs(a[i][p]) > Math.Abs(a[max][p]))
                    {
                        max = i;
                    }
                }
                float[] temp = a[p]; a[p] = a[max]; a[max] = temp;
                float t = q[p]; q[p] = q[max]; q[max] = t;

                if (Math.Abs(a[p][p]) <= eps)
                {
                    return b;
                }

                for (int i = p + 1; i < N; i++)
                {
                    float alpha = a[i][p] / a[p][p];
                    q[i] -= alpha * q[p];
                    for (int j = p; j < N; j++)
                    {
                        a[i][j] -= alpha * a[p][j];
                    }
                }
            }

            // Result
            float[] x = new float[N];
            for (int i = N - 1; i >= 0; i--)
            {
                float sum = 0;
                for (int j = i + 1; j < N; j++)
                {
                    sum += a[i][j] * x[j];
                }
                x[i] = (q[i] - sum) / a[i][i];
            }
            return x;
        }
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Extended matrix</param>
        /// <returns>Array</returns>
        public static Complex32[] Solve(this Complex32[,] A)
        {
            int N = A.GetLength(0);
            Complex32[,] Q = new Complex32[N, N];
            Complex32[] b = new Complex32[N];

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Q[i, j] = A[i, j];
                }
                b[i] = A[i, N];
            }

            return Q.Solve(b);
        }
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex32[] Solve(this Complex32[,] A, Complex32[] b)
        {
            // Input data
            if (!IsSquare(A))
                throw new ArgumentException("The matrix must be square");

            int M = A.GetLength(0);
            int N = b.GetLength(0);

            if (N != M)
                throw new ArgumentException("Vector length should be equal to the height of the matrix");

            Complex32[][] a = Jagged.ToJagged(A);
            Complex32[] q = (Complex32[])b.Clone();
            float eps = 1e-16f;

            // method of Gauss 
            for (int p = 0; p < N; p++)
            {

                int max = p;
                for (int i = p + 1; i < N; i++)
                {
                    if (Maths.Abs(a[i][p]) > Maths.Abs(a[max][p]))
                    {
                        max = i;
                    }
                }
                Complex32[] temp = a[p]; a[p] = a[max]; a[max] = temp;
                Complex32 t = q[p]; q[p] = q[max]; q[max] = t;

                if (Maths.Abs(a[p][p]) <= eps)
                {
                    return b;
                }

                for (int i = p + 1; i < N; i++)
                {
                    Complex32 alpha = a[i][p] / a[p][p];
                    q[i] -= alpha * q[p];
                    for (int j = p; j < N; j++)
                    {
                        a[i][j] -= alpha * a[p][j];
                    }
                }
            }

            // Result
            Complex32[] x = new Complex32[N];
            for (int i = N - 1; i >= 0; i--)
            {
                Complex32 sum = 0;
                for (int j = i + 1; j < N; j++)
                {
                    sum += a[i][j] * x[j];
                }
                x[i] = (q[i] - sum) / a[i][i];
            }
            return x;
        }
        #endregion
    }
}
