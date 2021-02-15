using System;
using System.Threading.Tasks;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to implement standard algebraic operations on matrices and vectors.
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
        public static bool IsEquals(this Complex[,] m, Complex[,] n)
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
                if (Matrice.IsEquals(m.Transponate(), m.Negate()))
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
        public static bool IsVector(this Complex[,] m)
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
        public static bool IsSquare(this Complex[,] m)
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
        public static bool IsSymmetric(this Complex[,] m)
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
        public static bool IsSkewSymmetric(this Complex[,] m)
        {
            if (Matrice.IsSquare(m))
            {
                // ?A' = -A
                if (Matrice.IsEquals(m.Hermitian(), m.Negate()))
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
        public static bool IsDiagonal(this Complex[,] m)
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

        #region Matrix tranform
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

            return Jagged.FromJagged(LinealgOptions.Invert(Jagged.ToJagged(m)));
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
        public static Complex[,] Invert(this Complex[,] m)
        {
            if (m.GetLength(0) != m.GetLength(1))
            {
                return m;
            }

            return Jagged.FromJagged(LinealgOptions.Invert(Jagged.ToJagged(m)));
        }
        /// <summary>
        /// Implements the transpose of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Transponate(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r1, r0];
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
        public static Complex[,] Conjugate(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static Complex[,] Hermitian(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r1, r0];
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
        #endregion

        #region Matrix properties
        /// <summary>
        /// Returns the trace value of a square matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>float precision floating point number</returns>
        public static float Trace(this float[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("The matrix must be square");

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
        /// <returns>float precision floating point number</returns>
        public static float Det(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new Exception("The matrix must be square");

            unsafe
            {
                // copy array
                float[,] n = (float[,])m.Clone();

                fixed (float* pm = &n[0, 0])
                    return LinealgOptions.Determinant(pm, mr);
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
                throw new Exception("The matrix must be square");

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
        public static Complex Trace(this Complex[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("The matrix must be square");

            int d = m.GetLength(0);
            int i;
            Complex kernel = 0;

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
        /// <returns>float precision floating point number</returns>
        public static Complex Det(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new Exception("The matrix must be square");

            unsafe
            {
                // copy array
                Complex[,] n = (Complex[,])m.Clone();

                fixed (Complex* pm = &n[0, 0])
                    return LinealgOptions.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Returns the P-norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Matrix</returns>
        public static float Norm(this Complex[,] m, float p)
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
        public static float Norm(this Complex[,] m)
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
        public static Complex[,] Round(this Complex[,] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
            Complex c;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    c = m[i, j];
                    H[i, j] = new Complex((float)Math.Round(c.Real, digits, mode), (float)Math.Round(c.Imag, digits, mode));
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
        public static Complex[,] Kronecker(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex[,] H = new Complex[ml * nl, mr * nr];
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
        public static Complex[,] Kronecker(this Complex[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex[,] H = new Complex[ml * nl, mr * nr];
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
        public static Complex[,] Kronecker(this float[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            Complex[,] H = new Complex[ml * nl, mr * nr];
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
        public static Complex[,] Add(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        public static Complex[,] Add(this Complex[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        public static Complex[,] Add(this float[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(this Complex[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(this Complex[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(this float[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(Complex a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(Complex a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Add(float a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static Complex[,] Sub(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        public static Complex[,] Sub(this Complex[,] m, float[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        public static Complex[,] Sub(this float[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(this Complex[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(this Complex[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(this float[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(Complex a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(Complex a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Sub(float a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static Complex[,] Mul(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        public static Complex[,] Mul(this Complex[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        public static Complex[,] Mul(this float[,] m, Complex[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(this Complex[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(this Complex[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(this float[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(Complex a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(Complex a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Mul(float a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static Complex[,] Div(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        public static Complex[,] Div(this Complex[,] m, float[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        public static Complex[,] Div(this float[,] m, Complex[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            Complex[,] H = new Complex[mr, ml];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(this Complex[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(this Complex[,] m, float a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(this float[,] m, Complex a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(Complex a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(Complex a, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Div(float a, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="pow">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Pow(this Complex[,] m, float pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="pow">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Pow(this float[,] m, Complex pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="pow">Number</param>
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
                    H[i, j] = (float)Maths.Pow(m[i, j], pow);
                }
            }

            return H;
        }

        /// <summary>
        /// Raises the number to the power of the matrix.
        /// </summary>
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Pow(Complex a, float[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        /// <param name="a">Number</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Pow(float a, Complex[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static float[,] Negate(this float[,] m)
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
        public static Complex[,] Negate(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
        public static Complex[,] ToComplex(this float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
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
            float max = Matrice.Max(Matrice.Max(m));
            float min = Matrice.Min(Matrice.Min(m));
            float range = max - min;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = (m[i, j] - min) / range;
                }
            }
            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[,] Abs(this float[,] m)
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
        public static float[,] Abs(this Complex[,] m)
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
        public static float[,] Angle(this Complex[,] m)
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
        public static float[,] Real(this Complex[,] m)
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
        public static float[,] Imag(this Complex[,] m)
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
                    //if (Maths.IsSingular(m[j, i])) continue;
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
        public static Complex[] Sum(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[] kernel = new Complex[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    //if (Maths.IsSingular(m[j, i])) continue;
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
                    //if (Maths.IsSingular(m[j, i])) continue;
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
        public static Complex[] Mul(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[] kernel = new Complex[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    //if (Maths.IsSingular(m[j, i])) continue;
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
                    //if (Maths.IsSingular(m[j, i])) continue;
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
        public static Complex[] Div(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[] kernel = new Complex[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                kernel[i] = 1;

                for (j = 0; j < ml; j++)
                {
                    //if (Maths.IsSingular(m[j, i])) continue;
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
        public static Complex[] Mode(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[] h = new Complex[mr];
            Complex[] kernel = new Complex[ml];
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
        /// Returns the matrix vector corresponding to the specified threshold value.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="threshold">Threshold</param>
        /// <returns>Array</returns>
        public static float[] Morph(this float[,] m, int threshold)
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
                u[i] = v[threshold];
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
        public static Complex[] Mean(this Complex[,] m)
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
        public static Complex[] Var(this Complex[,] m)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            Complex[] u = Matrice.Mean(m);
            Complex[] v = new Complex[mr];
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
        public static Complex[] Var(this Complex[,] m, Complex[,] n)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            Complex[] v = new Complex[mr];
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
        public static Complex[] StnDev(this Complex[,] m)
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
        public static Complex[] StnDev(this Complex[,] m, Complex[,] n)
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
        public static Complex[,] Cov(this Complex[,] m)
        {
            Complex[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            Complex[,] H = (Complex[,])m.Clone();
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
            return Jagged.FromJagged(LinealgOptions.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this Complex[,] m, Complex[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this Complex[,] m, float[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this float[,] m, Complex[,] n)
        {
            return Jagged.FromJagged(LinealgOptions.Mul(Jagged.ToJagged(m), Jagged.ToJagged(n)));
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
            return LinealgOptions.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, Complex[,] n, bool normalize = true)
        {
            return LinealgOptions.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, float[,] n, bool normalize = true)
        {
            return LinealgOptions.Conv(m, n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this float[,] m, Complex[,] n, bool normalize = true)
        {
            return LinealgOptions.Conv(m, n, normalize);
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
                return LinealgOptions.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvVertical(LinealgOptions.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this float[,] m, Complex[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvVertical(LinealgOptions.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, float[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvVertical(LinealgOptions.ConvHorizontal(m, n, normalize), n, normalize);
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, Complex[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.ConvHorizontal(m, n, normalize);
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.ConvVertical(m, n, normalize);
            }

            // both processing
            return LinealgOptions.ConvVertical(LinealgOptions.ConvHorizontal(m, n, normalize), n, normalize);
        }
        #endregion

        #region Matrix morphology (separable)
        /// <summary>
        /// Returns the matrix result of morphological minimum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static float[,] Min(this float[,] m, int r0, int r1)
        {
            return LinealgOptions.MinVertical(LinealgOptions.MinHorizontal(m, r1), r0);
        }
        /// <summary>
        /// Returns the matrix result of morphological maximum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static float[,] Max(this float[,] m, int r0, int r1)
        {
            return LinealgOptions.MaxVertical(LinealgOptions.MaxHorizontal(m, r1), r0);
        }
        /// <summary>
        /// Returns the matrix result of morphology.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        /// <param name="threshold">Threshold</param>
        public static float[,] Morph(this float[,] m, int r0, int r1, int threshold)
        {
            return LinealgOptions.MorphVertical(LinealgOptions.MorphHorizontal(m, r1, threshold), r0, threshold);
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
            // both processing
            return LinealgOptions.MeanVertical(LinealgOptions.MeanHorizontal(m, r1), r0);
        }
        /// <summary>
        /// Returns the result matrix of local averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static Complex[,] Mean(this Complex[,] m, int r0, int r1)
        {
            return LinealgOptions.MeanVertical(LinealgOptions.MeanHorizontal(m, r1), r0);
        }
        #endregion

        // Vector voids

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
        public static bool IsEquals(this Complex[] a, Complex[] b)
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
        public static bool IsCollinear(this Complex[] a, Complex[] b)
        {
            int N = a.Length, i, j;
            Complex k;

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
        public static bool IsCollinear(this Complex[] a, float[] b)
        {
            int N = a.Length, i, j;
            Complex k;

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
        public static bool IsCollinear(this float[] a, Complex[] b)
        {
            int N = a.Length, i, j;
            Complex k;

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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float Norm(this float[] a)
        {
            return Norm(a, 2);
        }
        /// <summary>
        /// Returns the P-norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="p">Parameter p</param>
        /// <returns>float precision floating point number</returns>
        public static float Norm(this Complex[] a, float p)
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
        /// <returns>float precision floating point number</returns>
        public static float Norm(this Complex[] a)
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
        public static Complex[] Round(this Complex[] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0);
            Complex[] H = new Complex[ml];
            Complex c;
            int i;

            for (i = 0; i < ml; i++)
            {
                c = m[i];
                H[i] = new Complex((float)Math.Round(c.Real, digits, mode), (float)Math.Round(c.Imag, digits, mode));
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
        /// <returns>float precision floating point number</returns>
        public static float Angle(this float[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Angle(this Complex[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Angle(this float[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Angle(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }

        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static float Proj(this float[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Proj(this Complex[] a, float[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Proj(this float[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>float precision floating point number</returns>
        public static Complex Proj(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
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
        public static Complex[] Cosines(this Complex[] v)
        {
            int length = v.Length, i;
            Complex abs = Matrice.Norm(v);
            Complex[] cos = new Complex[length];

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
        public static Complex[] Add(this Complex[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Add(this Complex[] a, float[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Add(this float[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(this Complex[] a, Complex b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(this Complex[] a, float b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(this float[] a, Complex b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(Complex b, Complex[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(float b, Complex[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Add(Complex b, float[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Sub(this Complex[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Sub(this Complex[] a, float[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Sub(this float[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(this Complex[] a, Complex b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(this Complex[] a, float b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(this float[] a, Complex b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(Complex b, Complex[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(Complex b, float[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="b">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Sub(float b, Complex[] a)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Mul(this Complex[] a, float[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Mul(this float[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Mul(this Complex[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(this float[] v, Complex a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(this Complex[] v, float a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(this Complex[] v, Complex a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(Complex a, Complex[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(Complex a, float[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Mul(float a, Complex[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        public static Complex[] Div(this Complex[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Div(this Complex[] a, float[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        public static Complex[] Div(this float[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex[] c = new Complex[length];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(this float[] v, Complex a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(this Complex[] v, float a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(this Complex[] v, Complex a)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(Complex a, Complex[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(float a, Complex[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Div(Complex a, float[] v)
        {
            int length = v.Length, i;
            Complex[] H = new Complex[length];

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
        public static Complex[] Pow(this Complex[] v, float power)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

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
        public static Complex[] Pow(this float[] v, Complex power)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Maths.Pow(v[i], power);
            }
            return H;
        }

        /// <summary>
        /// Raises the number to the power of the vector.
        /// </summary>
        /// <param name="a">Number</param>
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
        /// <param name="a">Number</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Pow(Complex a, float[] v)
        {
            int n = v.GetLength(0);
            Complex[] H = new Complex[n];
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
        /// <param name="a">Number</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Pow(float a, Complex[] v)
        {
            int n = v.GetLength(0);
            Complex[] H = new Complex[n];
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
        public static float[] Abs(this float[] v)
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
        public static float[] Negate(this float[] v)
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
        public static Complex[] Negate(this Complex[] v)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

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
        public static Complex[] ToComplex(this float[] v)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

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
        public static float[] Abs(this Complex[] v)
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
        public static float[] Angle(this Complex[] v)
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
        public static float[] Real(this Complex[] v)
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
        public static float[] Imag(this Complex[] v)
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
        public static Complex[] Conjugate(this Complex[] v)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

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
        /// <returns>float precision floating point number</returns>
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
        public static Complex Sum(this Complex[] v)
        {
            Complex total = 0;
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static Complex Mul(this Complex[] v)
        {
            Complex total = 1;
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static Complex Div(this Complex[] v)
        {
            Complex total = 1;
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
        /// <returns>float precision floating point number</returns>
        public static float Mean(this float[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Returns the average value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex Mean(this Complex[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Returns the variance value.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>float precision floating point number</returns>
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
        public static Complex Var(this Complex[] v)
        {
            int length = v.Length;
            Complex mean = Matrice.Mean(v);
            Complex sum = 0;

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
        /// <returns>float precision floating point number</returns>
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
        public static Complex Var(this Complex[] x, Complex[] y)
        {
            int length = x.Length;
            Complex sum = 0;

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
        /// <returns>float precision floating point number</returns>
        public static float StnDev(this float[] v)
        {
            return (float)Math.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Complex number</returns>
        public static Complex StnDev(this Complex[] v)
        {
            return Maths.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Returns the standard deviation.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="y">Array</param>
        /// <returns>float precision floating point number</returns>
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
        public static Complex StnDev(this Complex[] x, Complex[] y)
        {
            return Maths.Sqrt(Matrice.Var(x, y));
        }
        /// <summary>
        /// Returns the value of the vector mode.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>float precision floating point number</returns>
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
        public static Complex Mode(this Complex[] v)
        {
            int count = 0;
            int length = v.Length;
            Complex frequent = 0;
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
        /// <returns>float precision floating point number</returns>
        public static float Min(this float[] v)
        {
            return Min(v, out _);
        }
        /// <summary>
        /// Gets the value of the minimum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="index">Max index</param>
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static float Max(this float[] v)
        {
            return Max(v, out _);
        }
        /// <summary>
        /// Gets the value of the maximum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="index">Max index</param>
        /// <returns>float precision floating point number</returns>
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
        /// Gets the value of the vector element corresponding to the threshold value.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="threshold">Threshold</param>
        /// <returns>float precision floating point number</returns>
        public static float Morph(this float[] v, int threshold)
        {
            float[] u = (float[])v.Clone();
            Array.Sort(u);
            return u[threshold];
        }
        /// <summary>
        /// Returns the covariance value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>float precision floating point number</returns>
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
        public static Complex Cov(this Complex[] v)
        {
            int xlength = v.Length;
            Complex xv = Matrice.Mean(v);
            Complex total = 0;
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
        /// <returns>float precision floating point number</returns>
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
        public static Complex[,] Dot(this Complex[,] m, Complex[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        public static Complex[,] Dot(this Complex[,] m, float[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        public static Complex[,] Dot(this float[,] m, Complex[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        public static Complex[,] Dot(this Complex[] v, Complex[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        public static Complex[,] Dot(this Complex[] v, float[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        public static Complex[,] Dot(this float[] v, Complex[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
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
        /// <returns>float precision floating point number</returns>
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
        /// <returns>float precision floating point number</returns>
        public static Complex Dot(this Complex[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex sum = 0;

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
        /// <returns>float precision floating point number</returns>
        public static Complex Dot(this Complex[] a, float[] b)
        {
            int length = a.Length, i;
            Complex sum = 0;

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
        /// <returns>float precision floating point number</returns>
        public static Complex Dot(this float[] a, Complex[] b)
        {
            int length = a.Length, i;
            Complex sum = 0;

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
        public static Complex[] Dot(this float[] v, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[] temp = new Complex[r0];

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
        public static Complex[] Dot(this Complex[] v, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[] temp = new Complex[r0];

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
        public static Complex[] Dot(this Complex[] v, float[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[] temp = new Complex[r0];

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
        public static Complex[,] Dotp(this Complex[] a, Complex[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            Complex c;
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
        public static Complex[,] Dotp(this Complex[] a, float[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            Complex c;
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
        public static Complex[,] Dotp(this float[] a, Complex[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            Complex c;
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
        public static Complex[] Conv(this Complex[] v, Complex[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex sum, div;
            Complex[] uv = new Complex[n];

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
        public static Complex[] Conv(this Complex[] v, float[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex sum, div;
            Complex[] uv = new Complex[n];

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
        public static Complex[] Conv(this float[] v, Complex[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            Complex sum, div;
            Complex[] uv = new Complex[n];

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
        /// Returns the vector result of morphology maximum.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Max(this float[] v, int r)
        {
            int l = v.Length;
            float[] B = new float[l];
            int r2 = r / 2;
            int i, j, k, c;
            float max;

            for (i = 0; i < l; i++)
            {
                max = float.MinValue;
                k = i - r2;

                for (j = 0; j < r; j++)
                {
                    c = k + j;

                    if (c < 0 || c >= l) continue;

                    if (v[c] > max)
                    {
                        max = v[c];
                    }
                }
                B[i] = max;
            }

            return B;
        }
        /// <summary>
        /// Returns the vector result of morphology minimum.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static float[] Min(this float[] v, int r)
        {
            int l = v.Length;
            float[] B = new float[l];
            int r2 = r / 2;
            int i, j, k, c;
            float min;

            for (i = 0; i < l; i++)
            {
                min = float.MaxValue;
                k = i - r2;

                for (j = 0; j < r; j++)
                {
                    c = k + j;

                    if (c < 0 || c >= l) continue;

                    if (v[c] < min)
                    {
                        min = v[c];
                    }
                }
                B[i] = min;
            }

            return B;
        }
        /// <summary>
        /// Returns the vector result of morphology.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        /// <param name="threshold">Threshold</param>
        /// <returns>Array</returns>
        public static float[] Morph(this float[] v, int r, int threshold)
        {
            int n = v.Length;
            float[] u = new float[r];
            float[] uv = new float[n];
            int r2 = r / 2;
            int i, j, k, c, p;

            for (i = 0; i < n; i++)
            {
                k = i - r2;
                p = 0;

                for (j = 0; j < r; j++)
                {
                    c = k + j;

                    if (c >= 0 && c < n)
                    {
                        u[p] = v[c]; p++;
                    }
                }

                Array.Sort(u);
                uv[i] = u[threshold];
            }

            return uv;
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
            int l = v.Length;
            if (l == 1)
                return v;

            float[] output = new float[l];
            int r2 = r >> 1;
            int dl = l - r2;
            float s = 0;
            int i;

            for (i = 0; i < r; i++)
            {
                s += v[i];
            }

            for (i = 0; i < r2; i++)
            {
                output[i] = s / r;
            }

            for (i = r2; i < dl; i++)
            {
                s = s - v[i - r2] + v[i + r2];
                output[i] = s / r;
            }

            for (i = dl; i < l; i++)
            {
                s = s - v[i - r2] + v[i];
                output[i] = s / r;
            }

            return output;
        }
        /// <summary>
        /// Returns the result vector of local averaging.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Radius</param>
        public static Complex[] Mean(this Complex[] v, int r)
        {
            int l = v.Length;
            if (l == 1)
                return v;

            Complex[] output = new Complex[l];
            int r2 = r >> 1;
            int dl = l - r2;
            Complex s = 0;
            int i;


            for (i = 0; i < r; i++)
            {
                s += v[i];
            }

            for (i = 0; i < r2; i++)
            {
                output[i] = s / r;
            }

            for (i = r2; i < dl; i++)
            {
                s = s - v[i - r2] + v[i + r2];
                output[i] = s / r;
            }

            for (i = dl; i < l; i++)
            {
                s = s - v[i - r2] + v[i];
                output[i] = s / r;
            }

            return output;
        }
        #endregion

        // MATLAB voids

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
        public static Complex[] GetCol(this Complex[,] m, int r)
        {
            int w = m.GetLength(0);
            Complex[] U = new Complex[w];
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
        public static Complex[,] SetCol(this Complex[,] m, Complex[] n, int r)
        {
            int w = m.GetLength(0);
            Complex[,] U = (Complex[,])m.Clone();
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
        public static Complex[] GetRow(this Complex[,] m, int r)
        {
            int w = m.GetLength(1);
            Complex[] U = new Complex[w];
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
        public static Complex[,] SetRow(this Complex[,] m, Complex[] n, int r)
        {
            int w = m.GetLength(1);
            Complex[,] U = (Complex[,])m.Clone();
            for (int i = 0; i < w; i++)
            {
                U[r, i] = n[i];
            }
            return U;
        }
        #endregion

        #region Shift voids
        /// <summary>
        /// Implements block matrix rearrangement.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="m">The number of positions to which a shift in height occurs</param>
        /// <param name="l">The number of positions by which the shift occurs in width</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Shift(this Complex[,] a, int m, int l)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            Complex[,] temp = new Complex[l0, l1];
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
        /// Implements block matrix rearrangement.
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
        public static Complex[] Shift(this Complex[] v, int l)
        {
            int N = v.Length;
            Complex[] temp = new Complex[N];

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
        public static Complex[,] Flip(this Complex[,] m, Direction direction)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
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
        public static Complex[] Flip(this Complex[] v)
        {
            int mr = v.Length;
            Complex[] H = new Complex[mr];

            for (int j = 0; j < mr; j++)
            {
                H[j] = v[mr - j - 1];
            }

            return H;
        }
        #endregion

        #region Merge voids
        /// <summary>
        /// Implements vector merging.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Merge(this float[] a, float[] b)
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
        /// Implements vector merging.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Merge(this Complex[] a, Complex[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex[] v = new Complex[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        /// <summary>
        /// Implements vector merging.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Merge(this Complex[] a, float[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex[] v = new Complex[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        /// <summary>
        /// Implements vector merging.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Merge(this float[] a, Complex[] b)
        {
            int na = a.Length, nb = b.Length, i;
            Complex[] v = new Complex[na + nb];

            for (i = 0; i < na; i++)
                v[i] = a[i];

            for (i = 0; i < nb; i++)
                v[na + i] = b[i];

            return v;
        }
        #endregion

        #region Cut voids
        /// <summary>
        /// Returns the specified part of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="start">Starting position</param>
        /// <param name="length">Vector length</param>
        /// <returns>Array</returns>
        public static float[] Cut(this float[] a, int start, int length)
        {
            int na = a.Length, i;
            float[] v = new float[length];

            for (i = 0; i < length; i++)
                v[i] = a[Maths.Mod(start + i, na)];

            return v;
        }
        /// <summary>
        /// Returns the specified part of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="start">Starting position</param>
        /// <param name="length">Vector length</param>
        /// <returns>Array</returns>
        public static Complex[] Cut(this Complex[] a, int start, int length)
        {
            int na = a.Length, i;
            Complex[] v = new Complex[length];

            for (i = 0; i < length; i++)
                v[i] = a[Maths.Mod(start + i, na)];

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
        /// <returns>Matrix</returns>
        public static float[,] Cut(this float[,] m, int y, int x, int height, int width)
        {
            // params
            float[,] B = new float[height, width];
            int i, j;

            // do job
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    B[i, j] = m[i + y, j + x];
                }
            }

            return B;
        }
        /// <summary>
        /// Crops the matrix to the specified size.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="y">Starting position in height</param>
        /// <param name="x">Starting position in width</param>
        /// <param name="height">Height</param>
        /// <param name="width">Width</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Cut(this Complex[,] m, int y, int x, int height, int width)
        {
            // params
            Complex[,] B = new Complex[height, width];
            int i, j;

            // do job
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    B[i, j] = m[i + y, j + x];
                }
            }

            return B;
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
        public static Complex[,] Reshape(this Complex[] a, int height)
        {
            int n = a.Length;
            int width = (int)Math.Ceiling(n / (float)height);
            Complex[,] H = new Complex[height, width];
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
        public static Complex[] Reshape(this Complex[,] a, int length)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            Complex[] v = new Complex[length];

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
        public static Complex[,] Diag(this Complex[] v)
        {
            int n = v.Length, i;
            Complex[,] diag = new Complex[n, n];

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
        public static Complex[] Diag(this Complex[,] a)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            Complex[] v = new Complex[height];
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
            // by rows and colums
            else
            {
                // rows:
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
                // colums:
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
            return;
        }
        /// <summary>
        /// Implements a permutation of the vectors of the matrix.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="i">First row or column</param>
        /// <param name="j">Second row or column</param>
        /// <param name="direction">Processing direction</param>
        public static void Swap(this Complex[,] a, int i, int j, Direction direction = Direction.Horizontal)
        {
            // properties:
            Complex temp;
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
            // by rows and colums
            else
            {
                // rows:
                for (z = 0; z < col; z++)
                {
                    temp = a[i, z];
                    a[i, z] = a[j, z];
                    a[j, z] = temp;
                }
                // colums:
                for (z = 0; z < row; z++)
                {
                    temp = a[z, i];
                    a[z, i] = a[z, j];
                    a[z, j] = temp;
                }
            }
            return;
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
            return;
        }
        /// <summary>
        /// Implements a permutation of the elements of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="i">First element position</param>
        /// <param name="j">Second element position</param>
        public static void Swap(this Complex[] v, int i, int j)
        {
            // get elements:
            Complex e1 = v[i], e2 = v[j];
            // swapping vector elements:
            v[j] = e1; v[i] = e2;
            return;
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

            // by colums:
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
        public static Complex[,] Remove(this Complex[,] a, int i, int length, Direction direction = Direction.Horizontal)
        {
            // properties:
            Complex[,] H;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z, y;

            // by colums:
            if (direction == Direction.Vertical)
            {
                H = new Complex[col, length];

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
                H = new Complex[length, row];

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
                H = new Complex[length, length];
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
            int n = v.Length;
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
        public static Complex[] Remove(this Complex[] v, int i, int length)
        {
            int n = v.Length;
            Complex[] w = new Complex[length];

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
            if (height != width) throw new Exception("The matrix must be square");
            if (n >= height || n < 0) throw new Exception("Row and column number specified is invalid");

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
        public static Complex[,] Minor(this Complex[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new Exception("The matrix must be square");
            if (n >= height || n < 0) throw new Exception("Row and column number specified is invalid");

            // new matrix:
            Complex[,] H = new Complex[height - 1, width - 1];
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
        public static Complex[] Extend(this Complex[] v, int length)
        {
            int r0 = v.GetLength(0);
            int rr = (length - r0) / 2;
            int dr = length - rr;
            Complex[] b = new Complex[length];
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
        public static Complex[,] Extend(this Complex[,] m, int height, int width)
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
        private static Complex[,] ExtendVertical(Complex[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int rr = (length - r0) / 2;
            Complex[,] B = new Complex[length, c0];
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
        private static Complex[,] ExtendHorizontal(Complex[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int cc = (length - c0) / 2;
            Complex[,] B = new Complex[r0, length];
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
        public static Complex[] Compute(this Complex[] v, IComplex function)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

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
        public static Complex[,] Compute(this float[] x, Complex[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex[,] z = new Complex[xlength, ylength];
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
        public static Complex[,] Compute(this Complex[] x, float[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex[,] z = new Complex[xlength, ylength];
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
        public static Complex[,] Compute(this Complex[] x, Complex[] y, IComplexMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            Complex[,] z = new Complex[xlength, ylength];
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
        public static Complex[,] Compute(this Complex[,] m, IComplex function)
        {
            int i, j;
            int ml = m.GetLength(1), mr = m.GetLength(0);
            Complex[,] H = new Complex[mr, ml];

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
        public static Complex[,] Companion(this Complex[] v)
        {
            int n = v.Length, i;
            Complex[,] H = new Complex[n, n];

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
        public static Complex[,] Vander(this Complex[] v)
        {
            int n = v.Length;
            Complex[,] H = new Complex[n, n];
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
        public static Complex[,] Hankeli(this Complex[] v)
        {
            int n = v.Length;
            Complex[,] H = new Complex[n, n];
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
        public static Complex[,] Hankel(this Complex[] v)
        {
            int n = v.Length / 2;
            Complex[,] H = new Complex[n, n];
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
        public static Complex[,] Toeplitz(this Complex[] v)
        {
            int n = v.Length / 2;
            Complex[,] H = new Complex[n, n];
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
        public static Complex[,] Cauchy(this Complex[] x, Complex[] y)
        {
            int m = x.Length, l = y.Length;
            Complex[,] H = new Complex[m, l];
            Complex kern;
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
        public static Complex[,] Circulant(this Complex[] v)
        {
            int n = v.Length;
            Complex[,] H = new Complex[n, n];
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
        public static Complex[,] Symmetric(this Complex[] v)
        {
            int n = v.Length;
            Complex[,] H = new Complex[n, n];
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
                throw new Exception("Dimension of the matrix must be an odd number");

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
        private static Random rnd = new Random();

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
        public static Complex[] Randc(int n)
        {
            Complex[] v = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = new Complex((float)rnd.NextDouble(), (float)rnd.NextDouble());
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
        public static Complex[,] Randc(int m, int l)
        {
            Complex[,] H = new Complex[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = new Complex((float)rnd.NextDouble(), (float)rnd.NextDouble());
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
        public static Complex[] Randic(int n)
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
        public static Complex[] Randic(int n, int a, int b)
        {
            Complex[] v = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = new Complex(rnd.Next(a, b), rnd.Next(a, b));
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
        public static Complex[,] Randic(int m, int l)
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
        public static Complex[,] Randic(int m, int l, int a, int b)
        {
            Complex[,] H = new Complex[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = new Complex(rnd.Next(a, b), rnd.Next(a, b));
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
            string[] nums = rows[0].Split('|');
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
                nums = rows[i].Split('|');
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
        public static Complex[,] Parse(this Complex[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length, n = nums.Length, k;
            Complex[,] H = new Complex[r, n];
            int i, j;

            // first row
            for (j = 0; j < n; j++)
            {
                H[0, j] = StringOptions.Compar(nums[j]);
            }

            // other rows
            for (i = 1; i < r; i++)
            {
                nums = rows[i].Split('|');
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
        public static bool TryParse(string s, out Complex[,] result)
        {
            Complex[,] zero = null;
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
                string[] nums = rows[0].Split('|');
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
                throw new Exception("The input string was in the wrong format");
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
        public static Complex[] Parse(this Complex[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
                string[] nums = rows[0].Split('|');
                int n = nums.Length, i;
                Complex[] H = new Complex[n];

                // collecting rows:
                for (i = 0; i < n; i++)
                {
                    H[i] = StringOptions.Compar(nums[i]);
                }
                return H;
            }
            else
            {
                throw new Exception("The input string was in the wrong format");
            }
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of complex numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, out Complex[] result)
        {
            Complex[] zero = null;
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
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height + 1 != width)
                throw new Exception("Input matrix has invalid sizes");

            float[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            float[] x = new float[height];
            float[] v, w;
            float temp;

            for (i = 0; i < height; i++)
            {
                w = B[i];
                temp = w[i];

                for (j = 0; j < width; j++)
                {
                    w[j] /= temp;
                }

                for (k = i + 1; k < height; k++)
                {
                    v = B[k];
                    temp = v[i];

                    for (j = i; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    B[k] = v;
                }

                B[i] = w;
            }

            for (i = 0; i < height; i++)
            {
                l = (height - 1) - i;
                w = B[l];

                for (k = 0; k < l; k++)
                {
                    v = B[k];
                    temp = v[l];

                    for (j = l; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    B[k] = v;
                }

                B[l] = w;
            }

            for (k = 0; k < height; k++)
            {
                x[k] = B[k][height];
            }
            return x;
        }
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static float[] Solve(this float[,] A, float[] b)
        {
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height != width)
                throw new Exception("The matrix must be square");
            if (height != b.Length)
                throw new Exception("Vector length should be equal to the height of the matrix");

            float[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            float[] x = (float[])b.Clone();
            float[] v, w;
            float temp;


            for (i = 0; i < height; i++)
            {
                w = B[i];
                temp = w[i];

                for (j = 0; j < width; j++)
                {
                    w[j] /= temp;
                }
                x[i] /= temp;

                for (k = i + 1; k < height; k++)
                {
                    v = B[k];
                    temp = v[i];

                    for (j = i; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    x[k] -= x[i] * temp;
                    B[k] = v;
                }
            }

            for (i = 0; i < height; i++)
            {
                l = (height - 1) - i;
                w = B[l];

                for (k = 0; k < l; k++)
                {
                    v = B[k];
                    temp = v[l];

                    for (j = l; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    x[k] -= x[l] * temp;
                    B[k] = v;
                }

                B[l] = w;
            }

            return x;
        }

        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Extended matrix</param>
        /// <returns>Array</returns>
        public static Complex[] Solve(this Complex[,] A)
        {
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height + 1 != width)
                throw new Exception("Input matrix has invalid sizes");

            Complex[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            Complex[] x = new Complex[height];
            Complex[] v, w;
            Complex temp;

            for (i = 0; i < height; i++)
            {
                w = B[i];
                temp = w[i];

                for (j = 0; j < width; j++)
                {
                    w[j] /= temp;
                }

                for (k = i + 1; k < height; k++)
                {
                    v = B[k];
                    temp = v[i];

                    for (j = i; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    B[k] = v;
                }

                B[i] = w;
            }

            for (i = 0; i < height; i++)
            {
                l = (height - 1) - i;
                w = B[l];

                for (k = 0; k < l; k++)
                {
                    v = B[k];
                    temp = v[l];

                    for (j = l; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    B[k] = v;
                }

                B[l] = w;
            }

            for (k = 0; k < height; k++)
            {
                x[k] = B[k][height];
            }
            return x;
        }
        /// <summary>
        /// Returns a vector corresponding to the solution of a system of linear algebraic equations: Ax = b.
        /// </summary>
        /// <param name="A">Square matrix</param>
        /// <param name="b">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Solve(this Complex[,] A, Complex[] b)
        {
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height != width)
                throw new Exception("The matrix must be square");
            if (height != b.Length)
                throw new Exception("Vector length should be equal to the height of the matrix");

            Complex[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            Complex[] x = (Complex[])b.Clone();
            Complex[] v, w;
            Complex temp;

            for (i = 0; i < height; i++)
            {
                w = B[i];
                temp = w[i];

                for (j = 0; j < width; j++)
                {
                    w[j] /= temp;
                }
                x[i] /= temp;

                for (k = i + 1; k < height; k++)
                {
                    v = B[k];
                    temp = v[i];

                    for (j = i; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    x[k] -= x[i] * temp;
                    B[k] = v;
                }
            }

            for (i = 0; i < height; i++)
            {
                l = (height - 1) - i;
                w = B[l];

                for (k = 0; k < l; k++)
                {
                    v = B[k];
                    temp = v[l];

                    for (j = l; j < width; j++)
                    {
                        v[j] = v[j] - w[j] * temp;
                    }

                    x[k] -= x[l] * temp;
                    B[k] = v;
                }

                B[l] = w;
            }

            return x;
        }
        #endregion
    }
}
