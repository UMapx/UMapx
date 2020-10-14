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
        public static bool IsEquals(this double[,] m, double[,] n)
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
        public static bool IsVector(this double[,] m)
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
        public static bool IsSquare(this double[,] m)
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
        public static bool IsPositive(this double[,] m)
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
        public static bool IsSymmetric(this double[,] m)
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
        public static bool IsSkewSymmetric(this double[,] m)
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
        public static bool IsDiagonal(this double[,] m)
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
        public static double[,] Invert(this double[,] m)
        {
            int n = m.GetLength(0);

            if (n != m.GetLength(1))
                throw new Exception("The matrix must be square");

            double[][] z = m.ToJagged();
            int n2 = n * 2;
            double[] w;
            double[] v;
            double[][] a = new double[n][];
            int i, j;

            for (i = 0; i < n; i++)
            {
                w = z[i];
                v = new double[n2];

                for (j = 0; j < n; j++)
                {
                    v[j] = w[j];
                }

                a[i] = v;
                a[i][i + n] = 1;
            }

            const double epsilon = 1e-4;
            int k, c, l, t;
            double temp, factor, div1, div2;

            for (i = 0; i < n; i++)
            {
                v = a[i];

                // first decomposition:
                for (k = i + 1; k < n; k++)
                {
                    w = a[k];

                    if (Math.Abs(w[i]) > epsilon)
                    {
                        for (c = 0; c < n2; c++)
                        {
                            temp = v[c];
                            v[c] = w[c];
                            w[c] = temp;
                        }
                        break;
                    }
                }
                {
                    // second decomposition:
                    div1 = v[i];

                    for (j = 0; j < n2; j++)
                    {
                        if (j != i)
                        {
                            v[j] /= div1;
                        }
                    }

                    v[i] = 1;
                    div2 = v[i];

                    for (t = 0; t < n; t++)
                    {
                        if (t != i)
                        {
                            w = a[t];
                            factor = w[i] / div2;

                            for (l = 0; l < n2; l++)
                            {
                                w[l] -= factor * v[l];
                            }
                        }
                    }
                }
            }

            // building invert matrix:
            double[][] inv = new double[n][];

            for (i = 0; i < n; i++)
            {
                w = new double[n];
                v = a[i];

                for (j = 0; j < n; j++)
                {
                    w[j] = v[j + n];
                }

                inv[i] = w;
            }

            return inv.FromJagged();
        }
        /// <summary>
        /// Implements the transpose of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[,] Transponate(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r1, r0];
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
            int n = m.GetLength(0);

            if (n != m.GetLength(1))
                throw new Exception("The matrix must be square");

            Complex[][] z = m.ToJagged();
            int n2 = n * 2;
            Complex[] w;
            Complex[] v;
            Complex[][] a = new Complex[n][];
            int i, j;

            for (i = 0; i < n; i++)
            {
                w = z[i];
                v = new Complex[n2];

                for (j = 0; j < n; j++)
                {
                    v[j] = w[j];
                }

                a[i] = v;
                a[i][i + n] = 1;
            }

            const double epsilon = 1e-4;
            int k, c, l, t;
            Complex temp, factor, div1, div2;

            for (i = 0; i < n; i++)
            {
                v = a[i];

                // first decomposition:
                for (k = i + 1; k < n; k++)
                {
                    w = a[k];

                    if (Maths.Abs(w[i]) > epsilon)
                    {
                        for (c = 0; c < n2; c++)
                        {
                            temp = v[c];
                            v[c] = w[c];
                            w[c] = temp;
                        }
                        break;
                    }
                }
                {
                    // second decomposition:
                    div1 = v[i];

                    for (j = 0; j < n2; j++)
                    {
                        if (j != i)
                        {
                            v[j] /= div1;
                        }
                    }

                    v[i] = 1;
                    div2 = v[i];

                    for (t = 0; t < n; t++)
                    {
                        if (t != i)
                        {
                            w = a[t];
                            factor = w[i] / div2;

                            for (l = 0; l < n2; l++)
                            {
                                w[l] -= factor * v[l];
                            }
                        }
                    }
                }
            }

            // building invert matrix:
            Complex[][] inv = new Complex[n][];

            for (i = 0; i < n; i++)
            {
                w = new Complex[n];
                v = a[i];

                for (j = 0; j < n; j++)
                {
                    w[j] = v[j + n];
                }

                inv[i] = w;
            }

            return inv.FromJagged();
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
        /// <returns>Double precision floating point number</returns>
        public static double Trace(this double[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("The matrix must be square");

            int d = m.GetLength(0);
            int i;
            double kernel = 0;

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
        /// <returns>Double precision floating point number</returns>
        public static double Det(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new Exception("The matrix must be square");

            unsafe
            {
                // copy array
                double[,] n = (double[,])m.Clone();

                fixed (double* pm = &n[0, 0])
                    return LinealgOptions.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Returns the P-norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Matrix</returns>
        public static double Norm(this double[,] m, double p)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double norm = 0;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    norm += Math.Pow(Math.Abs(m[i, j]), p);
                }
            }

            return Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double Norm(this double[,] m)
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
        public static double[,] Round(this double[,] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = Math.Round(m[i, j], digits, mode);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a square permutation matrix.
        /// </summary>
        /// <param name="m">Square matrix</param>
        /// <returns>Square matrix</returns>
        public static double[,] Permutation(this double[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("The matrix must be square");

            int i, j, r, n = m.GetLength(0);
            double[] temp; double diagonal;
            double[][] perm = Jagged.ToJagged(Matrice.Eye(n, n));

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
        /// <returns>Double precision floating point number</returns>
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
        public static double Norm(this Complex[,] m, double p)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double norm = 0;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    norm += Math.Pow(m[i, j].Abs, p);
                }
            }

            return Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double Norm(this Complex[,] m)
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
                    H[i, j] = new Complex(Math.Round(c.Real, digits, mode), Math.Round(c.Imag, digits, mode));
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
        public static double[,] Kronecker(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            double[,] H = new double[ml * nl, mr * nr];
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
        public static Complex[,] Kronecker(this Complex[,] m, double[,] n)
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
        public static Complex[,] Kronecker(this double[,] m, Complex[,] n)
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
        public static double[,] Add(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[,] H = new double[ml, mr];
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
        public static Complex[,] Add(this Complex[,] m, double[,] n)
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
        public static Complex[,] Add(this double[,] m, Complex[,] n)
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
        public static double[,] Add(this double[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Add(this Complex[,] m, double a)
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
        public static Complex[,] Add(this double[,] m, Complex a)
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
        public static double[,] Add(double a, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Add(Complex a, double[,] m)
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
        public static Complex[,] Add(double a, Complex[,] m)
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
        public static double[,] Sub(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[,] H = new double[ml, mr];
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
        public static Complex[,] Sub(this Complex[,] m, double[,] n)
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
        public static Complex[,] Sub(this double[,] m, Complex[,] n)
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
        public static double[,] Sub(this double[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Sub(this Complex[,] m, double a)
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
        public static Complex[,] Sub(this double[,] m, Complex a)
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
        public static double[,] Sub(double a, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Sub(Complex a, double[,] m)
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
        public static Complex[,] Sub(double a, Complex[,] m)
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
        public static double[,] Mul(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            double[,] H = new double[mr, ml];

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
        public static Complex[,] Mul(this Complex[,] m, double[,] n)
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
        public static Complex[,] Mul(this double[,] m, Complex[,] n)
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
        public static double[,] Mul(this double[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Mul(this Complex[,] m, double a)
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
        public static Complex[,] Mul(this double[,] m, Complex a)
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
        public static double[,] Mul(double a, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Mul(Complex a, double[,] m)
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
        public static Complex[,] Mul(double a, Complex[,] m)
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
        public static double[,] Div(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(1), mr = m.GetLength(0);
            int i, j;
            double[,] H = new double[mr, ml];

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
        public static Complex[,] Div(this Complex[,] m, double[,] n)
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
        public static Complex[,] Div(this double[,] m, Complex[,] n)
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
        public static double[,] Div(this double[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Div(this Complex[,] m, double a)
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
        public static Complex[,] Div(this double[,] m, Complex a)
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
        public static double[,] Div(double a, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] Div(Complex a, double[,] m)
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
        public static Complex[,] Div(double a, Complex[,] m)
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
        public static Complex[,] Pow(this Complex[,] m, double pow)
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
        public static Complex[,] Pow(this double[,] m, Complex pow)
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
        public static double[,] Pow(this double[,] m, double pow)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        /// Raises the number to the power of the matrix.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[,] Pow(double a, double[,] m)
        {
            int r0 = m.GetLength(0);
            int r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Math.Pow(a, m[i, j]);
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
        public static Complex[,] Pow(Complex a, double[,] m)
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
        public static Complex[,] Pow(double a, Complex[,] m)
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

        #region Matrix log and exp
        /// <summary>
        /// Logarithms all elements of the matrix at the base.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static double[,] Log(this double[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Log(m[i, j], a);
                }
            }

            return H;
        }
        /// <summary>
        /// Takes an exponent from all matrix values.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[,] Exp(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Exp(m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Logarithms all elements of the matrix at the base.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="a">Number</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Log(this Complex[,] m, double a)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Log(m[i, j], a);
                }
            }

            return H;
        }
        /// <summary>
        /// Takes an exponent from all matrix values.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Exp(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Exp(m[i, j]);
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
        public static double[,] Negate(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static Complex[,] ToComplex(this double[,] m)
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
        public static double[,] ToByte(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static double[,] ToDouble(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            double max = Matrice.Max(Matrice.Max(m));
            double min = Matrice.Min(Matrice.Min(m));
            double range = max - min;
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
        public static double[,] Abs(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = Maths.Abs(m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[,] Abs(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static double[,] Angle(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static double[,] Real(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        public static double[,] Imag(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
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
        /// Sorts the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        public static void Sort(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int i, j, k;
            double temp;

            for (k = 0; k < mr; k++)
            {
                for (i = 0; i < ml; i++)
                {
                    for (j = i + 1; j < ml; j++)
                    {
                        if (m[i, k] > m[j, k])
                        {
                            temp = m[i, k];
                            m[i, k] = m[j, k];
                            m[j, k] = temp;
                        }
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Returns the vector of matrix sums.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static double[] Sum(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] kernel = new double[mr];
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
        public static double[] Mul(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] kernel = new double[mr];
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
        public static double[] Div(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] kernel = new double[mr];
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
        /// Returns the maximum matrix vector.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static double[] Max(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                v[i] = double.MinValue;
            }

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (m[j, i] > v[i])
                    {
                        v[i] = m[j, i];
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
        public static double[] Min(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                v[i] = double.MaxValue;
            }

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (m[j, i] < v[i])
                    {
                        v[i] = m[j, i];
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
        public static double[] Morph(this double[,] m, int threshold)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[ml];
            double[] u = new double[mr];
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
        public static double[] Mean(this double[,] m)
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
        public static double[] Var(this double[,] m)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            double[] u = Matrice.Mean(m);
            double[] v = new double[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += Math.Pow(m[j, i] - u[i], 2);
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
        public static double[] Var(this double[,] m, double[,] n)
        {
            int mr = m.GetLength(1), ml = m.GetLength(0);
            double[] v = new double[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[i] += Math.Pow(m[j, i] - n[j, i], 2);
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
        public static double[] StnDev(this double[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Array</returns>
        public static Complex[] StnDev(this Complex[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static double[] StnDev(this double[,] m, double[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5);
        }
        /// <summary>
        /// Returns the standard deviation vector of the matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Array</returns>
        public static Complex[] StnDev(this Complex[,] m, Complex[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5);
        }
        /// <summary>
        /// Returns the covariance matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[,] Cov(this double[,] m)
        {
            double[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            double[,] H = (double[,])m.Clone();
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
        public static double[] Entropy(this double[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            double[] v = new double[width];
            int i, j;

            for (i = 0; i < width; i++)
            {
                for (j = 0; j < height; j++)
                {
                    if (m[j, i] > 0)
                    {
                        v[i] += -m[j, i] * Maths.Log2(m[j, i]);
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
        public static double[,] Dot(this double[,] m, double[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this Complex[,] m, Complex[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this Complex[,] m, double[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Implements a scalar product of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Dot(this double[,] m, Complex[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
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
        public static double[,] Conv(this double[,] m, double[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
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
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, double[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this double[,] m, Complex[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
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
        public static double[,] Conv(this double[,] m, double[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvHorizontal(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }

            // both processing
            float[] nn = LinealgOptions.ToJagged(n);
            float[][] mm = LinealgOptions.ToJagged(m);

            return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                                             LinealgOptions.ConvHorizontal(mm, nn, normalize), nn, normalize));
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this double[,] m, Complex[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvHorizontal(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }

            // both processing
            LinealgOptions.Complex32[] nn = LinealgOptions.ToJagged(n);
            float[][] mm = LinealgOptions.ToJagged(m);

            return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                                             LinealgOptions.ConvHorizontal(mm, nn, normalize), nn, normalize));
        }
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="n">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] Conv(this Complex[,] m, double[] n, Direction direction, bool normalize = true)
        {
            // direction of processing
            if (direction == Direction.Horizontal)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvHorizontal(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }

            // both processing
            float[] nn = LinealgOptions.ToJagged(n);
            LinealgOptions.Complex32[][] mm = LinealgOptions.ToJagged(m);

            return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                                             LinealgOptions.ConvHorizontal(mm, nn, normalize), nn, normalize));
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
                return LinealgOptions.FromJagged(LinealgOptions.ConvHorizontal(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }
            else if (direction == Direction.Vertical)
            {
                return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                    LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
            }

            // both processing
            LinealgOptions.Complex32[] nn = LinealgOptions.ToJagged(n);
            LinealgOptions.Complex32[][] mm = LinealgOptions.ToJagged(m);

            return LinealgOptions.FromJagged(LinealgOptions.ConvVertical(
                                             LinealgOptions.ConvHorizontal(mm, nn, normalize), nn, normalize));
        }
        #endregion

        #region Matrix morphology (separable)
        /// <summary>
        /// Returns the matrix result of morphological minimum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static double[,] Min(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MinVertical(
                                             LinealgOptions.MinHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Returns the matrix result of morphological maximum.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static double[,] Max(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MaxVertical(
                                             LinealgOptions.MaxHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Returns the matrix result of morphology.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        /// <param name="threshold">Threshold</param>
        public static double[,] Morph(this double[,] m, int r0, int r1, int threshold)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MorphVertical(
                                             LinealgOptions.MorphHorizontal(mm, r1, threshold), r0, threshold));
        }
        #endregion

        #region Matrix mean (separable)
        /// <summary>
        /// Returns the result matrix of local averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static double[,] Mean(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MeanVertical(
                                             LinealgOptions.MeanHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Returns the result matrix of local averaging.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="r0">Height radius</param>
        /// <param name="r1">Width radius</param>
        public static Complex[,] Mean(this Complex[,] m, int r0, int r1)
        {
            // both processing
            LinealgOptions.Complex32[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MeanVertical(
                                             LinealgOptions.MeanHorizontal(mm, r1), r0));
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
        public static bool IsEquals(this double[] a, double[] b)
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
        public static bool IsPositive(this double[] v)
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
        public static bool IsCollinear(this double[] a, double[] b)
        {
            int N = a.Length, i, j;
            double k;

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
        public static bool IsCollinear(this Complex[] a, double[] b)
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
        public static bool IsCollinear(this double[] a, Complex[] b)
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
        /// <returns>Double precision floating point number</returns>
        public static double Norm(this double[] a, double p)
        {
            int length = a.Length, i;
            double norm = 0;

            for (i = 0; i < length; i++)
            {
                norm += Math.Pow(Math.Abs(a[i]), p);
            }
            return Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Norm(this double[] a)
        {
            return Norm(a, 2);
        }
        /// <summary>
        /// Returns the P-norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="p">Parameter p</param>
        /// <returns>Double precision floating point number</returns>
        public static double Norm(this Complex[] a, double p)
        {
            int length = a.Length, i;
            double norm = 0;

            for (i = 0; i < length; i++)
            {
                norm += Maths.Pow(Maths.Abs(a[i]), p);
            }
            return Maths.Sqrt(norm, p);
        }
        /// <summary>
        /// Returns the norm of the vector.
        /// </summary>
        /// <param name="a">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Norm(this Complex[] a)
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
        public static double[] Round(this double[] m, int digits, MidpointRounding mode)
        {
            int ml = m.GetLength(0);
            double[] H = new double[ml];
            int i;

            for (i = 0; i < ml; i++)
            {
                H[i] = Math.Round(H[i], digits, mode);
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
                H[i] = new Complex(Math.Round(c.Real, digits, mode), Math.Round(c.Imag, digits, mode));
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
        /// <returns>Double precision floating point number</returns>
        public static double Angle(this double[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Angle(this Complex[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Angle(this double[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the angle between two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Angle(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }

        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Proj(this double[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Proj(this Complex[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Proj(this double[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Returns the projection of two vectors.
        /// </summary>
        /// <param name="a">Array</param>
        /// <param name="b">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static Complex Proj(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }

        /// <summary>
        /// Returns the direction cosines of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static double[] Cosines(this double[] v)
        {
            int length = v.Length, i;
            double abs = Matrice.Norm(v);
            double[] cos = new double[length];

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
        public static double[] Add(this double[] a, double[] b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Add(this Complex[] a, double[] b)
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
        public static Complex[] Add(this double[] a, Complex[] b)
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
        public static double[] Add(this double[] a, double b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Add(this Complex[] a, double b)
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
        public static Complex[] Add(this double[] a, Complex b)
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
        public static double[] Add(double b, double[] a)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Add(double b, Complex[] a)
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
        public static Complex[] Add(Complex b, double[] a)
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
        public static double[] Sub(this double[] a, double[] b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Sub(this Complex[] a, double[] b)
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
        public static Complex[] Sub(this double[] a, Complex[] b)
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
        public static double[] Sub(this double[] a, double b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Sub(this Complex[] a, double b)
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
        public static Complex[] Sub(this double[] a, Complex b)
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
        public static double[] Sub(double b, double[] a)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Sub(Complex b, double[] a)
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
        public static Complex[] Sub(double b, Complex[] a)
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
        public static double[] Mul(this double[] a, double[] b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Mul(this Complex[] a, double[] b)
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
        public static Complex[] Mul(this double[] a, Complex[] b)
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
        public static double[] Mul(this double[] v, double a)
        {
            int length = v.Length, i;
            double[] H = new double[length];

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
        public static Complex[] Mul(this double[] v, Complex a)
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
        public static Complex[] Mul(this Complex[] v, double a)
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
        public static double[] Mul(double a, double[] v)
        {
            int length = v.Length, i;
            double[] H = new double[length];

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
        public static Complex[] Mul(Complex a, double[] v)
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
        public static Complex[] Mul(double a, Complex[] v)
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
        public static double[] Div(this double[] a, double[] b)
        {
            int length = a.Length, i;
            double[] c = new double[length];

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
        public static Complex[] Div(this Complex[] a, double[] b)
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
        public static Complex[] Div(this double[] a, Complex[] b)
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
        public static double[] Div(this double[] v, double a)
        {
            int length = v.Length, i;
            double[] H = new double[length];

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
        public static Complex[] Div(this double[] v, Complex a)
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
        public static Complex[] Div(this Complex[] v, double a)
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
        public static double[] Div(double a, double[] v)
        {
            int length = v.Length, i;
            double[] H = new double[length];

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
        public static Complex[] Div(double a, Complex[] v)
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
        public static Complex[] Div(Complex a, double[] v)
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
        public static double[] Pow(this double[] v, double power)
        {
            int length = v.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Math.Pow(v[i], power);
            }
            return H;
        }
        /// <summary>
        /// Raises the elements of a vector to a power.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="power">Power</param>
        /// <returns>Array</returns>
        public static Complex[] Pow(this Complex[] v, double power)
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
        public static Complex[] Pow(this double[] v, Complex power)
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
        public static double[] Pow(double a, double[] v)
        {
            int n = v.GetLength(0);
            double[] H = new double[n];
            int i;

            for (i = 0; i < n; i++)
            {
                H[i] = Math.Pow(a, v[i]);
            }

            return H;
        }
        /// <summary>
        /// Raises the number to the power of the vector.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Pow(Complex a, double[] v)
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
        public static Complex[] Pow(double a, Complex[] v)
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

        #region Vector exp/log
        /// <summary>
        /// Logarithms the elements of the vector base.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static double[] Log(this double[] v, double a = Math.E)
        {
            int length = v.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Math.Log(v[i], a);
            }
            return H;
        }
        /// <summary>
        /// Logarithms the elements of the vector base.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="a">Number</param>
        /// <returns>Array</returns>
        public static Complex[] Log(this Complex[] v, double a = Math.E)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Maths.Log(v[i], a);
            }
            return H;
        }
        /// <summary>
        /// Takes an exponent from all vector values.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static double[] Exp(this double[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Math.Exp(v[i]);
            }
            return H;
        }
        /// <summary>
        /// Takes an exponent from all vector values.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static Complex[] Exp(this Complex[] v)
        {
            int length = v.Length;
            Complex[] H = new Complex[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Maths.Exp(v[i]);
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
        public static double[] ToDouble(this double[] v)
        {
            int length = v.Length;
            double max = Matrice.Max(v);
            double min = Matrice.Min(v);
            double range = max - min;
            double[] u = new double[length];

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
        public static double[] ToByte(this double[] v)
        {
            int length = v.Length;
            double[] u = new double[length];

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
        public static double[] Abs(this double[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static double[] Negate(this double[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static Complex[] ToComplex(this double[] v)
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
        public static double[] Abs(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static double[] Angle(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static double[] Real(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static double[] Imag(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        /// Returns the covariance value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cov(this double[] v)
        {
            int xlength = v.Length;
            double xv = Matrice.Mean(v);
            double total = 0;
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
        /// Returns the total value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sum(this double[] v)
        {
            double total = 0;
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
        /// <returns>Double precision floating point number</returns>
        public static double Mul(this double[] v)
        {
            double total = 1;
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
        /// <returns>Double precision floating point number</returns>
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
        /// <returns>Double precision floating point number</returns>
        public static double Div(this double[] v)
        {
            double total = 1;
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
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the entropy of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Entropy(this double[] v)
        {
            double H = 0;
            int length = v.Length, i;

            for (i = 0; i < length; i++)
            {
                if (v[i] > 0)
                {
                    H += -v[i] * Maths.Log2(v[i]);
                }
            }
            return H;
        }
        /// <summary>
        /// Returns the average value of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Mean(this double[] v)
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
        /// <returns>Double precision floating point number</returns>
        public static double Var(this double[] v)
        {
            int length = v.Length;
            double mean = Matrice.Mean(v);
            double sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += Math.Pow(v[i] - mean, 2);
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
        /// <returns>Double precision floating point number</returns>
        public static double Var(this double[] x, double[] y)
        {
            int length = x.Length;
            double sum = 0;

            for (int i = 0; i < length; i++)
            {
                sum += Math.Pow(x[i] - y[i], 2);
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
        /// <returns>Double precision floating point number</returns>
        public static double StnDev(this double[] v)
        {
            return Math.Sqrt(Matrice.Var(v));
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
        /// <returns>Double precision floating point number</returns>
        public static double StnDev(this double[] x, double[] y)
        {
            return Math.Sqrt(Matrice.Var(x, y));
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
        /// <returns>Double precision floating point number</returns>
        public static double Mode(this double[] v)
        {
            int count = 0;
            int length = v.Length;
            double frequent = 0;
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
        /// Returns the minimum and maximum values of a vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Pair of numbers</returns>
        public static RangeDouble Extremum(this double[] v)
        {
            double max = double.MinValue, min = double.MaxValue;
            int length = v.Length;
            double c;

            for (int i = 0; i < length; i++)
            {
                c = v[i];

                if (c > max)
                    max = c;
                if (c < min)
                    min = c;
            }
            return new RangeDouble(min, max);
        }
        /// <summary>
        /// Gets the value of the minimum element of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Double precision floating point number</returns>
        public static double Min(this double[] v)
        {
            int index = 0, length = v.Length;
            double minimum = double.MaxValue;
            double c;

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
        /// <returns>Double precision floating point number</returns>
        public static double Max(this double[] v)
        {
            int index = 0, length = v.Length;
            double maximum = double.MinValue;
            double c;

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
        /// <returns>Double precision floating point number</returns>
        public static double Morph(this double[] v, int threshold)
        {
            double[] u = (double[])v.Clone();
            Array.Sort(u);
            return u[threshold];
        }
        /// <summary>
        /// Sorts the vector.
        /// </summary>
        /// <param name="v">Array</param>
        public static void Sort(this double[] v)
        {
            int i, j, length = v.Length;
            double temp;
            for (i = 0; i < length; i++)
            {
                for (j = i + 1; j < length; j++)
                {
                    if (v[i] > v[j])
                    {
                        temp = v[i];
                        v[i] = v[j];
                        v[j] = temp;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Sorts the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="r">Start</param>
        /// <param name="l">End</param>
        public static void Sort(this double[] v, int r, int l)
        {
            double temp;
            double x = v[r + (l - r) / 2];
            int i = r;
            int j = l;

            while (i <= j)
            {
                while (v[i] < x) i++;
                while (v[j] > x) j--;
                if (i <= j)
                {
                    temp = v[i];
                    v[i] = v[j];
                    v[j] = temp;
                    i++;
                    j--;
                }
            }
            if (i < l)
                Sort(v, i, l);

            if (r < j)
                Sort(v, r, j);
            return;
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
        public static double[,] Dot(this double[,] m, double[] v, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] temp = new double[r0, r1];
            double alpha;
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
        public static Complex[,] Dot(this Complex[,] m, double[] v, bool inverse = false)
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
        public static Complex[,] Dot(this double[,] m, Complex[] v, bool inverse = false)
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
        public static double[,] Dot(this double[] v, double[,] m, bool inverse = false)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] temp = new double[r0, r1];
            double alpha;
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
        public static Complex[,] Dot(this Complex[] v, double[,] m, bool inverse = false)
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
        public static Complex[,] Dot(this double[] v, Complex[,] m, bool inverse = false)
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
        /// <returns>Double precision floating point number</returns>
        public static double Dot(this double[] a, double[] b)
        {
            int length = a.Length, i;
            double sum = 0;

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
        /// <returns>Double precision floating point number</returns>
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
        /// <returns>Double precision floating point number</returns>
        public static Complex Dot(this Complex[] a, double[] b)
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
        /// <returns>Double precision floating point number</returns>
        public static Complex Dot(this double[] a, Complex[] b)
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
        public static double[] Dot(this double[] v, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[] temp = new double[r0];

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
        public static Complex[] Dot(this double[] v, Complex[,] m)
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
        public static Complex[] Dot(this Complex[] v, double[,] m)
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
        public static double[,] Dotp(this double[] a, double[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            double[,] H = new double[l0, l1];
            double c;
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
        public static Complex[,] Dotp(this Complex[] a, double[] b)
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
        public static Complex[,] Dotp(this double[] a, Complex[] b)
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
        public static double[] Conv(this double[] v, double[] u, bool normalize = true)
        {
            int n = v.Length;
            int r = u.Length;
            int r2 = r / 2;
            int i, j, k, c;
            double sum, div;
            double[] uv = new double[n];

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
        public static Complex[] Conv(this Complex[] v, double[] u, bool normalize = true)
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
        public static Complex[] Conv(this double[] v, Complex[] u, bool normalize = true)
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
        public static double[] Max(this double[] v, int r)
        {
            int l = v.Length;
            double[] B = new double[l];
            int r2 = r / 2;
            int i, j, k, c;
            double max;

            for (i = 0; i < l; i++)
            {
                max = double.MinValue;
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
        public static double[] Min(this double[] v, int r)
        {
            int l = v.Length;
            double[] B = new double[l];
            int r2 = r / 2;
            int i, j, k, c;
            double min;

            for (i = 0; i < l; i++)
            {
                min = double.MaxValue;
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
        public static double[] Morph(this double[] v, int r, int threshold)
        {
            int n = v.Length;
            double[] u = new double[r];
            double[] uv = new double[n];
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
        public static double[] Mean(this double[] v, int r)
        {
            int l = v.Length;
            double[] output = new double[l];
            int r2 = r >> 1;
            int dl = l - r2;
            double s = 0;
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
        public static double[] GetCol(this double[,] m, int r)
        {
            int w = m.GetLength(0);
            double[] U = new double[w];
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
        public static double[,] SetCol(this double[,] m, double[] n, int r)
        {
            int w = m.GetLength(0);
            double[,] U = (double[,])m.Clone();
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
        public static double[] GetRow(this double[,] m, int r)
        {
            int w = m.GetLength(1);
            double[] U = new double[w];
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
        public static double[,] SetRow(this double[,] m, double[] n, int r)
        {
            int w = m.GetLength(1);
            double[,] U = (double[,])m.Clone();
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
        public static double[,] Shift(this double[,] a, int m, int l)
        {
            int l0 = a.GetLength(0), l1 = a.GetLength(1);
            double[,] temp = new double[l0, l1];
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
        public static double[] Shift(this double[] v, int l)
        {
            int N = v.Length;
            double[] temp = new double[N];

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
        public static double[,] Flip(this double[,] m, Direction direction)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[,] H = new double[ml, mr];
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
        public static double[] Flip(this double[] v)
        {
            int mr = v.Length;
            double[] H = new double[mr];

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
        public static double[] Merge(this double[] a, double[] b)
        {
            int na = a.Length, nb = b.Length, i;
            double[] v = new double[na + nb];

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
        public static Complex[] Merge(this Complex[] a, double[] b)
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
        public static Complex[] Merge(this double[] a, Complex[] b)
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
        public static double[] Cut(this double[] a, int start, int length)
        {
            int na = a.Length, i;
            double[] v = new double[length];

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
        public static double[,] Cut(this double[,] m, int y, int x, int height, int width)
        {
            // params
            double[,] B = new double[height, width];
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
        public static double[,] Reshape(this double[] a, int height)
        {
            int n = a.Length;
            int width = (int)Math.Ceiling(n / (double)height);
            double[,] H = new double[height, width];
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
            int width = (int)Math.Ceiling(n / (double)height);
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
        public static double[] Reshape(this double[,] a, int length)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            double[] v = new double[length];

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
        public static double[,] Diag(this double[] v)
        {
            int n = v.Length, i;
            double[,] diag = new double[n, n];

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
        public static double[] Diag(this double[,] a)
        {
            int height = a.GetLength(0), width = a.GetLength(1);
            double[] v = new double[height];
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
        public static void Swap(this double[,] a, int i, int j, Direction direction = Direction.Horizontal)
        {
            // properties:
            double temp;
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
        public static void Swap(this double[] v, int i, int j)
        {
            // get elements:
            double e1 = v[i], e2 = v[j];
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
        public static double[,] Remove(this double[,] a, int i, int length, Direction direction = Direction.Horizontal)
        {
            // properties:
            double[,] H;
            int col = a.GetLength(0);
            int row = a.GetLength(1);
            int z, y;

            // by colums:
            if (direction == Direction.Vertical)
            {
                H = new double[col, length];

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
                H = new double[length, row];

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
                H = new double[length, length];
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
        public static double[] Remove(this double[] v, int i, int length)
        {
            int n = v.Length;
            double[] w = new double[length];

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
        public static double[,] Minor(this double[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new Exception("The matrix must be square");
            if (n >= height || n < 0) throw new Exception("Row and column number specified is invalid");

            // new matrix:
            double[,] H = new double[height - 1, width - 1];
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

        #region Diff voids
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Order</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        public static double[,] Diff(double[,] a, int n, Direction direction, bool reverse = false)
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
        private static double[,] DiffVertical(double[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new double[0, m];

            // new array
            double[,] y = new double[r, m];
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
        private static double[,] DiffHorizontal(double[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new double[m, 0];

            // new array
            double[,] y = new double[m, c];
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
        public static Complex[,] Diff(Complex[,] a, int n, Direction direction, bool reverse = false)
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
        private static Complex[,] DiffVertical(Complex[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new Complex[0, m];

            // new array
            Complex[,] y = new Complex[r, m];
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
        private static Complex[,] DiffHorizontal(Complex[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new Complex[m, 0];

            // new array
            Complex[,] y = new Complex[m, c];
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
        public static double[] Diff(this double[] v, int n, bool reverse = false)
        {
            // start
            double[] z;
            double[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new double[0];

                y = new double[length];

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
        public static Complex[] Diff(this Complex[] v, int n, bool reverse = false)
        {
            // start
            Complex[] z;
            Complex[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new Complex[0];

                y = new Complex[length];

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

        #region Extend voids
        /// <summary>
        /// Extends the vector to the specified length.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Extend(this double[] v, int length)
        {
            int r0 = v.GetLength(0);
            int rr = (length - r0) / 2;
            int dr = length - rr;
            double[] b = new double[length];
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
        public static double[,] Extend(this double[,] m, int height, int width)
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
        private static double[,] ExtendVertical(double[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int rr = (length - r0) / 2;
            double[,] B = new double[length, c0];
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
        private static double[,] ExtendHorizontal(double[,] m, int length)
        {
            // params
            int r0 = m.GetLength(0);
            int c0 = m.GetLength(1);
            int cc = (length - c0) / 2;
            double[,] B = new double[r0, length];
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
        public static double[] Compute(double min, double max, double step)
        {
            // ******************************
            // MATLAB vector computing void
            // designed by Valery Asiryan
            // ******************************

            // shifts and variables:
            double dy = max - min + step, i;
            int dx = (int)Maths.Fix(dy / step), j;

            // C# has a significant bug, which you can check with:
            // min = 0.5, max = 1, step = 0.001,
            // maxz = max, j = 63.

            // output vector and eps:
            double[] x = new double[dx];
            double eps = max / 1e8, maxz;

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
        public static double[] Compute(this double[] v, IDouble function)
        {
            int length = v.Length;
            double[] H = new double[length];

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
        public static double[,] Compute(this double[] x, double[] y, IDoubleMesh function)
        {
            int xlength = x.Length, ylength = y.Length;
            double[,] z = new double[xlength, ylength];
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
        public static Complex[,] Compute(this double[] x, Complex[] y, IComplexMesh function)
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
        public static Complex[,] Compute(this Complex[] x, double[] y, IComplexMesh function)
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
        public static double[,] Compute(this double[,] m, IDouble function)
        {
            int i, j;
            int ml = m.GetLength(1), mr = m.GetLength(0);
            double[,] H = new double[mr, ml];

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
        public static double[] One(int n)
        {
            double[] v = new double[n];

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
        public static double[] Zero(int n)
        {
            return new double[n];
        }
        #endregion

        #region Matrix products
        /// <summary>
        /// Returns the Householder vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Array</returns>
        public static double[] Householder(this double[] v)
        {
            int length = v.Length;

            // length checking:
            if (length > 0)
            {
                double norm = v.Norm();
                double[] u = new double[length];

                // householder vector:
                if (norm != 0)
                {
                    u[0] = v[0] / norm;
                    u[0] = u[0] + Maths.Sign(u[0]);
                    u[0] = u[0] / Math.Sqrt(Math.Abs(u[0]));

                    for (int i = 1; i < length; i++)
                    {
                        u[i] = v[i] / norm;
                        u[i] = u[i] / u[0];
                    }
                }
                else
                {
                    u = v;
                    u[0] = Maths.Sqrt2;
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
        public static double[,] Companion(this double[] v)
        {
            int n = v.Length, i;
            double[,] H = new double[n, n];

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
        public static double[,] Vander(this double[] v)
        {
            int n = v.Length;
            double[,] H = new double[n, n];
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    H[i, j] = Math.Pow(v[i], j);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of an incomplete Hankel matrix.
        /// </summary>
        /// <param name="v">Array</param>
        /// <returns>Matrix</returns>
        public static double[,] Hankeli(this double[] v)
        {
            int n = v.Length;
            double[,] H = new double[n, n];
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
        public static double[,] Hankel(this double[] v)
        {
            int n = v.Length / 2;
            double[,] H = new double[n, n];
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
        public static double[,] Toeplitz(this double[] v)
        {
            int n = v.Length / 2;
            double[,] H = new double[n, n];
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
        public static double[,] Cauchy(this double[] x, double[] y)
        {
            int m = x.Length, l = y.Length;
            double[,] H = new double[m, l];
            double kern;
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
        public static double[,] Circulant(this double[] v)
        {
            int n = v.Length;
            double[,] H = new double[n, n];
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
        public static double[,] Symmetric(this double[] v)
        {
            int n = v.Length;
            double[,] H = new double[n, n];
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
        public static double[,] Zero(int m, int l)
        {
            return new double[m, l];
        }
        /// <summary>
        /// Implements the construction of a eye matrix.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[,] Eye(int m, int l)
        {
            double[,] H = new double[m, l];
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
        public static double[,] One(int m, int l)
        {
            double[,] H = new double[m, l];
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
        public static double[,] Exchange(int n)
        {
            double[,] H = new double[n, n];
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
        public static double[,] Lehmer(int n)
        {
            double[,] H = new double[n, n];
            int i, j;
            double x, y;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    x = (double)i; y = (double)j;
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
        public static double[,] Redheffer(int n)
        {
            double[,] H = new double[n, n];
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
        public static double[,] Hilbert(int n)
        {
            double[,] H = new double[n, n];
            int i, j;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    H[i - 1, j - 1] = 1.0 / (i + j - 1);
                }
            }
            return H;
        }
        /// <summary>
        /// Implements the construction of a cyclic matrix.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static double[,] Circulant(int n)
        {
            double[,] H = new double[n, n];
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
        public static double[,] Symmetric(int n)
        {
            double[,] H = new double[n, n];
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
        public static double[,] GCD(int n)
        {
            double[,] H = new double[n, n];
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
        public static double[,] Stirling(int n, bool second = false)
        {
            // Stirling's matrix 
            // of the second kind
            double[,] S = new double[n, n];
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
        public static double[,] Magic(int n)
        {
            if (Maths.Mod(n, 2) != 1)
                throw new Exception("Dimension of the matrix must be an odd number");

            double[,] m = new double[n, n];
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
        /// Implements the construction of a vector of random numbers, the values of which are distributed according to a uniform distribution.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static double[] Rand(int n)
        {
            double[] v = new double[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = rnd.NextDouble();
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a vector of random numbers, the values of which are distributed according to a uniform distribution.
        /// </summary>
        /// <param name="n">Dimension</param>
        /// <returns>Array</returns>
        public static Complex[] Randc(int n)
        {
            Complex[] v = new Complex[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = new Complex(rnd.NextDouble(), rnd.NextDouble());
            }

            return v;
        }
        /// <summary>
        /// Implements the construction of a matrix of random numbers, the values of which are distributed according to a uniform distribution.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[,] Rand(int m, int l)
        {
            double[,] H = new double[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = rnd.NextDouble();
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of random numbers, the values of which are distributed according to a uniform distribution.
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
                    H[i, j] = new Complex(rnd.NextDouble(), rnd.NextDouble());
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
        public static double[] Randi(int n)
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
        public static double[] Randi(int n, int a, int b)
        {
            double[] v = new double[n];

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
        public static double[,] Randi(int m, int l)
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
        public static double[,] Randi(int m, int l, int a, int b)
        {
            double[,] H = new double[m, l];
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
        /// Parses the original string into a matrix of double numbers.
        /// <remarks>
        /// Example: "[1, 2, 3; 4, 5, 6; 7, 8, 9]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static double[,] Parse(this double[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length, n = nums.Length, k;
            double[,] H = new double[r, n];
            int i, j;

            // first row
            for (j = 0; j < n; j++)
            {
                H[0, j] = double.Parse(nums[j]);
            }

            // other rows
            for (i = 1; i < r; i++)
            {
                nums = rows[i].Split('|');
                k = Math.Min(n, nums.Length);

                for (j = 0; j < k; j++)
                {
                    H[i, j] = double.Parse(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of double numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, ref double[,] result)
        {
            double[,] zero = null;
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
        public static bool TryParse(string s, ref Complex[,] result)
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
        /// Parses the original string into a vector of double numbers.
        /// <remarks>
        /// Example: "[1, 2, 3, 4]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static double[] Parse(this double[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
                string[] nums = rows[0].Split('|');
                int n = nums.Length, i;
                double[] H = new double[n];

                // collecting rows:
                for (i = 0; i < n; i++)
                {
                    H[i] = double.Parse(nums[i]);
                }
                return H;
            }
            else
            {
                throw new Exception("The input string was in the wrong format");
            }
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of double numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, ref double[] result)
        {
            double[] zero = null;
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
        public static bool TryParse(string s, ref Complex[] result)
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
        public static double[] Solve(this double[,] A)
        {
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height + 1 != width)
                throw new Exception("Input matrix has invalid sizes");

            double[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            double[] x = new double[height];
            double[] v, w;
            double temp;

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
        public static double[] Solve(this double[,] A, double[] b)
        {
            int height = A.GetLength(0);
            int width = A.GetLength(1);

            if (height != width)
                throw new Exception("The matrix must be square");
            if (height != b.Length)
                throw new Exception("Vector length should be equal to the height of the matrix");

            double[][] B = Jagged.ToJagged(A);
            int i, j, k, l;
            double[] x = (double[])b.Clone();
            double[] v, w;
            double temp;


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
