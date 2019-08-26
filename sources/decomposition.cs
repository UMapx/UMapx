// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                            UMAPX.DECOMPOSITION
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Matrix decompositions
    /// <summary>
    /// Определяет разложение с приведением к форме Хессенберга.
    /// <remarks>
    /// Это представление квадратной матрицы в виде произведения трех матриц: A = P * H * P', где H - форма Хессенберга, а P - унитарная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Hessenberg_matrix
    /// </remarks>
    /// </summary>
    public class Hessenberg
    {
        #region Private data
        private double[][] matrices;
        private double[][] hessenberg;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует разложение с приведением к форме Хессенберга.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public Hessenberg(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // Reduce to Hessenberg form.
            orthes(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает унитарную матрицу.
        /// </summary>
        public double[,] P
        {
            get { return Jagged.FromJagged(matrices); }
        }
        /// <summary>
        /// Получает форму Хессенберга.
        /// </summary>
        public double[,] H
        {
            get { return Jagged.FromJagged(hessenberg); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Несимметричное сокращение до формы Хессенберга.
        /// </summary>
        /// <param name="A">Матрица</param>
        private void orthes(double[,] A)
        {
            // Properties
            int n = A.GetLength(0);
            this.matrices = Jagged.Zero(n, n);
            this.hessenberg = Jagged.ToJagged(A);
            double[] orthogonal = new double[n];

            // Nonsymmetric reduction to Hessenberg form.
            // This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
            int low = 0;
            int high = n - 1;
            int m, i, j;
            double scale, h, g, f;

            for (m = low + 1; m <= high - 1; m++)
            {
                // Scale column.

                scale = 0;
                for (i = m; i <= high; i++)
                    scale = scale + System.Math.Abs(hessenberg[i][m - 1]);

                if (scale != 0)
                {
                    // Compute Householder transformation.
                    h = 0;
                    for (i = high; i >= m; i--)
                    {
                        orthogonal[i] = hessenberg[i][m - 1] / scale;
                        h += orthogonal[i] * orthogonal[i];
                    }

                    g = (Double)System.Math.Sqrt(h);
                    if (orthogonal[m] > 0) g = -g;

                    h = h - orthogonal[m] * g;
                    orthogonal[m] = orthogonal[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u') / h)
                    for (j = m; j < n; j++)
                    {
                        f = 0;
                        for (i = high; i >= m; i--)
                            f += orthogonal[i] * hessenberg[i][j];

                        f = f / h;
                        for (i = m; i <= high; i++)
                            hessenberg[i][j] -= f * orthogonal[i];
                    }

                    for (i = 0; i <= high; i++)
                    {
                        f = 0;
                        for (j = high; j >= m; j--)
                            f += orthogonal[j] * hessenberg[i][j];

                        f = f / h;
                        for (j = m; j <= high; j++)
                            hessenberg[i][j] -= f * orthogonal[j];
                    }

                    orthogonal[m] = scale * orthogonal[m];
                    hessenberg[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    matrices[i][j] = (i == j ? 1 : 0);

            for (m = high - 1; m >= low + 1; m--)
            {
                if (hessenberg[m][m - 1] != 0)
                {
                    for (i = m + 1; i <= high; i++)
                        orthogonal[i] = hessenberg[i][m - 1];

                    for (j = m; j <= high; j++)
                    {
                        g = 0;
                        for (i = m; i <= high; i++)
                            g += orthogonal[i] * matrices[i][j];

                        // Double division avoids possible underflow.
                        g = (g / orthogonal[m]) / hessenberg[m][m - 1];
                        for (i = m; i <= high; i++)
                            matrices[i][j] += g * orthogonal[i];
                    }
                }
            }

            // final reduction:
            if (n > 2)
            {
                for (i = 0; i < n - 2; i++)
                {
                    for (j = i + 2; j < n; j++)
                    {
                        hessenberg[j][i] = 0;
                    }
                }
            }

        }
        #endregion
    }
    /// <summary>
    /// Определяет обобщенное спектральное разложение квадратных матриц A и B.
    /// <remarks>
    /// Это задача нахождение такого вектора значений V, что справедливо представление: A * V = B * V * D. 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem
    /// </remarks>
    /// </summary>
    public class GEVD
    {
        #region Private data
        private int n;
        private double[] ar;
        private double[] ai;
        private double[] beta;
        private double[][] Z;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует обобщенное спектральное разложение квадратных матриц A и B.
        /// </summary>
        /// <param name="a">Матрица A</param>
        /// <param name="b">Матрица B</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public GEVD(double[,] a, double[,] b, double eps = 1e-16)
        {
            if (a.GetLength(0) != a.GetLength(1))
                throw new ArgumentException("Матрица должна быть квадратной");

            if (b.GetLength(0) != b.GetLength(1))
                throw new ArgumentException("Матрица должна быть квадратной");

            if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1))
                throw new ArgumentException("Матрицы должны быть одинаковых размеров");

            // params
            this.n = a.GetLength(0);
            this.ar = new double[n];
            this.ai = new double[n];
            this.beta = new double[n];
            this.Z = Jagged.Zero(n, n);
            var A = Jagged.ToJagged(a);
            var B = Jagged.ToJagged(b);
            bool matz = true;
            int ierr = 0;


            // reduces A to upper Hessenberg form and B to upper
            // triangular form using orthogonal transformations
            qzhes(n, A, B, matz, Z);

            // reduces the Hessenberg matrix A to quasi-triangular form
            // using orthogonal transformations while maintaining the
            // triangular form of the B matrix.
            qzit(n, A, B, Maths.Double(eps), matz, Z, ref ierr);

            // reduces the quasi-triangular matrix further, so that any
            // remaining 2-by-2 blocks correspond to pairs of complex
            // eigenvalues, and returns quantities whose ratios give the
            // generalized eigenvalues.
            qzval(n, A, B, ar, ai, beta, matz, Z);

            // computes the eigenvectors of the triangular problem and
            // transforms the results back to the original coordinate system.
            qzvec(n, A, B, ar, ai, beta, Z);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает вектор α.
        /// </summary>
        public Complex[] Alpha
        {
            get
            {
                Complex[] a = new Complex[n];

                for (int i = 0; i < n; i++)
                    a[i] = new Complex(ar[i], ai[i]);

                return a;
            }
        }
        /// <summary>
        /// Возвращает вектор β.
        /// </summary>
        public double[] Beta
        {
            get { return beta; }
        }
        /// <summary>
        /// Возвращает вектор собственных значений.
        /// </summary>
        public Complex[] Eigenvalues
        {
            get
            {
                // c = (ar + j*ai) / beta
                Complex[] c = new Complex[n];
                for (int i = 0; i < n; i++)
                    c[i] = new Complex(ar[i], ai[i]) / beta[i];

                return c;
            }
        }
        /// <summary>
        /// Возвращает матрицу собственных значений D.
        /// </summary>
        public double[,] D
        {
            get
            {
                double[,] x = new double[n, n];

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                        x[i, j] = 0.0;

                    x[i, i] = ar[i] / beta[i];
                    if (ai[i] > 0)
                        x[i, i + 1] = ai[i] / beta[i];
                    else if (ai[i] < 0)
                        x[i, i - 1] = ai[i] / beta[i];
                }

                return x;
            }
        }
        /// <summary>
        /// Возвращает матрицу значений V.
        /// </summary>
        public double[,] V
        {
            get
            {
                return Z.FromJagged();
            }
        }
        /// <summary>
        /// Проверяет является ли одна из матриц сингулярной или нет.
        /// </summary>
        public bool IsSingular
        {
            get
            {
                for (int i = 0; i < n; i++)
                    if (beta[i] == 0)
                        return true;
                return false;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        ///   Adaptation of the original Fortran QZHES routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the first step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///
        ///   This subroutine accepts a pair of real general matrices and
        ///   reduces one of them to upper Hessenberg form and the other
        ///   to upper triangular form using orthogonal transformations.
        ///   it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="matz"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        private static void qzhes(int n, double[][] a, double[][] b, bool matz, double[][] z)
        {
            int i, j, k, l;
            double r, s, t;
            int l1;
            double u1, u2, v1, v2;
            int lb, nk1;
            double rho;

            // Reduce b to upper triangular form
            if (n <= 1) return;
            for (l = 0; l < n - 1; ++l)
            {
                l1 = l + 1;
                s = 0.0;

                for (i = l1; i < n; ++i)
                    s += (System.Math.Abs(b[i][l]));

                if (s == 0.0) continue;
                s += (System.Math.Abs(b[l][l]));
                r = 0.0;

                for (i = l; i < n; ++i)
                {
                    // Computing 2nd power
                    b[i][l] /= s;
                    r += b[i][l] * b[i][l];
                }

                r = Sign(System.Math.Sqrt(r), b[l][l]);
                b[l][l] += r;
                rho = r * b[l][l];

                for (j = l1; j < n; ++j)
                {
                    t = 0.0;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * b[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        b[i][j] += t * b[i][l];
                }

                for (j = 0; j < n; ++j)
                {
                    t = 0.0;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * a[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        a[i][j] += t * b[i][l];
                }

                b[l][l] = -s * r;
                for (i = l1; i < n; ++i)
                    b[i][l] = 0.0;
            }

            // Reduce a to upper Hessenberg form, while keeping b triangular
            if (n == 2) return;
            for (k = 0; k < n - 2; ++k)
            {
                nk1 = n - 2 - k;

                // for l=n-1 step -1 until k+1 do
                for (lb = 0; lb < nk1; ++lb)
                {
                    l = n - lb - 2;
                    l1 = l + 1;

                    // Zero a(l+1,k)
                    s = (System.Math.Abs(a[l][k])) + (System.Math.Abs(a[l1][k]));

                    if (s == 0.0) continue;
                    u1 = a[l][k] / s;
                    u2 = a[l1][k] / s;
                    r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (j = k; j < n; ++j)
                    {
                        t = a[l][j] + u2 * a[l1][j];
                        a[l][j] += t * v1;
                        a[l1][j] += t * v2;
                    }

                    a[l1][k] = 0.0;

                    for (j = l; j < n; ++j)
                    {
                        t = b[l][j] + u2 * b[l1][j];
                        b[l][j] += t * v1;
                        b[l1][j] += t * v2;
                    }

                    // Zero b(l+1,l)
                    s = (System.Math.Abs(b[l1][l1])) + (System.Math.Abs(b[l1][l]));

                    if (s == 0.0) continue;
                    u1 = b[l1][l1] / s;
                    u2 = b[l1][l] / s;
                    r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (i = 0; i <= l1; ++i)
                    {
                        t = b[i][l1] + u2 * b[i][l];
                        b[i][l1] += t * v1;
                        b[i][l] += t * v2;
                    }

                    b[l1][l] = 0.0;

                    for (i = 0; i < n; ++i)
                    {
                        t = a[i][l1] + u2 * a[i][l];
                        a[i][l1] += t * v1;
                        a[i][l] += t * v2;
                    }

                    if (matz)
                    {
                        for (i = 0; i < n; ++i)
                        {
                            t = z[i][l1] + u2 * z[i][l];
                            z[i][l1] += t * v1;
                            z[i][l] += t * v2;
                        }
                    }
                }
            }

            return;
        }
        /// <summary>
        ///   Adaptation of the original Fortran QZIT routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the second step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart,
        ///   as modified in technical note nasa tn d-7305(1973) by ward.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them
        ///   in upper Hessenberg form and the other in upper triangular form.
        ///   it reduces the Hessenberg matrix to quasi-triangular form using
        ///   orthogonal transformations while maintaining the triangular form
        ///   of the other matrix.  it is usually preceded by  qzhes  and
        ///   followed by  qzval  and, possibly,  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="eps1"></param>
        /// <param name="matz"></param>
        /// <param name="z"></param>
        /// <param name="ierr"></param>
        /// <returns></returns>
        private static void qzit(int n, double[][] a, double[][] b, double eps1, bool matz, double[][] z, ref int ierr)
        {

            int i, j, k, l = 0;
            double r, s, t, a1, a2, a3 = 0;
            int k1, k2, l1, ll;
            double u1, u2, u3;
            double v1, v2, v3;
            double a11, a12, a21, a22, a33, a34, a43, a44;
            double b11, b12, b22, b33, b34, b44;
            int na, en, ld;
            double ep;
            double sh = 0;
            int km1, lm1 = 0;
            double ani, bni;
            int ish, itn, its, enm2, lor1;
            double epsa, epsb, anorm = 0, bnorm = 0;
            int enorn;
            bool notlas;


            ierr = 0;

            #region Compute epsa and epsb
            for (i = 0; i < n; ++i)
            {
                ani = 0.0;
                bni = 0.0;

                if (i != 0)
                    ani = (Math.Abs(a[i][(i - 1)]));

                for (j = i; j < n; ++j)
                {
                    ani += Math.Abs(a[i][j]);
                    bni += Math.Abs(b[i][j]);
                }

                if (ani > anorm) anorm = ani;
                if (bni > bnorm) bnorm = bni;
            }

            if (anorm == 0.0) anorm = 1.0;
            if (bnorm == 0.0) bnorm = 1.0;

            ep = eps1;
            if (ep == 0.0)
            {
                // Use round-off level if eps1 is zero
                ep = Epsilon(1.0);
            }

            epsa = ep * anorm;
            epsb = ep * bnorm;
            #endregion


            // Reduce a to quasi-triangular form, while keeping b triangular
            lor1 = 0;
            enorn = n;
            en = n - 1;
            itn = n * 30;

        // Begin QZ step
        L60:
            if (en <= 1) goto L1001;
            if (!matz) enorn = en + 1;

            its = 0;
            na = en - 1;
            enm2 = na;

        L70:
            ish = 2;
            // Check for convergence or reducibility.
            for (ll = 0; ll <= en; ++ll)
            {
                lm1 = en - ll - 1;
                l = lm1 + 1;

                if (l + 1 == 1)
                    goto L95;

                if ((Math.Abs(a[l][lm1])) <= epsa)
                    break;
            }

        L90:
            a[l][lm1] = 0.0;
            if (l < na) goto L95;

            // 1-by-1 or 2-by-2 block isolated
            en = lm1;
            goto L60;

        // Check for small top of b 
        L95:
            ld = l;

        L100:
            l1 = l + 1;
            b11 = b[l][l];

            if (Math.Abs(b11) > epsb) goto L120;

            b[l][l] = 0.0;
            s = (Math.Abs(a[l][l]) + Math.Abs(a[l1][l]));
            u1 = a[l][l] / s;
            u2 = a[l1][l] / s;
            r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
            v1 = -(u1 + r) / r;
            v2 = -u2 / r;
            u2 = v2 / v1;

            for (j = l; j < enorn; ++j)
            {
                t = a[l][j] + u2 * a[l1][j];
                a[l][j] += t * v1;
                a[l1][j] += t * v2;

                t = b[l][j] + u2 * b[l1][j];
                b[l][j] += t * v1;
                b[l1][j] += t * v2;
            }

            if (l != 0)
                a[l][lm1] = -a[l][lm1];

            lm1 = l;
            l = l1;
            goto L90;

        L120:
            a11 = a[l][l] / b11;
            a21 = a[l1][ l] / b11;
            if (ish == 1) goto L140;

            // Iteration strategy
            if (itn == 0) goto L1000;
            if (its == 10) goto L155;

            // Determine type of shift
            b22 = b[l1][l1];
            if (Math.Abs(b22) < epsb) b22 = epsb;
            b33 = b[na][na];
            if (Math.Abs(b33) < epsb) b33 = epsb;
            b44 = b[en][en];
            if (Math.Abs(b44) < epsb) b44 = epsb;
            a33 = a[na][na] / b33;
            a34 = a[na][en] / b44;
            a43 = a[en][na] / b33;
            a44 = a[en][en] / b44;
            b34 = b[na][en] / b44;
            t = (a43 * b34 - a33 - a44) * .5;
            r = t * t + a34 * a43 - a33 * a44;
            if (r < 0.0) goto L150;

            // Determine single shift zero-th column of a
            ish = 1;
            r = Math.Sqrt(r);
            sh = -t + r;
            s = -t - r;
            if (Math.Abs(s - a44) < Math.Abs(sh - a44))
                sh = s;

            // Look for two consecutive small sub-diagonal elements of a.
            for (ll = ld; ll + 1 <= enm2; ++ll)
            {
                l = enm2 + ld - ll - 1;

                if (l == ld)
                    goto L140;

                lm1 = l - 1;
                l1 = l + 1;
                t = a[l + 1][l + 1];

                if (Math.Abs(b[l][l]) > epsb)
                    t -= sh * b[l][l];

                if (Math.Abs(a[l][lm1]) <= (Math.Abs(t / a[l1][l])) * epsa)
                    goto L100;
            }

        L140:
            a1 = a11 - sh;
            a2 = a21;
            if (l != ld)
                a[l][lm1] = -a[l][lm1];
            goto L160;

        // Determine double shift zero-th column of a
        L150:
            a12 = a[l][l1] / b22;
            a22 = a[l1][l1] / b22;
            b12 = b[l][l1] / b22;
            a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + a12 - a11 * b12;
            a2 = a22 - a11 - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
            a3 = a[l1 + 1][l1] / b22;
            goto L160;

        // Ad hoc shift
        L155:
            a1 = 0.0;
            a2 = 1.0;
            a3 = 1.1605;

        L160:
            ++its;
            --itn;

            if (!matz) lor1 = ld;

            // Main loop
            for (k = l; k <= na; ++k)
            {
                notlas = k != na && ish == 2;
                k1 = k + 1;
                k2 = k + 2;

                km1 = Math.Max(k, l + 1) - 1; // Computing MAX
                ll = Math.Min(en, k1 + ish);  // Computing MIN

                if (notlas) goto L190;

                // Zero a(k+1,k-1)
                if (k == l) goto L170;
                a1 = a[k][km1];
                a2 = a[k1][km1];

            L170:
                s = Math.Abs(a1) + Math.Abs(a2);
                if (s == 0.0) goto L70;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;

                    t = b[k][j] + u2 * b[k1][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                }

                if (k != l)
                    a[k1][km1] = 0.0;
                goto L240;

                // Zero a(k+1,k-1) and a(k+2,k-1)
            L190:
                if (k == l) goto L200;
                a1 = a[k][ km1];
                a2 = a[k1][km1];
                a3 = a[k2][km1];

            L200:
                s = Math.Abs(a1) + Math.Abs(a2) + Math.Abs(a3);
                if (s == 0.0) goto L260;
                u1 = a1 / s;
                u2 = a2 / s;
                u3 = a3 / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j] + u3 * a[k2][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;
                    a[k2][j] += t * v3;

                    t = b[k][j] + u2 * b[k1][j] + u3 * b[k2][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                    b[k2][j] += t * v3;
                }

                if (k == l) goto L220;
                a[k1][km1] = 0.0;
                a[k2][km1] = 0.0;

            // Zero b(k+2,k+1) and b(k+2,k)
            L220:
                s = (Math.Abs(b[k2][k2])) + (Math.Abs(b[k2][k1])) + (Math.Abs(b[k2][k]));
                if (s == 0.0) goto L240;
                u1 = b[k2][k2] / s;
                u2 = b[k2][k1] / s;
                u3 = b[k2][k] / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k2] + u2 * a[i][k1] + u3 * a[i][k];
                    a[i][k2] += t * v1;
                    a[i][k1] += t * v2;
                    a[i][k] += t * v3;

                    t = b[i][k2] + u2 * b[i][k1] + u3 * b[i][k];
                    b[i][k2] += t * v1;
                    b[i][k1] += t * v2;
                    b[i][k] += t * v3;
                }

                b[k2][k] = 0.0;
                b[k2][k1] = 0.0;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k2] + u2 * z[i][k1] + u3 * z[i][k];
                        z[i][k2] += t * v1;
                        z[i][k1] += t * v2;
                        z[i][k] += t * v3;
                    }
                }

            // Zero b(k+1,k)
            L240:
                s = (Math.Abs(b[k1][k1])) + (Math.Abs(b[k1][k]));
                if (s == 0.0) goto L260;
                u1 = b[k1][k1] / s;
                u2 = b[k1][k] / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k1] + u2 * a[i][k];
                    a[i][k1] += t * v1;
                    a[i][k] += t * v2;

                    t = b[i][k1] + u2 * b[i][k];
                    b[i][k1] += t * v1;
                    b[i][k] += t * v2;
                }

                b[k1][k] = 0.0;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k1] + u2 * z[i][k];
                        z[i][k1] += t * v1;
                        z[i][k] += t * v2;
                    }
                }

            L260:
                ;
            }

            goto L70; // End QZ step

        // Set error -- all eigenvalues have not converged after 30*n iterations
        L1000:
            ierr = en + 1;

        // Save epsb for use by qzval and qzvec
        L1001:
            if (n > 1)
                b[n - 1][0] = epsb;
            return;
        }
        /// <summary>
        ///   Adaptation of the original Fortran QZVAL routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the third step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them
        ///   in quasi-triangular form and the other in upper triangular form.
        ///   it reduces the quasi-triangular matrix further, so that any
        ///   remaining 2-by-2 blocks correspond to pairs of complex
        ///   Eigenvalues, and returns quantities whose ratios give the
        ///   generalized eigenvalues.  it is usually preceded by  qzhes
        ///   and  qzit  and may be followed by  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="alfr"></param>
        /// <param name="alfi"></param>
        /// <param name="beta"></param>
        /// <param name="matz"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        private static void qzval(int n, double[][] a, double[][] b, double[] alfr, double[] alfi, double[] beta, bool matz, double[][] z)
        {
            int i, j;
            int na, en, nn;
            double c, d, e = 0;
            double r, s, t;
            double a1, a2, u1, u2, v1, v2;
            double a11, a12, a21, a22;
            double b11, b12, b22;
            double di, ei;
            double an = 0, bn;
            double cq, dr;
            double cz, ti, tr;
            double a1i, a2i, a11i, a12i, a22i, a11r, a12r, a22r;
            double sqi, ssi, sqr, szi, ssr, szr;

            double epsb = b[n - 1][0];
            int isw = 1;


            // Find eigenvalues of quasi-triangular matrices.
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;

                if (isw == 2) goto L505;
                if (en == 0) goto L410;
                if (a[en][na] != 0.0) goto L420;

            // 1-by-1 block, one real root
            L410:
                alfr[en] = a[en][en];
                if (b[en][ en] < 0.0)
                {
                    alfr[en] = -alfr[en];
                }
                beta[en] = (Math.Abs(b[en][en]));
                alfi[en] = 0.0;
                goto L510;

            // 2-by-2 block
            L420:
                if (Math.Abs(b[na][na]) <= epsb) goto L455;
                if (Math.Abs(b[en][ en]) > epsb) goto L430;
                a1 = a[en][en];
                a2 = a[en][na];
                bn = 0.0;
                goto L435;

            L430:
                an = Math.Abs(a[na][na]) + Math.Abs(a[na][en]) + Math.Abs(a[en][na]) + Math.Abs(a[en][en]);
                bn = Math.Abs(b[na][na]) + Math.Abs(b[na][en]) + Math.Abs(b[en][en]);
                a11 = a[na][na] / an;
                a12 = a[na][en] / an;
                a21 = a[en][na] / an;
                a22 = a[en][en] / an;
                b11 = b[na][na] / bn;
                b12 = b[na][en] / bn;
                b22 = b[en][en] / bn;
                e = a11 / b11;
                ei = a22 / b22;
                s = a21 / (b11 * b22);
                t = (a22 - e * b22) / b22;

                if (Math.Abs(e) <= Math.Abs(ei))
                    goto L431;

                e = ei;
                t = (a11 - e * b11) / b11;

            L431:
                c = (t - s * b12) * .5;
                d = c * c + s * (a12 - e * b12);
                if (d < 0.0) goto L480;

                // Two real roots. Zero both a(en,na) and b(en,na)
                e += c + Sign(Math.Sqrt(d), c);
                a11 -= e * b11;
                a12 -= e * b12;
                a22 -= e * b22;

                if (Math.Abs(a11) + Math.Abs(a12) < Math.Abs(a21) + Math.Abs(a22))
                    goto L432;

                a1 = a12;
                a2 = a11;
                goto L435;

            L432:
                a1 = a22;
                a2 = a21;

            // Choose and apply real z
            L435:
                s = Math.Abs(a1) + Math.Abs(a2);
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = 0; i <= en; ++i)
                {
                    t = a[i][en] + u2 * a[i][na];
                    a[i][en] += t * v1;
                    a[i][na] += t * v2;

                    t = b[i][en] + u2 * b[i][na];
                    b[i][en] += t * v1;
                    b[i][na] += t * v2;
                }

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][en] + u2 * z[i][na];
                        z[i][en] += t * v1;
                        z[i][na] += t * v2;
                    }
                }

                if (bn == 0.0) goto L475;
                if (an < System.Math.Abs(e) * bn) goto L455;
                a1 = b[na][na];
                a2 = b[en][na];
                goto L460;

            L455:
                a1 = a[na][na];
                a2 = a[en][na];

            // Choose and apply real q
            L460:
                s = System.Math.Abs(a1) + System.Math.Abs(a2);
                if (s == 0.0) goto L475;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = na; j < n; ++j)
                {
                    t = a[na][j] + u2 * a[en][j];
                    a[na][j] += t * v1;
                    a[en][j] += t * v2;

                    t = b[na][j] + u2 * b[en][j];
                    b[na][j] += t * v1;
                    b[en][j] += t * v2;
                }

            L475:
                a[en][na] = 0.0;
                b[en][na] = 0.0;
                alfr[na] = a[na][na];
                alfr[en] = a[en][en];

                if (b[na][na] < 0.0)
                    alfr[na] = -alfr[na];

                if (b[en][en] < 0.0)
                    alfr[en] = -alfr[en];

                beta[na] = (System.Math.Abs(b[na][na]));
                beta[en] = (System.Math.Abs(b[en][en]));
                alfi[en] = 0.0;
                alfi[na] = 0.0;
                goto L505;

                // Two complex roots
            L480:
                e += c;
                ei = System.Math.Sqrt(-d);
                a11r = a11 - e * b11;
                a11i = ei * b11;
                a12r = a12 - e * b12;
                a12i = ei * b12;
                a22r = a22 - e * b22;
                a22i = ei * b22;

                if (System.Math.Abs(a11r) + System.Math.Abs(a11i) +
                    System.Math.Abs(a12r) + System.Math.Abs(a12i) <
                    System.Math.Abs(a21) + System.Math.Abs(a22r)
                    + System.Math.Abs(a22i))
                    goto L482;

                a1 = a12r;
                a1i = a12i;
                a2 = -a11r;
                a2i = -a11i;
                goto L485;

            L482:
                a1 = a22r;
                a1i = a22i;
                a2 = -a21;
                a2i = 0.0;

            // Choose complex z
            L485:
                cz = System.Math.Sqrt(a1 * a1 + a1i * a1i);
                if (cz == 0.0) goto L487;
                szr = (a1 * a2 + a1i * a2i) / cz;
                szi = (a1 * a2i - a1i * a2) / cz;
                r = System.Math.Sqrt(cz * cz + szr * szr + szi * szi);
                cz /= r;
                szr /= r;
                szi /= r;
                goto L490;

            L487:
                szr = 1.0;
                szi = 0.0;

            L490:
                if (an < (System.Math.Abs(e) + ei) * bn) goto L492;
                a1 = cz * b11 + szr * b12;
                a1i = szi * b12;
                a2 = szr * b22;
                a2i = szi * b22;
                goto L495;

            L492:
                a1 = cz * a11 + szr * a12;
                a1i = szi * a12;
                a2 = cz * a21 + szr * a22;
                a2i = szi * a22;

            // Choose complex q
            L495:
                cq = System.Math.Sqrt(a1 * a1 + a1i * a1i);
                if (cq == 0.0) goto L497;
                sqr = (a1 * a2 + a1i * a2i) / cq;
                sqi = (a1 * a2i - a1i * a2) / cq;
                r = System.Math.Sqrt(cq * cq + sqr * sqr + sqi * sqi);
                cq /= r;
                sqr /= r;
                sqi /= r;
                goto L500;

            L497:
                sqr = 1.0;
                sqi = 0.0;

            // Compute diagonal elements that would result if transformations were applied
            L500:
                ssr = sqr * szr + sqi * szi;
                ssi = sqr * szi - sqi * szr;
                i = 0;
                tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
                ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
                dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
                di = cq * szi * b12 + ssi * b22;
                goto L503;

            L502:
                i = 1;
                tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
                ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
                dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
                di = -ssi * b11 - sqi * cz * b12;

            L503:
                t = ti * dr - tr * di;
                j = na;

                if (t < 0.0)
                    j = en;

                r = Math.Sqrt(dr * dr + di * di);
                beta[j] = bn * r;
                alfr[j] = an * (tr * dr + ti * di) / r;
                alfi[j] = an * t / r;
                if (i == 0) goto L502;

            L505:
                isw = 3 - isw;

            L510:
                ;
            }

            b[n - 1][0] = epsb;

            return;
        }
        /// <summary>
        ///   Adaptation of the original Fortran QZVEC routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the optional fourth step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them in
        ///   quasi-triangular form (in which each 2-by-2 block corresponds to
        ///   a pair of complex eigenvalues) and the other in upper triangular
        ///   form.  It computes the eigenvectors of the triangular problem and
        ///   transforms the results back to the original coordinate system.
        ///   it is usually preceded by  qzhes,  qzit, and  qzval.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="alfr"></param>
        /// <param name="alfi"></param>
        /// <param name="beta"></param>
        /// <param name="z"></param>
        /// <returns></returns>
        private static void qzvec(int n, double[][] a, double[][] b, double[] alfr, double[] alfi, double[] beta, double[][] z)
        {
            int i, j, k, m;
            int na, ii, en, jj, nn, enm2;
            double d, q;
            double r = 0, s = 0, t, w, x = 0, y, t1, t2, w1, x1 = 0, z1 = 0, di;
            double ra, dr, sa;
            double ti, rr, tr, zz = 0;
            double alfm, almi, betm, almr;

            double epsb = b[n - 1][0];
            int isw = 1;


            // for en=n step -1 until 1 do --
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;
                if (isw == 2) goto L795;
                if (alfi[en] != 0.0) goto L710;

                // Real vector
                m = en;
                b[en][en] = 1.0;
                if (na == -1) goto L800;
                alfm = alfr[m];
                betm = beta[m];

                // for i=en-1 step -1 until 1 do --
                for (ii = 0; ii <= na; ++ii)
                {
                    i = en - ii - 1;
                    w = betm * a[i][i] - alfm * b[i][i];
                    r = 0.0;

                    for (j = m; j <= en; ++j)
                        r += (betm * a[i][j] - alfm * b[i][j]) * b[j][en];

                    if (i == 0 || isw == 2)
                        goto L630;

                    if (betm * a[i][i - 1] == 0.0)
                        goto L630;

                    zz = w;
                    s = r;
                    goto L690;

                L630:
                    m = i;
                    if (isw == 2) goto L640;

                    // Real 1-by-1 block
                    t = w;
                    if (w == 0.0)
                        t = epsb;
                    b[i][en] = -r / t;
                    goto L700;

                // Real 2-by-2 block
                L640:
                    x = betm * a[i][i + 1] - alfm * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    q = w * zz - x * y;
                    t = (x * s - zz * r) / q;
                    b[i][en] = t;
                    if (Math.Abs(x) <= Math.Abs(zz)) goto L650;
                    b[i + 1][en] = (-r - w * t) / x;
                    goto L690;

                L650:
                    b[i + 1][en] = (-s - y * t) / zz;

                L690:
                    isw = 3 - isw;

                L700:
                    ;
                }
                // End real vector
                goto L800;

            // Complex vector
            L710:
                m = na;
                almr = alfr[m];
                almi = alfi[m];
                betm = beta[m];

                // last vector component chosen imaginary so that eigenvector matrix is triangular
                y = betm * a[en][na];
                b[na][na] = -almi * b[en][en] / y;
                b[na][en] = (almr * b[en][en] - betm * a[en][en]) / y;
                b[en][na] = 0.0;
                b[en][en] = 1.0;
                enm2 = na;
                if (enm2 == 0) goto L795;

                // for i=en-2 step -1 until 1 do --
                for (ii = 0; ii < enm2; ++ii)
                {
                    i = na - ii - 1;
                    w = betm * a[i][i] - almr * b[i][i];
                    w1 = -almi * b[i][i];
                    ra = 0.0;
                    sa = 0.0;

                    for (j = m; j <= en; ++j)
                    {
                        x = betm * a[i][j] - almr * b[i][j];
                        x1 = -almi * b[i][j];
                        ra = ra + x * b[j][na] - x1 * b[j][en];
                        sa = sa + x * b[j][en] + x1 * b[j][na];
                    }

                    if (i == 0 || isw == 2) goto L770;
                    if (betm * a[i][i - 1] == 0.0) goto L770;

                    zz = w;
                    z1 = w1;
                    r = ra;
                    s = sa;
                    isw = 2;
                    goto L790;

                L770:
                    m = i;
                    if (isw == 2) goto L780;

                    // Complex 1-by-1 block 
                    tr = -ra;
                    ti = -sa;

                L773:
                    dr = w;
                    di = w1;

                    // Complex divide (t1,t2) = (tr,ti) / (dr,di)
                L775:
                    if (Math.Abs(di) > Math.Abs(dr)) goto L777;
                    rr = di / dr;
                    d = dr + di * rr;
                    t1 = (tr + ti * rr) / d;
                    t2 = (ti - tr * rr) / d;

                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                L777:
                    rr = dr / di;
                    d = dr * rr + di;
                    t1 = (tr * rr + ti) / d;
                    t2 = (ti * rr - tr) / d;
                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                   // Complex 2-by-2 block 
                L780:
                    x = betm * a[i][i + 1] - almr * b[i][i + 1];
                    x1 = -almi * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    tr = y * ra - w * r + w1 * s;
                    ti = y * sa - w * s - w1 * r;
                    dr = w * zz - w1 * z1 - x * y;
                    di = w * z1 + w1 * zz - x1 * y;
                    if (dr == 0.0 && di == 0.0)
                        dr = epsb;
                    goto L775;

                L782:
                    b[i + 1][na] = t1;
                    b[i + 1][en] = t2;
                    isw = 1;
                    if (Math.Abs(y) > Math.Abs(w) + Math.Abs(w1))
                        goto L785;
                    tr = -ra - x * b[(i + 1)][na] + x1 * b[(i + 1)][en];
                    ti = -sa - x * b[(i + 1)][en] - x1 * b[(i + 1)][na];
                    goto L773;

                L785:
                    t1 = (-r - zz * b[(i + 1)][na] + z1 * b[(i + 1)][en]) / y;
                    t2 = (-s - zz * b[(i + 1)][en] - z1 * b[(i + 1)][na]) / y;

                L787:
                    b[i][na] = t1;
                    b[i][en] = t2;

                L790:
                    ;
                }

                // End complex vector
            L795:
                isw = 3 - isw;

            L800:
                ;
            }

            // End back substitution. Transform to original coordinate system.
            for (jj = 0; jj < n; ++jj)
            {
                j = n - jj - 1;

                for (i = 0; i < n; ++i)
                {
                    zz = 0.0;
                    for (k = 0; k <= j; ++k)
                        zz += z[i][k] * b[k][j];
                    z[i][j] = zz;
                }
            }

            // Normalize so that modulus of largest component of each vector is 1.
            // (isw is 1 initially from before)
            for (j = 0; j < n; ++j)
            {
                d = 0.0;
                if (isw == 2) goto L920;
                if (alfi[j] != 0.0) goto L945;

                for (i = 0; i < n; ++i)
                {
                    if ((Math.Abs(z[i][j])) > d)
                        d = (Math.Abs(z[i][j]));
                }

                for (i = 0; i < n; ++i)
                    z[i][j] /= d;

                goto L950;

            L920:
                for (i = 0; i < n; ++i)
                {
                    r = System.Math.Abs(z[i][j - 1]) + System.Math.Abs(z[i][j]);
                    if (r != 0.0)
                    {
                        // Computing 2nd power
                        double u1 = z[i][j - 1] / r;
                        double u2 = z[i][j] / r;
                        r *= Math.Sqrt(u1 * u1 + u2 * u2);
                    }
                    if (r > d)
                        d = r;
                }

                for (i = 0; i < n; ++i)
                {
                    z[i][j - 1] /= d;
                    z[i][j] /= d;
                }

            L945:
                isw = 3 - isw;

            L950:
                ;
            }

            return;
        }
        /// <summary>
        ///   Estimates unit round-off in quantities of size x.
        /// </summary>
        /// <remarks>
        ///   This is a port of the epslon function from EISPACK.
        /// </remarks>
        /// 
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Epsilon(double x)
        {
            double a, b, c, eps;

            a = 1.3333333333333333;

        L10:
            b = a - 1.0;
            c = b + b + b;
            eps = System.Math.Abs(c - 1.0);

            if (eps == 0.0)
                goto L10;

            return eps * System.Math.Abs(x);
        }
        /// <summary>
        ///   Returns <paramref name="a"/> with the sign of <paramref name="b"/>. 
        /// </summary>
        /// 
        /// <remarks>
        ///   This is a port of the sign transfer function from EISPACK,
        ///   and is is equivalent to C++'s std::copysign function.
        /// </remarks>
        /// 
        /// <returns>If B > 0 then the result is ABS(A), else it is -ABS(A).</returns>
        /// 
        /// <param name="a"></param>
        /// <param name="b"></param>
        private static double Sign(double a, double b)
        {
            double x = (a >= 0 ? a : -a);
            return (b >= 0 ? x : -x);
        }
        #endregion
    }
    /// <summary>
    /// Определяет спектральное разложение квадратной матрицы.
    /// <remarks>
    /// Спектральное разложение - это представление квадратной матрицы A в виде произведения трёх матриц A = V * D * inv(V), где V - матрица спектральных векторов, а D - диагональная (в общем виде комплексная) матрица собственных значений.
    /// Матрица A может быть также представлена в виде произведения трех матриц: A = V * R * inv(V), где R - вещественная почти диагональная матрица собственных значений.
    /// Не все матрицы могут быть представлены в таком виде, а только те, которые обладают полным набором собственных векторов. 
    /// Спектральное разложение может использоваться для нахождения собственных значений и собственных векторов матрицы, решения систем линейных уравнений, обращения матрицы, нахождения определителя матрицы и вычисления аналитических функций от матриц.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix
    /// </remarks>
    /// </summary>
    public class EVD
    {
        #region Private data
        private int n;                    // Размерность матрицы
        private double[] Re, Im;          // Спектральные значения [Re, Im]
        private double[][] matrices;      // Спектральные вектора
        private double[][] hessenberg;    // Несимметричная форма Гейзенберга
        private double[] orthogonal;      // Ортогональные вектора
        private double eps;               // Погрешность
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует спектральное разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        /// <param name="eps">Погрешность [0, 1]</param>
        public EVD(double[,] A, double eps = 1e-16)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            this.n = A.GetLength(0);
            this.Re = new double[n];
            this.Im = new double[n];
            this.eps = Maths.Double(eps);

            // for symmetric matrices eigen-value decomposition
            // without Hessenberg form.
            if (Matrice.IsSymmetric(A))
            {
                hessenberg = Jagged.Zero(n, n);
                matrices = Jagged.ToJagged(A);

                tred2(); // Tridiagonalize.
                tql2();  // Diagonalize.
            }
            // with Hessenberg form.
            else
            {
                matrices = Jagged.Zero(n, n);
                hessenberg = Jagged.ToJagged(A);
                orthogonal = new double[n];

                orthes(); // Reduce to Hessenberg form.
                hqr2();   // Reduce Hessenberg to real Schur form.
            }
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает спектральные вектора.
        /// </summary>
        public double[,] V
        {
            get { return Jagged.FromJagged(matrices); }
        }
        /// <summary>
        /// Получает собственные значения.
        /// </summary>
        public Complex[] D
        {
            get
            {
                Complex[] D = new Complex[n];

                for (int i = 0; i < n; i++)
                {
                    D[i] = new Complex(Re[i], Im[i]);
                }

                return D;
            }
        }
        /// <summary>
        /// Получает вещественную диагональную матрицу собственных значений.
        /// </summary>
        public double[,] R
        {
            get
            {
                double[,] D = new double[n, n];
                int i, j;

                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        D[i, j] = 0;
                    }

                    D[i, i] = Re[i];

                    if (Im[i] > 0)
                    {
                        D[i, i + 1] = Im[i];
                    }
                    else if (Im[i] < 0)
                    {
                        D[i, i - 1] = Im[i];
                    }
                }

                return D;
            }
        }
        /// <summary>
        /// Получает форму Хессенберга.
        /// </summary>
        public double[,] H
        {
            get { return Jagged.FromJagged(hessenberg); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Сокращение до трехдиагональной формы. Для симметричного случая.
        /// </summary>
        private void tred2()
        {
            int i, j, k;
            // Symmetric Householder reduction to tridiagonal form.
            // This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
            }

            double scale, h, f, g, hh;

            // Householder reduction to tridiagonal form.
            for (i = n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                scale = 0;
                h = 0;
                for (k = 0; k < i; k++)
                    scale = scale + System.Math.Abs(Re[k]);

                if (scale == 0)
                {
                    Im[i] = Re[i - 1];
                    for (j = 0; j < i; j++)
                    {
                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                        matrices[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder Matrice.
                    for (k = 0; k < i; k++)
                    {
                        Re[k] /= scale;
                        h += Re[k] * Re[k];
                    }

                    f = Re[i - 1];
                    g = (Double)System.Math.Sqrt(h);
                    if (f > 0) g = -g;

                    Im[i] = scale * g;
                    h = h - f * g;
                    Re[i - 1] = f - g;
                    for (j = 0; j < i; j++)
                        Im[j] = 0;

                    // Apply similarity transformation to remaining columns.
                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        matrices[j][i] = f;
                        g = Im[j] + matrices[j][j] * f;
                        for (k = j + 1; k <= i - 1; k++)
                        {
                            g += matrices[k][j] * Re[k];
                            Im[k] += matrices[k][j] * f;
                        }
                        Im[j] = g;
                    }

                    f = 0;
                    for (j = 0; j < i; j++)
                    {
                        Im[j] /= h;
                        f += Im[j] * Re[j];
                    }

                    hh = f / (h + h);
                    for (j = 0; j < i; j++)
                        Im[j] -= hh * Re[j];

                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        g = Im[j];
                        for (k = j; k <= i - 1; k++)
                            matrices[k][j] -= (f * Im[k] + g * Re[k]);

                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                    }
                }
                Re[i] = h;
            }

            // Accumulate transformations.
            for (i = 0; i < n - 1; i++)
            {
                matrices[n - 1][i] = matrices[i][i];
                matrices[i][i] = 1;
                h = Re[i + 1];
                if (h != 0)
                {
                    for (k = 0; k <= i; k++)
                        Re[k] = matrices[k][i + 1] / h;

                    for (j = 0; j <= i; j++)
                    {
                        g = 0;
                        for (k = 0; k <= i; k++)
                            g += matrices[k][i + 1] * matrices[k][j];
                        for (k = 0; k <= i; k++)
                            matrices[k][j] -= g * Re[k];
                    }
                }

                for (k = 0; k <= i; k++)
                    matrices[k][i + 1] = 0;
            }

            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
                matrices[n - 1][j] = 0;
            }

            matrices[n - 1][n - 1] = 1;
            Im[0] = 0;
        }
        /// <summary>
        /// Приведение к трехдиагональной QL-форме. Для симметричного случая.
        /// </summary>
        private void tql2()
        {
            double f = 0;
            double tst1 = 0;
            int i, l, j, k, iter, m;
            double g, p, r, dl1, h;
            double c, c2, c3, el1, s, s2;

            // Symmetric tridiagonal QL algorithm.
            // This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (i = 1; i < n; i++)
                Im[i - 1] = Im[i];

            Im[n - 1] = 0;

            for (l = 0; l < n; l++)
            {
                // Find small subdiagonal element.
                tst1 = System.Math.Max(tst1, System.Math.Abs(Re[l]) + System.Math.Abs(Im[l]));
                m = l;
                while (m < n)
                {
                    if (System.Math.Abs(Im[m]) <= eps * tst1)
                        break;
                    m++;
                }

                // If m == l, d[l] is an eigenvalue, otherwise, iterate.
                if (m > l)
                {
                    iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        g = Re[l];
                        p = (Re[l + 1] - g) / (2 * Im[l]);
                        r = Maths.Hypotenuse(p, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        Re[l] = Im[l] / (p + r);
                        Re[l + 1] = Im[l] * (p + r);
                        dl1 = Re[l + 1];
                        h = g - Re[l];
                        for (i = l + 2; i < n; i++)
                        {
                            Re[i] -= h;
                        }

                        f = f + h;

                        // Implicit QL transformation.
                        p = Re[m];
                        c = 1;
                        c2 = c;
                        c3 = c;
                        el1 = Im[l + 1];
                        s = 0;
                        s2 = 0;

                        for (i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * Im[i];
                            h = c * p;
                            r = Maths.Hypotenuse(p, Im[i]);
                            Im[i + 1] = s * r;
                            s = Im[i] / r;
                            c = p / r;
                            p = c * Re[i] - s * g;
                            Re[i + 1] = h + s * (c * g + s * Re[i]);

                            // Accumulate transformation.
                            for (k = 0; k < n; k++)
                            {
                                h = matrices[k][i + 1];
                                matrices[k][i + 1] = s * matrices[k][i] + c * h;
                                matrices[k][i] = c * matrices[k][i] - s * h;
                            }
                        }

                        p = -s * s2 * c3 * el1 * Im[l] / dl1;
                        Im[l] = s * p;
                        Re[l] = c * p;

                        // Check for convergence.
                    }
                    while (System.Math.Abs(Im[l]) > eps * tst1);
                }
                Re[l] = Re[l] + f;
                Im[l] = 0;
            }

            // Sort eigenvalues and corresponding Matrices.
            for (i = 0; i < n - 1; i++)
            {
                k = i;
                p = Re[i];
                for (j = i + 1; j < n; j++)
                {
                    if (Re[j] < p)
                    {
                        k = j;
                        p = Re[j];
                    }
                }

                if (k != i)
                {
                    Re[k] = Re[i];
                    Re[i] = p;
                    for (j = 0; j < n; j++)
                    {
                        p = matrices[j][i];
                        matrices[j][i] = matrices[j][k];
                        matrices[j][k] = p;
                    }
                }
            }
        }
        /// <summary>
        /// Несимметричное сокращение до формы Хессенберга.
        /// </summary>
        private void orthes()
        {
            // Nonsymmetric reduction to Hessenberg form.
            // This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
            int low = 0;
            int high = n - 1;
            int m, i, j;
            double scale, h, g, f;

            for (m = low + 1; m <= high - 1; m++)
            {
                // Scale column.

                scale = 0;
                for (i = m; i <= high; i++)
                    scale = scale + System.Math.Abs(hessenberg[i][m - 1]);

                if (scale != 0)
                {
                    // Compute Householder transformation.
                    h = 0;
                    for (i = high; i >= m; i--)
                    {
                        orthogonal[i] = hessenberg[i][m - 1] / scale;
                        h += orthogonal[i] * orthogonal[i];
                    }

                    g = (Double)System.Math.Sqrt(h);
                    if (orthogonal[m] > 0) g = -g;

                    h = h - orthogonal[m] * g;
                    orthogonal[m] = orthogonal[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u') / h)
                    for (j = m; j < n; j++)
                    {
                        f = 0;
                        for (i = high; i >= m; i--)
                            f += orthogonal[i] * hessenberg[i][j];

                        f = f / h;
                        for (i = m; i <= high; i++)
                            hessenberg[i][j] -= f * orthogonal[i];
                    }

                    for (i = 0; i <= high; i++)
                    {
                        f = 0;
                        for (j = high; j >= m; j--)
                            f += orthogonal[j] * hessenberg[i][j];

                        f = f / h;
                        for (j = m; j <= high; j++)
                            hessenberg[i][j] -= f * orthogonal[j];
                    }

                    orthogonal[m] = scale * orthogonal[m];
                    hessenberg[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    matrices[i][j] = (i == j ? 1 : 0);

            for (m = high - 1; m >= low + 1; m--)
            {
                if (hessenberg[m][m - 1] != 0)
                {
                    for (i = m + 1; i <= high; i++)
                        orthogonal[i] = hessenberg[i][m - 1];

                    for (j = m; j <= high; j++)
                    {
                        g = 0;
                        for (i = m; i <= high; i++)
                            g += orthogonal[i] * matrices[i][j];

                        // Double division avoids possible underflow.
                        g = (g / orthogonal[m]) / hessenberg[m][m - 1];
                        for (i = m; i <= high; i++)
                            matrices[i][j] += g * orthogonal[i];
                    }
                }
            }
        }
        /// <summary>
        /// Сокращение Хессенберга до реальной формы Шура.
        /// </summary>
        private void hqr2()
        {
            // Nonsymmetric reduction from Hessenberg to real Schur form.   
            // This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
            // Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
            int nn = this.n;
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            //double eps = 2 * double.Epsilon;
            double exshift = 0;
            double p = 0;
            double q = 0;
            double r = 0;
            double s = 0;
            double z = 0;
            double t;
            double w;
            double x;
            double y;
            int i, j, k, m;
            bool notlast;

            // Store roots isolated by balanc and compute matrix norm
            double norm = 0;
            for (i = 0; i < nn; i++)
            {
                if (i < low | i > high)
                {
                    Re[i] = hessenberg[i][i];
                    Im[i] = 0;
                }

                for (j = System.Math.Max(i - 1, 0); j < nn; j++)
                    norm = norm + System.Math.Abs(hessenberg[i][j]);
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n >= low)
            {
                // Look for single small sub-diagonal element
                int l = n;
                while (l > low)
                {
                    s = System.Math.Abs(hessenberg[l - 1][l - 1]) + System.Math.Abs(hessenberg[l][l]);

                    if (s == 0)
                        s = norm;

                    if (double.IsNaN(s))
                        break;

                    if (System.Math.Abs(hessenberg[l][l - 1]) < eps * s)
                        break;

                    l--;
                }

                // Check for convergence
                if (l == n)
                {
                    // One root found
                    hessenberg[n][n] = hessenberg[n][n] + exshift;
                    Re[n] = hessenberg[n][n];
                    Im[n] = 0;
                    n--;
                    iter = 0;
                }
                else if (l == n - 1)
                {
                    // Two roots found
                    w = hessenberg[n][n - 1] * hessenberg[n - 1][n];
                    p = (hessenberg[n - 1][n - 1] - hessenberg[n][n]) / 2;
                    q = p * p + w;
                    z = (double)System.Math.Sqrt(System.Math.Abs(q));
                    hessenberg[n][n] = hessenberg[n][n] + exshift;
                    hessenberg[n - 1][n - 1] = hessenberg[n - 1][n - 1] + exshift;
                    x = hessenberg[n][n];

                    if (q >= 0)
                    {
                        // Real pair
                        z = (p >= 0) ? (p + z) : (p - z);
                        Re[n - 1] = x + z;
                        Re[n] = Re[n - 1];
                        if (z != 0)
                            Re[n] = x - w / z;
                        Im[n - 1] = 0;
                        Im[n] = 0;
                        x = hessenberg[n][n - 1];
                        s = System.Math.Abs(x) + System.Math.Abs(z);
                        p = x / s;
                        q = z / s;
                        r = (Double)System.Math.Sqrt(p * p + q * q);
                        p = p / r;
                        q = q / r;

                        // Row modification
                        for (j = n - 1; j < nn; j++)
                        {
                            z = hessenberg[n - 1][j];
                            hessenberg[n - 1][j] = q * z + p * hessenberg[n][j];
                            hessenberg[n][j] = q * hessenberg[n][j] - p * z;
                        }

                        // Column modification
                        for (i = 0; i <= n; i++)
                        {
                            z = hessenberg[i][n - 1];
                            hessenberg[i][n - 1] = q * z + p * hessenberg[i][n];
                            hessenberg[i][n] = q * hessenberg[i][n] - p * z;
                        }

                        // Accumulate transformations
                        for (i = low; i <= high; i++)
                        {
                            z = matrices[i][n - 1];
                            matrices[i][n - 1] = q * z + p * matrices[i][n];
                            matrices[i][n] = q * matrices[i][n] - p * z;
                        }
                    }
                    else
                    {
                        // Complex pair
                        Re[n - 1] = x + p;
                        Re[n] = x + p;
                        Im[n - 1] = z;
                        Im[n] = -z;
                    }

                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet     

                    // Form shift
                    x = hessenberg[n][n];
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = hessenberg[n - 1][n - 1];
                        w = hessenberg[n][n - 1] * hessenberg[n - 1][n];
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10)
                    {
                        exshift += x;
                        for (i = low; i <= n; i++)
                            hessenberg[i][i] -= x;

                        s = System.Math.Abs(hessenberg[n][n - 1]) + System.Math.Abs(hessenberg[n - 1][n - 2]);
                        x = y = (double)0.75 * s;
                        w = (double)(-0.4375) * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = (double)System.Math.Sqrt(s);
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (i = low; i <= n; i++)
                                hessenberg[i][i] -= s;
                            exshift += s;
                            x = y = w = (Double)0.964;
                        }
                    }

                    iter = iter + 1;

                    // Look for two consecutive small sub-diagonal elements
                    m = n - 2;
                    while (m >= l)
                    {
                        z = hessenberg[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / hessenberg[m + 1][m] + hessenberg[m][m + 1];
                        q = hessenberg[m + 1][m + 1] - z - r - s;
                        r = hessenberg[m + 2][m + 1];
                        s = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                        p = p / s;
                        q = q / s;
                        r = r / s;
                        if (m == l)
                            break;
                        if (System.Math.Abs(hessenberg[m][m - 1]) * (System.Math.Abs(q) + System.Math.Abs(r)) < eps * (System.Math.Abs(p) * (System.Math.Abs(hessenberg[m - 1][m - 1]) + System.Math.Abs(z) + System.Math.Abs(hessenberg[m + 1][m + 1]))))
                            break;
                        m--;
                    }

                    for (i = m + 2; i <= n; i++)
                    {
                        hessenberg[i][i - 2] = 0;
                        if (i > m + 2)
                            hessenberg[i][i - 3] = 0;
                    }

                    // Double QR step involving rows l:n and columns m:n
                    for (k = m; k <= n - 1; k++)
                    {
                        notlast = (k != n - 1);
                        if (k != m)
                        {
                            p = hessenberg[k][k - 1];
                            q = hessenberg[k + 1][k - 1];
                            r = (notlast ? hessenberg[k + 2][k - 1] : 0);
                            x = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                            if (x != 0)
                            {
                                p = p / x;
                                q = q / x;
                                r = r / x;
                            }
                        }

                        if (x == 0) break;

                        s = (Double)System.Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0) s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                hessenberg[k][k - 1] = -s * x;
                            else
                                if (l != m)
                                    hessenberg[k][k - 1] = -hessenberg[k][k - 1];

                            p = p + s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q = q / p;
                            r = r / p;

                            // Row modification
                            for (j = k; j < nn; j++)
                            {
                                p = hessenberg[k][j] + q * hessenberg[k + 1][j];
                                if (notlast)
                                {
                                    p = p + r * hessenberg[k + 2][j];
                                    hessenberg[k + 2][j] = hessenberg[k + 2][j] - p * z;
                                }

                                hessenberg[k][j] = hessenberg[k][j] - p * x;
                                hessenberg[k + 1][j] = hessenberg[k + 1][j] - p * y;
                            }

                            // Column modification
                            for (i = 0; i <= System.Math.Min(n, k + 3); i++)
                            {
                                p = x * hessenberg[i][k] + y * hessenberg[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * hessenberg[i][k + 2];
                                    hessenberg[i][k + 2] = hessenberg[i][k + 2] - p * r;
                                }

                                hessenberg[i][k] = hessenberg[i][k] - p;
                                hessenberg[i][k + 1] = hessenberg[i][k + 1] - p * q;
                            }

                            // Accumulate transformations
                            for (i = low; i <= high; i++)
                            {
                                p = x * matrices[i][k] + y * matrices[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * matrices[i][k + 2];
                                    matrices[i][k + 2] = matrices[i][k + 2] - p * r;
                                }

                                matrices[i][k] = matrices[i][k] - p;
                                matrices[i][k + 1] = matrices[i][k + 1] - p * q;
                            }
                        }
                    }
                }
            }

            // Backsubstitute to find Matrices of upper triangular form
            if (norm == 0)
            {
                return;
            }

            for (n = nn - 1; n >= 0; n--)
            {
                p = Re[n];
                q = Im[n];

                // Real Matrice
                if (q == 0)
                {
                    int l = n;
                    hessenberg[n][n] = 1;
                    for (i = n - 1; i >= 0; i--)
                    {
                        w = hessenberg[i][i] - p;
                        r = 0;
                        for (j = l; j <= n; j++)
                            r = r + hessenberg[i][j] * hessenberg[j][n];

                        if (Im[i] < 0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (Im[i] == 0)
                            {
                                hessenberg[i][n] = (w != 0) ? (-r / w) : (-r / (eps * norm));
                            }
                            else
                            {
                                // Solve real equations
                                x = hessenberg[i][i + 1];
                                y = hessenberg[i + 1][i];
                                q = (Re[i] - p) * (Re[i] - p) + Im[i] * Im[i];
                                t = (x * s - z * r) / q;
                                hessenberg[i][n] = t;
                                hessenberg[i + 1][n] = (System.Math.Abs(x) > System.Math.Abs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                            }

                            // Overflow control
                            t = System.Math.Abs(hessenberg[i][n]);
                            if ((eps * t) * t > 1)
                                for (j = i; j <= n; j++)
                                    hessenberg[j][n] = hessenberg[j][n] / t;
                        }
                    }
                }
                else if (q < 0)
                {
                    // Complex Matrice
                    int l = n - 1;

                    // Last Matrice component imaginary so matrix is triangular
                    if (System.Math.Abs(hessenberg[n][n - 1]) > System.Math.Abs(hessenberg[n - 1][n]))
                    {
                        hessenberg[n - 1][n - 1] = q / hessenberg[n][n - 1];
                        hessenberg[n - 1][n] = -(hessenberg[n][n] - p) / hessenberg[n][n - 1];
                    }
                    else
                    {
                        cdiv(0, -hessenberg[n - 1][n], hessenberg[n - 1][n - 1] - p, q, ref hessenberg[n - 1][n - 1], ref hessenberg[n - 1][n]);
                    }

                    hessenberg[n][n - 1] = 0;
                    hessenberg[n][n] = 1;
                    for (i = n - 2; i >= 0; i--)
                    {
                        double ra, sa, vr, vi;
                        ra = 0;
                        sa = 0;
                        for (j = l; j <= n; j++)
                        {
                            ra = ra + hessenberg[i][j] * hessenberg[j][n - 1];
                            sa = sa + hessenberg[i][j] * hessenberg[j][n];
                        }

                        w = hessenberg[i][i] - p;

                        if (Im[i] < 0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (Im[i] == 0)
                            {
                                cdiv(-ra, -sa, w, q, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
                            }
                            else
                            {
                                // Solve complex equations
                                x = hessenberg[i][i + 1];
                                y = hessenberg[i + 1][i];
                                vr = (Re[i] - p) * (Re[i] - p) + Im[i] * Im[i] - q * q;
                                vi = (Re[i] - p) * 2 * q;
                                if (vr == 0 & vi == 0)
                                    vr = eps * norm * (System.Math.Abs(w) + System.Math.Abs(q) + System.Math.Abs(x) + System.Math.Abs(y) + System.Math.Abs(z));
                                cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, ref hessenberg[i][n - 1], ref hessenberg[i][n]);
                                if (System.Math.Abs(x) > (System.Math.Abs(z) + System.Math.Abs(q)))
                                {
                                    hessenberg[i + 1][n - 1] = (-ra - w * hessenberg[i][n - 1] + q * hessenberg[i][n]) / x;
                                    hessenberg[i + 1][n] = (-sa - w * hessenberg[i][n] - q * hessenberg[i][n - 1]) / x;
                                }
                                else
                                {
                                    cdiv(-r - y * hessenberg[i][n - 1], -s - y * hessenberg[i][n], z, q, ref hessenberg[i + 1][n - 1], ref hessenberg[i + 1][n]);
                                }
                            }

                            // Overflow control
                            t = System.Math.Max(System.Math.Abs(hessenberg[i][n - 1]), System.Math.Abs(hessenberg[i][n]));
                            if ((eps * t) * t > 1)
                            {
                                for (j = i; j <= n; j++)
                                {
                                    hessenberg[j][n - 1] = hessenberg[j][n - 1] / t;
                                    hessenberg[j][n] = hessenberg[j][n] / t;
                                }
                            }
                        }
                    }
                }
            }

            // Matrices of isolated roots
            for (i = 0; i < nn; i++)
                if (i < low | i > high)
                    for (j = i; j < nn; j++)
                        matrices[i][j] = hessenberg[i][j];

            // Back transformation to get eigenMatrices of original matrix
            for (j = nn - 1; j >= low; j--)
            {
                for (i = low; i <= high; i++)
                {
                    z = 0;
                    for (k = low; k <= System.Math.Min(j, high); k++)
                        z = z + matrices[i][k] * hessenberg[k][j];
                    matrices[i][j] = z;
                }
            }
        }
        /// <summary>
        /// Скалярное деление комплексных чисел.
        /// </summary>
        /// <param name="xr"></param>
        /// <param name="xi"></param>
        /// <param name="yr"></param>
        /// <param name="yi"></param>
        /// <param name="cdivr"></param>
        /// <param name="cdivi"></param>
        private static void cdiv(double xr, double xi, double yr, double yi, ref double cdivr, ref double cdivi)
        {
            // Complex scalar division.
            double r;
            double d;

            if (System.Math.Abs(yr) > System.Math.Abs(yi))
            {
                r = yi / yr;
                d = yr + r * yi;
                cdivr = (xr + r * xi) / d;
                cdivi = (xi - r * xr) / d;
            }
            else
            {
                r = yr / yi;
                d = yi + r * yr;
                cdivr = (r * xr + xi) / d;
                cdivi = (r * xi - xr) / d;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет LDU-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление квадратной матрицы A в виде произведения трех матриц: A = L * D * U, где L - нижняя треугольная матрица, D - диагональная матрица, а U - верхняя треугольная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/LU_decomposition
    /// </remarks>
    /// </summary>
    public class LDU
    {
        #region Private data
        private LU ludecomp;
        private Diagonal diagdecomp;
        private double[,] lower;
        private double[,] upper;
        private double[]     diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует LDU-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public LDU(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // LDU algorithm:
            // LU-decomposition:
            ludecomp = new LU(A);
            lower = ludecomp.L;
            upper = ludecomp.U;

            // Diagonal decomposition:
            diagdecomp = new Diagonal(lower);
            lower = diagdecomp.B;
            diag = diagdecomp.D;
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает нижнюю треугольную матрицу.
        /// </summary>
        public double[,] L
        {
            get { return lower; }
        }
        /// <summary>
        /// Получает верхнюю треугольную матрицу.
        /// </summary>
        public double[,] U
        {
            get { return upper; }
        }
        /// <summary>
        /// Получает вектор элментов диагонали.
        /// </summary>
        public double[] D
        {
            get { return diag; }
        }
        #endregion
    }
    /// <summary>
    /// Определяет диагональное разложение квадратной матрицы.
    /// <remarks>
    /// Это представление квадратной матрицы A в виде произведения двух матриц: A = B * D, где B - квадратная матрица, а D - диагональная матрица.
    /// Данное разложение используется для выделения диагональных матриц в других разложениях (например, LDU-, LDL- разложения).
    /// </remarks>
    /// </summary>
    public class Diagonal
    {
        #region Private data
        private double[,] matrix;
        private double[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует диагональное разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public Diagonal(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            int n = A.GetLength(0), i;
            this.diag = new double[n];

            for (i = 0; i < n; i++)
                diag[i] = A[i, i];

            this.matrix = Matrice.Dot(A, diag, true);
            return;
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает квадратную матрицу.
        /// </summary>
        public double[,] B
        {
            get { return matrix; }
        }
        /// <summary>
        /// Получает вектор элментов диагонали.
        /// </summary>
        public double[] D
        {
            get { return diag; }
        }
        #endregion
    }
    /// <summary>
    /// Определяет LU-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление квадратной матрицы A в виде произведения двух матриц: A = L * U, где L - нижняя треугольная матрица, U - верхняя треугольная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/LU_decomposition
    /// </remarks>
    /// </summary>
    public class LU
    {
        #region Private data
        private double[][] lower;
        private double[][] upper;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует LU-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public LU(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // LU-decomposition:
            ludecomp(A.ToJagged());
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает нижнюю треугольную матрицу.
        /// </summary>
        public double[,] L
        {
            get { return Jagged.FromJagged(lower); }
        }
        /// <summary>
        /// Получает верхнюю треугольную матрицу.
        /// </summary>
        public double[,] U
        {
            get { return Jagged.FromJagged(upper); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует LU-разложение квадратной матрицы.
        /// </summary>
        /// <param name="a">Квадратная матрица</param>
        private void ludecomp(double[][] a)
        {
            int i, j, k;
            int n = a.GetLength(0);
            double alpha, beta;
            this.upper = Jagged.Zero(n, n);
            this.lower = Jagged.Zero(n, n);

            for (i = 0; i < n; i++)
            {
                this.upper[i][i] = 1;
            }

            for (j = 0; j < n; j++)
            {
                for (i = j; i < n; i++)
                {
                    alpha = 0;
                    for (k = 0; k < j; k++)
                    {
                        alpha = alpha + this.lower[i][k] * this.upper[k][j];
                    }
                    this.lower[i][j] = a[i][j] - alpha;
                }

                beta = lower[j][j];

                for (i = j; i < n; i++)
                {
                    alpha = 0;
                    for (k = 0; k < j; k++)
                    {
                        alpha = alpha + this.lower[j][k] * this.upper[k][i];
                    }

                    if (beta != 0)
                    {
                        this.upper[j][i] = (a[j][i] - alpha) / beta;
                    }
                }
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет разложение Холецкого квадратной матрицы.
    /// <remarks>
    /// Это представление симметричной положительно-определённой квадратной матрицы в виде произведения: A = L * L', где L - нижняя треугольная матрица со строго положительными элементами на диагонали. 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Cholesky_decomposition
    /// </remarks>
    /// </summary>
    public class Cholesky
    {
        #region Private data
        double[][] lower;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует разложение Холецкого квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная симметричная положительно-определенная матрица</param>
        public Cholesky(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // Cholesky decomposition:
            chol(A.ToJagged());
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает нижнюю треугольную матрицу L.
        /// </summary>
        public double[,] L
        {
            get { return Jagged.FromJagged(lower); }
        }
        /// <summary>
        /// Получает верхнюю треугольную матрицу U.
        /// </summary>
        public double[,] U
        {
            get { return Matrice.Transponate(L); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует разложение Холецкого.
        /// </summary>
        /// <param name="a">Квадратная матрица</param>
        private void chol(double[][] a)
        {
            // Cholesky decomposition
            int n = a.GetLength(0);
            this.lower = new double[n][];
            double[] v, w, z, d = new double[n];
            double alpha;
            int j, i, k;

            // get diagonal elements
            for (i = 0; i < n; i++)
            {
                d[i] = a[i][i];
            }

            // do job
            for (j = 0; j < n; j++)
            {
                v = lower[j] = new double[n];
                z = a[j];

                for (i = 0; i <= j; i++)
                {
                    w = lower[i];
                    alpha = 0;

                    if (i == j)
                    {
                        for (k = 0; k < i; k++)
                        {
                            alpha += w[k] * w[k];
                        }

                        w[i] = Math.Sqrt(d[i] - alpha);
                        lower[i] = w;
                    }
                    else
                    {
                        for (k = 0; k < i; k++)
                        {
                            alpha += w[k] * v[k];
                        }

                        v[i] = (z[i] - alpha) / w[i];
                    }
                }

                lower[j] = v;
            }
            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет UDL-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление симметричной квадратной матрицы в виде произведения трех матриц: A = U * D * L, где U - верхняя треугольная матрица, D - диагональная матрица, а L - нижняя треугольная матрица.
    /// Данное разложение представляет собой специфичную форму разложения Холецкого.
    /// </remarks>
    /// </summary>
    public class UDL
    {
        #region Private data
        private double[,] upper;
        private double[] diag;
        #endregion

        #region UDL components
        /// <summary>
        /// Инициализирует UDL-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная симметричная матрица</param>
        public UDL(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            udldecomp(A);
        }
        /// <summary>
        /// Возвращает верхнюю треугольную матрицу.
        /// </summary>
        public double[,] U
        {
            get
            {
                return this.upper;
            }
        }
        /// <summary>
        /// Возвращает диагональную матрицу.
        /// </summary>
        public double[] D
        {
            get
            {
                return this.diag;
            }
        }
        /// <summary>
        /// Возвращает нижнюю треугольную матрицу.
        /// </summary>
        public double[,] L
        {
            get
            {
                return this.upper.Transponate();
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// UDL-разложение.
        /// </summary>
        /// <param name="a">Симметричная матрица</param>
        private void udldecomp(double[,] a)
        {
            int i, j, k;
            int n = a.GetLength(0);
            this.upper = new double[n, n];
            this.diag = new double[n];
            double[][] p = Jagged.ToJagged(a);
            double alpha, beta, gamma;

            // Mathematics in science and engineering, v.128,
            // Factorization methods for discrete sequential estimation, Gerald J. Bierman.
            // UDU* factorization aglorithm.
            // 
            for (j = n - 1; j >= 1; j--)
            {
                gamma = p[j][j];
                diag[j] = gamma;
                alpha = 1.0 / gamma;

                for (k = 0; k < j; k++)
                {
                    beta = p[k][j];
                    upper[k, j] = alpha * beta;

                    for (i = 0; i <= k; i++)
                    {
                        p[i][k] -= beta * upper[i, j];
                    }
                }
            }
            diag[0] = p[0][0];

            // diagonal eyes:
            for (i = 0; i < n; i++)
            {
                upper[i, i] = 1.0;
            }
            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет LDL'-разложение квадратной матрицы.
    /// Это представление симметричной положительно-определённой квадратной матрицы в виде произведения трех матриц: A = L * D * L', где L - нижняя треугольная матрица со строго положительными элементами на диагонали, 
    /// а D - диагональная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2
    /// <remarks></remarks>
    /// </summary>
    public class LDL
    {
        #region Private data
        private Cholesky choldecomp;
        private Diagonal diagdecomp;
        private double[,] lower;
        private double[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует LDL-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная симметричная положительно-определенная матрица</param>
        public LDL(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // LDL'-decomposition algorithm
            // Cholesky decomposition:
            choldecomp = new Cholesky(A);
            lower = choldecomp.L;
            
            // Diagonal decomposition:
            diagdecomp = new Diagonal(lower);
            lower = diagdecomp.B;
            diag = diagdecomp.D;

            // D = d^2:
            diag = Matrice.Mul(diag, diag);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает нижнюю треугольную матрицу L.
        /// </summary>
        public double[,] L
        {
            get { return lower; }
        }
        /// <summary>
        /// Получает верхнюю треугольную матрицу U.
        /// </summary>
        public double[,] U
        {
            get { return Matrice.Transponate(lower); }
        }
        /// <summary>
        /// Получает диагональную матрицу.
        /// </summary>
        public double[] D
        {
            get { return diag; }
        }
        #endregion
    }
    /// <summary>
    /// Определяет сингулярное разложение матрицы. 
    /// <remarks>
    /// Это представление прямоугольной матрицы A в виде произведения трёх матриц A = U * S * V', где U - левые вектора, V - правые вектора, а S - сингулярные значения.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Singular_value_decomposition
    /// </remarks>
    /// </summary>
    public class SVD
    {
        #region Private data
        private int n, m;
        private int iterations;
        private double[][] Ur;
        private double[][] Vr;
        private double[] Sr;
        private bool reversed;
        #endregion

        #region Initialize
        /// <summary>
        /// Реализует сигнулярное разложение матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="iterations">Количество итераций</param>
        public SVD(double[,] A, int iterations = 10)
        {
            // set:
            this.iterations = iterations;
            this.n = A.GetLength(0);
            this.m = A.GetLength(1);

            // options:
            if (n < m)
            {
                this.reversed = true;
                this.n = A.GetLength(1);
                this.m = A.GetLength(0);
                this.svdcmp(A.Transponate());
            }
            else
            {
                this.reversed = false;
                this.svdcmp(A);
            }
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает левые вектора.
        /// </summary>
        public double[,] U
        {
            get 
            {
                return reversed ? Jagged.FromJagged(Vr) : Jagged.FromJagged(Ur); 
            }
        }
        /// <summary>
        /// Получает сингулярные значения.
        /// </summary>
        public double[] S
        {
            get { return Sr; }
        }
        /// <summary>
        /// Получает правые вектора.
        /// </summary>
        public double[,] V
        {
            get
            {
                return reversed ? Jagged.FromJagged(Ur) : Jagged.FromJagged(Vr);
            }
        }
        /// <summary>
        /// Получает псевдообратную матрицу к исходной.
        /// </summary>
        public double[,] P
        {
            get 
            {
                // Moore–Penrose inverse:
                // P = V * (I / S) * U'
                return V.Dot(Matrice.One(m).Div(S)).Dot(U.Transponate());
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует сингулярное разложение матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        private void svdcmp(double[,] A)
        {
            this.Ur = Jagged.ToJagged(A);
            this.Sr = new double[m];
            this.Vr = Jagged.Zero(m, m);
            double[] rv1 = new double[m];

            int flag, i, its, j, jj, k, l = 0, nm = 0;
            double anorm, c, f, g, h, e, scale, x, y, z;


            // householder reduction to bidiagonal form
            g = scale = anorm = 0.0;

            for (i = 0; i < m; i++)
            {
                l = i + 1;
                rv1[i] = scale * g;
                g = e = scale = 0;

                if (i < n)
                {
                    for (k = i; k < n; k++)
                    {
                        scale += Math.Abs(Ur[k][i]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = i; k < n; k++)
                        {
                            Ur[k][i] /= scale;
                            e += Ur[k][i] * Ur[k][i];
                        }

                        f = Ur[i][i];
                        g = -Sign(Math.Sqrt(e), f);
                        h = f * g - e;
                        Ur[i][i] = f - g;

                        if (i != m - 1)
                        {
                            for (j = l; j < m; j++)
                            {
                                for (e = 0.0, k = i; k < n; k++)
                                {
                                    e += Ur[k][i] * Ur[k][j];
                                }

                                f = e / h;

                                for (k = i; k < n; k++)
                                {
                                    Ur[k][j] += f * Ur[k][i];
                                }
                            }
                        }

                        for (k = i; k < n; k++)
                        {
                            Ur[k][i] *= scale;
                        }
                    }
                }

                Sr[i] = scale * g;
                g = e = scale = 0.0;

                if ((i < n) && (i != m - 1))
                {
                    for (k = l; k < m; k++)
                    {
                        scale += Math.Abs(Ur[i][k]);
                    }

                    if (scale != 0.0)
                    {
                        for (k = l; k < m; k++)
                        {
                            Ur[i][k] /= scale;
                            e += Ur[i][k] * Ur[i][k];
                        }

                        f = Ur[i][l];
                        g = -Sign(Math.Sqrt(e), f);
                        h = f * g - e;
                        Ur[i][l] = f - g;

                        for (k = l; k < m; k++)
                        {
                            rv1[k] = Ur[i][k] / h;
                        }

                        if (i != n - 1)
                        {
                            for (j = l; j < n; j++)
                            {
                                for (e = 0.0, k = l; k < m; k++)
                                {
                                    e += Ur[j][k] * Ur[i][k];
                                }
                                for (k = l; k < m; k++)
                                {
                                    Ur[j][k] += e * rv1[k];
                                }
                            }
                        }

                        for (k = l; k < m; k++)
                        {
                            Ur[i][k] *= scale;
                        }
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(Sr[i]) + Math.Abs(rv1[i])));
            }

            // accumulation of right-hand transformations
            for (i = m - 1; i >= 0; i--)
            {
                if (i < m - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < m; j++)
                        {
                            Vr[j][i] = (Ur[i][j] / Ur[i][l]) / g;
                        }

                        for (j = l; j < m; j++)
                        {
                            for (e = 0, k = l; k < m; k++)
                            {
                                e += Ur[i][k] * Vr[k][j];
                            }
                            for (k = l; k < m; k++)
                            {
                                Vr[k][j] += e * Vr[k][i];
                            }
                        }
                    }
                    for (j = l; j < m; j++)
                    {
                        Vr[i][j] = Vr[j][i] = 0;
                    }
                }
                Vr[i][i] = 1;
                g = rv1[i];
                l = i;
            }

            // accumulation of left-hand transformations
            for (i = m - 1; i >= 0; i--)
            {
                l = i + 1;
                g = Sr[i];

                if (i < m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        Ur[i][j] = 0.0;
                    }
                }

                if (g != 0)
                {
                    g = 1.0 / g;

                    if (i != m - 1)
                    {
                        for (j = l; j < m; j++)
                        {
                            for (e = 0, k = l; k < n; k++)
                            {
                                e += Ur[k][i] * Ur[k][j];
                            }

                            f = (e / Ur[i][i]) * g;

                            for (k = i; k < n; k++)
                            {
                                Ur[k][j] += f * Ur[k][i];
                            }
                        }
                    }

                    for (j = i; j < n; j++)
                    {
                        Ur[j][i] *= g;
                    }
                }
                else
                {
                    for (j = i; j < n; j++)
                    {
                        Ur[j][i] = 0;
                    }
                }
                ++Ur[i][i];
            }

            // diagonalization of the bidiagonal form: Loop over singular values
            // and over allowed iterations
            for (k = m - 1; k >= 0; k--)
            {
                for (its = 1; its <= iterations; its++)
                {
                    flag = 1;

                    for (l = k; l >= 0; l--)
                    {
                        // test for splitting
                        nm = l - 1;

                        if (Math.Abs(rv1[l]) + anorm == anorm)
                        {
                            flag = 0;
                            break;
                        }

                        if (Math.Abs(Sr[nm]) + anorm == anorm)
                            break;
                    }

                    if (flag != 0)
                    {
                        c = 0.0;
                        e = 1.0;
                        for (i = l; i <= k; i++)
                        {
                            f = e * rv1[i];

                            if (Math.Abs(f) + anorm != anorm)
                            {
                                g = Sr[i];
                                h = Maths.Hypotenuse(f, g);
                                Sr[i] = h;
                                h = 1.0 / h;
                                c = g * h;
                                e = -f * h;

                                //for (j = 1; j <= m; j++)
                                for (j = 1; j < n; j++)
                                {
                                    y = Ur[j][nm];
                                    z = Ur[j][i];
                                    Ur[j][nm] = y * c + z * e;
                                    Ur[j][i] = z * c - y * e;
                                }
                            }
                        }
                    }

                    z = Sr[k];

                    if (l == k)
                    {
                        // convergence
                        if (z < 0.0)
                        {
                            // singular value is made nonnegative
                            Sr[k] = -z;

                            for (j = 0; j < m; j++)
                            {
                                Vr[j][k] = -Vr[j][k];
                            }
                        }
                        break;
                    }

                    if (its == iterations)
                    {
                        throw new ApplicationException("Нет сходимости в " + iterations.ToString() + " итерациях сингулярного разложения");
                    }

                    // shift from bottom 2-by-2 minor
                    x = Sr[l];
                    nm = k - 1;
                    y = Sr[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = Maths.Hypotenuse(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + Sign(g, f))) - h)) / x;

                    // next QR transformation
                    c = e = 1.0;

                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = Sr[i];
                        h = e * g;
                        g = c * g;
                        z = Maths.Hypotenuse(f, h);
                        rv1[j] = z;
                        c = f / z;
                        e = h / z;
                        f = x * c + g * e;
                        g = g * c - x * e;
                        h = y * e;
                        y *= c;

                        for (jj = 0; jj < m; jj++)
                        {
                            x = Vr[jj][j];
                            z = Vr[jj][i];
                            Vr[jj][j] = x * c + z * e;
                            Vr[jj][i] = z * c - x * e;
                        }

                        z = Maths.Hypotenuse(f, h);
                        Sr[j] = z;

                        if (z != 0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            e = h * z;
                        }

                        f = c * g + e * y;
                        x = c * y - e * g;

                        for (jj = 0; jj < n; jj++)
                        {
                            y = Ur[jj][j];
                            z = Ur[jj][i];
                            Ur[jj][j] = y * c + z * e;
                            Ur[jj][i] = z * c - y * e;
                        }
                    }

                    rv1[l] = 0.0;
                    rv1[k] = f;
                    Sr[k] = x;
                }
            }
        }
        /// <summary>
        /// Определение знака числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Знак</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private static double Sign(double a, double b)
        {
            return (b >= 0.0) ? System.Math.Abs(a) : -System.Math.Abs(a);
        }
        #endregion
    }
    /// <summary>
    /// Определяет полярное разложение матрицы.
    /// <remarks>
    /// Это представление прямоугольной матрицы A в виде произведения двух матриц: A = U * P, где U - унитарная матрица, P - положительно-определенная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Polar_decomposition
    /// </remarks>
    /// </summary>
    public class Polar
    {
        #region Private data
        private SVD svd;
        double[,] u;
        double[,] p;
        #endregion

        #region Initialize
        /// <summary>
        /// Реализует полярное разложение матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="iterations">Количество итераций</param>
        public Polar(double[,] A, int iterations = 10)
        {
            // Сингулярное разложение матрицы:
            svd = new SVD(A, iterations);
            double[,] U = svd.U, V = svd.V, H = V.Transponate();
            double[] S = svd.S;

            // Определение матриц U и P:
            u = U.Dot(H); p = V.Dot(S).Dot(H);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает унитарную матрицу.
        /// </summary>
        public double[,] U
        {
            get
            {
                return this.u;
            }
        }
        /// <summary>
        /// Получает положительно-определенную матрицу
        /// </summary>
        public double[,] P
        {
            get
            {
                return this.p;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет неотрицательную матричную факторизацию.
    /// <remarks>
    /// Это представление прямоугольной матрицы A в виде произведения двух матриц: A = W * H.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Non-negative_matrix_factorization
    /// </remarks>
    /// </summary>
    public class NMF
    {
        #region Private data
        private double[,] w;  // W is m x r (weights)
        private double[,] h;  // H is r x n (transformed data) (transposed)
        private int n;   // number of input data vectors
        private int m;   // dimension of input vector
        private int r;   // dimension of output vector (reduced dimension)
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует неотрицательную матричную факторизацию.
        /// </summary>
        /// <param name="A">Неотрицательная матрица</param>
        /// <param name="r">Размерность новых матриц</param>
        /// <param name="iterations">Количество итераций</param>
        public NMF(double[,] A, int r, int iterations = 100)
        {
            this.m = A.GetLength(0);
            this.n = A.GetLength(1);

            if (n < m)
                throw new Exception("Ширина матрицы должна быть больше высоты");

            this.r = r;

            // decompose
            nnmf(A, iterations);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает левую матрицу.
        /// </summary>
        public double[,] W
        {
            get { return w; }
        }
        /// <summary>
        /// Получает правую матрицу.
        /// </summary>
        public double[,] H
        {
            get { return h; }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Представляет NNMF, основанный на мультипликативном методе.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="iterations">Количество итераций</param>
        private void nnmf(double[,] A, int iterations)
        {
            // chose W and H randomly, W with unit norm
            w = Matrice.Rand(m, r);
            h = Matrice.Rand(r, n);
            var Z = new double[r, r];

            // a small epsilon is added to the
            //  denominator to avoid overflow.
            double eps = 10e-9;
            int i, j, l, t;
            double s, d;

            for (t = 0; t < iterations; t++)
            {
                var newW = new double[m, r];
                var newH = new double[r, n];

                // Update H using the multiplicative
                // H = H .* (W'*A) ./ (W'*W*H + eps) 
                for (i = 0; i < r; i++)
                {
                    for (j = i; j < r; j++)
                    {
                        s = 0.0;
                        for (l = 0; l < m; l++)
                            s += w[l, i] * w[l, j];
                        Z[i, j] = Z[j, i] = s;
                    }

                    for (j = 0; j < n; j++)
                    {
                        d = 0.0;
                        for (l = 0; l < r; l++)
                            d += Z[i, l] * h[l, j];

                        s = 0.0;
                        for (l = 0; l < m; l++)
                            s += w[l, i] * A[l, j];

                        newH[i, j] = h[i, j] * s / (d + eps);
                    }
                }

                // Update W using the multiplicative
                //   W = W .* (A*H') ./ (W*H*H' + eps)
                for (j = 0; j < r; j++)
                {
                    for (i = j; i < r; i++)
                    {
                        s = 0.0;
                        for (l = 0; l < m; l++)
                            s += newH[i, l] * newH[j, l];
                        Z[i, j] = Z[j, i] = s;
                    }

                    for (i = 0; i < m; i++)
                    {
                        d = 0.0;
                        for (l = 0; l < r; l++)
                            d += w[i, l] * Z[j, l];

                        s = 0.0;
                        for (l = 0; l < n; l++)
                            s += A[i, l] * newH[j, l];

                        newW[i, j] = w[i, j] * s / (d + eps);
                    }
                }

                w = newW;
                h = newH;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет QR-разложение матрицы.
    /// <remarks>
    /// Это представление матрицы в виде произведения двух матриц: A = Q * R, где Q - унитарная (или ортогональная) матрица, а R - верхняя треугольная матрица. 
    /// QR-разложение является основой одного из методов поиска собственных векторов и чисел матрицы — QR-алгоритма.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    public class QR
    {
        #region Private data
        private int m, n;
        private double[][] qr;
        private double[] diag;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует QR-разложение матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        public QR(double[,] A)
        {
            qrdecomp(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает матрицу, содержащую векторы отражений Хаусхолдера.
        /// </summary>
        public double[,] H
        {
            get
            {
                double[,] H = new double[m, n];
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (i >= j)
                        {
                            H[i, j] = qr[i][j];
                        }
                        else
                        {
                            H[i, j] = 0.0f;
                        }
                    }
                }
                return H;
            }

        }
        /// <summary>
        /// Возвращает верхнюю треугольную матрицу R.
        /// </summary>
        public double[,] R
        {
            get
            {
                var r = new double[n, n]; // GeneralMatrix X = new GeneralMatrix(n, n);

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (i < j)
                        {
                            r[i, j] = qr[i][j];
                        }
                        else if (i == j)
                        {
                            r[i, j] = diag[i];
                        }
                        else
                        {
                            r[i, j] = 0.0f;
                        }
                    }
                }
                return r;
            }
        }
        /// <summary>
        /// Возвращает ортогональную матрицу Q.
        /// </summary>
        public double[,] Q
        {
            get
            {
                double[,] q = new double[m, n];
                int i, j, k;
                double s;

                for (k = n - 1; k >= 0; k--)
                {
                    for (i = 0; i < m; i++)
                    {
                        q[i, k] = 0.0f;
                    }
                    q[k, k] = 1.0f;
                    for (j = k; j < n; j++)
                    {
                        if (qr[k][k] != 0)
                        {
                            s = 0.0f;
                            for (i = k; i < m; i++)
                            {
                                s += qr[i][k] * q[i, j];
                            }
                            s = (-s) / qr[k][k];

                            for (i = k; i < m; i++)
                            {
                                q[i, j] += s * qr[i][k];
                            }
                        }
                    }
                }
                return q;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// QR-decomposition.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        private void qrdecomp(double[,] A)
        {
            // params
            this.m = A.GetLength(0);
            this.n = A.GetLength(1);
            this.diag = new double[n];
            this.qr = Jagged.ToJagged(A);
            double nrm, s;
            int k, i, j;

            // Main loop.
            for (k = 0; k < n; k++)
            {
                // Compute 2-norm of k-th column without under/overflow.
                nrm = 0;

                for (i = k; i < m; i++)
                {
                    nrm = Maths.Hypotenuse(nrm, qr[i][k]);
                }

                if (nrm != 0.0)
                {
                    // Form k-th Householder vector.
                    if (qr[k][k] < 0)
                    {
                        nrm = -nrm;
                    }
                    for (i = k; i < m; i++)
                    {
                        qr[i][k] /= nrm; // Make v a unit vector
                    }
                    qr[k][k] += 1.0f; // + the (e)kth vector

                    // Apply transformation to remaining columns.
                    for (j = k + 1; j < n; j++) // For each column
                    {
                        s = 0.0f;
                        for (i = k; i < m; i++) // For each row
                        {
                            s += qr[i][k] * qr[i][j];
                        }

                        s = (-s) / qr[k][k]; // Unit vector product

                        for (i = k; i < m; i++) // For each row
                        {
                            qr[i][j] += s * qr[i][k];
                        }
                    }
                }

                diag[k] = -nrm;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет RQ-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление матрицы в виде произведения двух матриц: A = R * Q, где Q - унитарная (или ортогональная) матрица, а R - верхняя треугольная матрица. 
    /// RQ-разложение является одной из модификаций QR-алгоритма.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    public class RQ
    {
        #region Private data
        private QR qr;
        private double[,] r;
        private double[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует RQ-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public RQ(double[,] A)
        {
            // QR-разложение матрицы A':
            qr = new QR(A.Flip(Direction.Vertical).Transponate());

            // Вычисление матриц R и Q:
            r = qr.R.Transponate();
            q = qr.Q.Transponate();

            // Отображение матриц:
            r = r.Flip(Direction.Both);
            q = q.Flip(Direction.Vertical);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает нижнюю треугольную матрицу R.
        /// </summary>
        public double[,] R
        {
            get
            {
                return this.r;
            }
        }
        /// <summary>
        /// Возвращает ортогональную матрицу Q.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет QL-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление матрицы в виде произведения двух матриц: A = Q * L, где Q - унитарная (или ортогональная) матрица, а L - нижняя треугольная матрица. 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    public class QL
    {
        #region Private data
        private QR qr;
        private double[,] l;
        private double[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует QL-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public QL(double[,] A)
        {
            // QR-разложение матрицы A':
            qr = new QR(A.Flip(Direction.Horizontal));

            // Вычисление матриц L и Q:
            q = qr.Q.Flip(Direction.Horizontal);
            l = qr.R.Flip(Direction.Both);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает нижнюю треугольную матрицу L.
        /// </summary>
        public double[,] L
        {
            get
            {
                return this.l;
            }
        }
        /// <summary>
        /// Возвращает ортогональную матрицу Q.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет LQ-разложение квадратной матрицы.
    /// <remarks>
    /// Это представление матрицы в виде произведения двух матриц: A = L * Q, где Q - унитарная (или ортогональная) матрица, а L - нижняя треугольная матрица. 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/QR_decomposition
    /// </remarks>
    /// </summary>
    public class LQ
    {
        #region Private data
        private QR qr;
        private double[,] l;
        private double[,] q;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует LQ-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public LQ(double[,] A)
        {
            // QR-разложение матрицы A':
            qr = new QR(A.Transponate());

            // Вычисление матриц L и Q:
            l = qr.R.Transponate();
            q = qr.Q.Transponate();
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает нижнюю треугольную матрицу L.
        /// </summary>
        public double[,] L
        {
            get
            {
                return this.l;
            }
        }
        /// <summary>
        /// Возвращает ортогональную матрицу Q.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет процесс ортогонализации Грама-Шмидта.
    /// <remarks>
    /// В математике, в частности линейной алгебре и численном анализе, процесс Грама-Шмидта является методом ортонормирования множества векторов 
    /// в пространстве внутренних произведений. Данная процедура активно используется для ортогонализации базисов.
    /// Более подробную информацию можно найти на сайте: 
    /// https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    /// </remarks>
    /// </summary>
    public class GramSchmidt
    {
        #region Private data
        private double[,] q;
        private double[] v1, v2;
        private double[] u;
        private int n;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует процесс ортогонализации Грама-Шмидта.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public GramSchmidt(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // UMapx.NET
            // gram-schmidt result matrix:
            n = A.GetLength(0);
            q = new double[n, n];
            int i, j;

            for (j = 0; j < n; j++)
            {
                u = Matrice.GetCol(A, j); // get j-column of matrix A,
                v2 = u;                   // copy this column for the second Matrice.

                for (i = 0; i < j; i++)
                {
                    v1 = Matrice.GetCol(q, i); // get i-column of matrix Q
                    u = Matrice.Sub(u, GramSchmidt.Proj(v1, v2)); // calculate: u - proj'<v1, v2>, 
                    // where ' - means transponate operator for projection.
                }

                q = Matrice.SetCol(q, Matrice.Div(u, Matrice.Norm(u)), j); // set j-column of matrix Q.
            }
            return;
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает ортогональную матрицу Q.
        /// </summary>
        public double[,] Q
        {
            get { return q; }
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Возвращает проекцию горизонтальных векторов.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Одномерный массив</param>
        /// <param name="a">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Proj(double[] e, double[] a)
        {
            int length = e.Length;
            double[] proj = new double[length];
            int i;
            double ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * a[i];
                ee += e[i] * e[i];
            }

            double div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }
        /// <summary>
        /// Возвращает проекцию горизонтальных векторов.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Одномерный массив</param>
        /// <param name="a">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Proj(Complex[] e, Complex[] a)
        {
            int length = e.Length;
            Complex[] proj = new Complex[length];
            int i;
            Complex ea = 0, ee = 0;

            for (i = 0; i < length; i++)
            {
                ea += e[i] * (a[i].Conjugate);
                ee += e[i] * (e[i].Conjugate);
            }

            Complex div = ea / ee;

            for (i = 0; i < length; i++)
            {
                proj[i] = e[i] * div;
            }

            return proj;
        }
        #endregion
    }
    /// <summary>
    /// Определяет преобразование Хаусхолдера.
    /// <remarks>
    /// Это линейное преобразование H(u) векторного пространства V, которое описывает его отображение относительно гиперплоскости, 
    /// которая проходит через начало координат. Было предложено в 1958 американским математиком Элстоном Скоттом Хаусхолдером. Широко применяется в линейной алгебре для QR разложения матрицы.
    /// Кроме того, преобразование Хаусхолдера активно используется для ортогонализации базисов, в конечном счете матрица Хаусхолдера обладает свойствами: 
    /// H = H', H' * H = I; det(H) = -1.
    /// В данном классе реализовано два вида преобразования Хаусхолдера: редукция до трехдиагональной матрицы и построение матрицы Хаусхолдера по заданному вектору.
    /// В первом случае исходная квадратная матрица определяется как: A = H * T * H'.
    /// Более подробную информацию можно найти на сайте: 
    /// https://en.wikipedia.org/wiki/Householder_transformation
    /// </remarks>
    /// </summary>
    public class Householder
    {
        #region Private data
        private int n;
        private double[] Re, Im;
        private double[][] matrices;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует преобразование Хаусхолдера.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        public Householder(double[] v)
        {
            // properties:
            this.n = v.Length;
            this.Re = Matrice.One(n);
            this.Im = new double[n];

            // reflection to 
            // Householder matrix:
            hmatx(v);
        }
        /// <summary>
        /// Инициализирует преобразование Хаусхолдера.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public Householder(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // properties:
            this.n = A.GetLength(0);
            this.Re = new double[n];
            this.Im = new double[n];

            // reduction to 
            // tridiagonalization matrix:
            tred2(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает матрицу Хаусхолдера.
        /// </summary>
        public double[,] H
        {
            get
            {
                return Jagged.FromJagged(matrices);
            }
        }
        /// <summary>
        /// Получает диагональную матрицу.
        /// </summary>
        public double[,] T
        {
            get
            {
                double[,] D = new double[n, n];
                int i;

                // diagonal:
                for (i = 0; i < n; i++)
                {
                    D[i, i] = Re[i];
                }
                // diagonal left and right 
                // sides:
                for (i = 1; i < n; i++)
                {
                    D[i - 1, i] = Im[i];
                    D[i, i - 1] = Im[i];
                }

                return D;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует генерацию матрицы Хаусхолдера.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        private void hmatx(double[] v)
        {
            // [1] Alston S. Householder, "Unitary Triangularization of a Nonsymmetric Matrix", 
            // Journal of the ACM 5, 339-242, 1958;
            // [2] G. W. Stewart, Matrix Algorithms: Volume 1: Basic Decompositions, SIAM, xix+458, 1998.
            // 
            // Get Householder vector:
            double[] w = Matrice.Householder(v);

            // Get Householder matrix:
            int n = w.Length, i, j;
            double[] z;
            this.matrices = new double[n][];
            
            // M = I - w * w':
            for (i = 0; i < n; i++)
            {
                // eye vector
                z = new double[n];
                z[i] = 1.0;

                for (j = 0; j < n; j++)
                {
                    z[j] -= w[i] * w[j];
                }

                matrices[i] = z;
            }

            return;
        }
        /// <summary>
        /// Реализует сокращение Хаусхолдера до трехдиагональной формы.
        /// </summary>
        /// <param name="a">Матрица</param>
        private void tred2(double[,] a)
        {
            int i, j, k;
            this.matrices = Jagged.ToJagged(a);

            // Symmetric Householder reduction to tridiagonal form.
            // This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
            }

            double scale, h, f, g, hh;

            // Householder reduction to tridiagonal form.
            for (i = n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                scale = 0;
                h = 0;
                for (k = 0; k < i; k++)
                    scale = scale + System.Math.Abs(Re[k]);

                if (scale == 0)
                {
                    Im[i] = Re[i - 1];
                    for (j = 0; j < i; j++)
                    {
                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                        matrices[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder Matrice.
                    for (k = 0; k < i; k++)
                    {
                        Re[k] /= scale;
                        h += Re[k] * Re[k];
                    }

                    f = Re[i - 1];
                    g = (Double)System.Math.Sqrt(h);
                    if (f > 0) g = -g;

                    Im[i] = scale * g;
                    h = h - f * g;
                    Re[i - 1] = f - g;
                    for (j = 0; j < i; j++)
                        Im[j] = 0;

                    // Apply similarity transformation to remaining columns.
                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        matrices[j][i] = f;
                        g = Im[j] + matrices[j][j] * f;
                        for (k = j + 1; k <= i - 1; k++)
                        {
                            g += matrices[k][j] * Re[k];
                            Im[k] += matrices[k][j] * f;
                        }
                        Im[j] = g;
                    }

                    f = 0;
                    for (j = 0; j < i; j++)
                    {
                        Im[j] /= h;
                        f += Im[j] * Re[j];
                    }

                    hh = f / (h + h);
                    for (j = 0; j < i; j++)
                        Im[j] -= hh * Re[j];

                    for (j = 0; j < i; j++)
                    {
                        f = Re[j];
                        g = Im[j];
                        for (k = j; k <= i - 1; k++)
                            matrices[k][j] -= (f * Im[k] + g * Re[k]);

                        Re[j] = matrices[i - 1][j];
                        matrices[i][j] = 0;
                    }
                }
                Re[i] = h;
            }

            // Accumulate transformations.
            for (i = 0; i < n - 1; i++)
            {
                matrices[n - 1][i] = matrices[i][i];
                matrices[i][i] = 1;
                h = Re[i + 1];
                if (h != 0)
                {
                    for (k = 0; k <= i; k++)
                        Re[k] = matrices[k][i + 1] / h;

                    for (j = 0; j <= i; j++)
                    {
                        g = 0;
                        for (k = 0; k <= i; k++)
                            g += matrices[k][i + 1] * matrices[k][j];
                        for (k = 0; k <= i; k++)
                            matrices[k][j] -= g * Re[k];
                    }
                }

                for (k = 0; k <= i; k++)
                    matrices[k][i + 1] = 0;
            }

            for (j = 0; j < n; j++)
            {
                Re[j] = matrices[n - 1][j];
                matrices[n - 1][j] = 0;
            }

            matrices[n - 1][n - 1] = 1;
            Im[0] = 0;
            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет итерационный алгоритм вычисления собственных значений матрицы.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Power_iteration
    /// </remarks>
    /// </summary>
    public class Power
    {
        #region Private data
        private double[] v;
        #endregion

        #region Power iteration components
        /// <summary>
        /// Инициалазирует итерационный алгоритм вычисления собственных значений матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="iterations">Количество итераций</param>
        public Power(double[,] A, int iterations = 10)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // eigenvalue power algorithm:
            int n = A.GetLength(0);
            this.v = Matrice.Rand(n);
            double[] w;
            double beta;

            // power iteration:
            for (int i = 0; i < iterations; i++)
            {
                // formula:
                // v[j] = (v[j-1] * A) / || v[j-1] * A ||
                w = Matrice.Dot(v, A);
                beta = Matrice.Norm(w);
                v = Matrice.Div(w, beta);
            }
            return;
        }
        /// <summary>
        /// Возвращает вектор собственных значений.
        /// </summary>
        public double[] V
        {
            get
            {
                return v;
            }
        }
        /// <summary>
        /// Возвращает диагонализированную матрицу собственных значений.
        /// </summary>
        public double[,] J
        {
            get
            {
                return Matrice.Diag(v);
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет преобразование Арнольди.
    /// <remarks>
    /// Данное преобразование испольузется для приведения квадратной матрицы к форме Хессенберга.
    /// Матрица A представляется в виде произведения трех матриц: A = Q * H * Q', где H - верхняя треугольная матрица Хессенберга, Q - ортогональная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Arnoldi_iteration
    /// </remarks>
    /// </summary>
    public class Arnoldi
    {
        #region Private data
        private double[,] q;
        private double[,] h;
        #endregion

        #region Arnoldi components
        /// <summary>
        /// Инициализирует преобразование Арнольди.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public Arnoldi(double[,] A)
        {
            // matrix properties
            int n = A.GetLength(0);
            int m = A.GetLength(1);

            if (n != m)
                throw new Exception("Матрица должна быть квадратной");

            // arnoldi decomposition
            arnoldi(A, n, m);
        }
        /// <summary>
        /// Возвращает ортогональную матрицу.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        /// <summary>
        /// Возвращает верхнюю треугольную матрицу Хессенберга.
        /// </summary>
        public double[,] H
        {
            get
            {
                return this.h;
            }
        }
        #endregion

        #region Private data
        /// <summary>
        /// Arnoldi decomposition.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Height</param>
        /// <param name="m">Width</param>
        private void arnoldi(double[,] a, int n, int m)
        {
            // vectors and matrices:
            this.q = new double[n, m];
            this.h = new double[n, m];
            double[,] p = new double[n, m + 1];
            double[] v, w;
            double alpha = 0, beta = 0;
            int i, j, k;

            // random 0-vector and norm:
            v = Matrice.Rand(n);
            p = Matrice.SetCol(p, v, 0);
            p = Matrice.Div(p, Matrice.Norm(v));

            // Start calculating
            // Arnoldi decomposition:
            for (k = 1; k <= m; k++)
            {
                // previous k-1-vector:
                v = Matrice.Dot(Matrice.GetCol(p, k - 1), a);

                for (j = 0; j < k; j++)
                {
                    // calculating α:
                    w = Matrice.GetCol(p, j);
                    alpha = Matrice.Dot(w, v);
                    h[j, k - 1] = alpha;

                    // finding k-vector:
                    for (i = 0; i < n; i++)
                    {
                        v[i] -= w[i] * alpha;
                    }
                }

                // transform:
                if (k < m)
                {
                    // calculating β:
                    beta = Matrice.Norm(v);
                    p = Matrice.SetCol(p, Matrice.Div(v, beta), k);
                    h[k, k - 1] = beta;
                }
            }

            // result Q-matrix:
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < m; j++)
                {
                    q[i, j] = p[i, j];
                }
            }
            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет преобразование Ланцоша.
    /// <remarks>
    /// Данное преобразование испольузется для представления симметричной матрицы A в виде произведения 
    /// трех матриц: A = Q * T * Q', где T - трехдиагональная матрица, а Q - ортогональная матрица.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Lanczos_algorithm
    /// </remarks>
    /// </summary>
    public class Lanczos
    {
        #region Private data
        private double[,] q;
        private double[,] t;
        #endregion

        #region Lanczos components
        /// <summary>
        /// Инициализирует преобразование Ланцоша.
        /// </summary>
        /// <param name="A">Симметричная матрица</param>
        /// <param name="full">Полная реортогонализация или нет</param>
        public Lanczos(double[,] A, bool full = false)
        {
            // exception
            if (!Matrice.IsSymmetric(A))
                throw new Exception("Матрица должна быть симметричной");

            // lanczos decomposition
            int n = A.GetLength(0);
            lanczos(A, n, full);
        }
        /// <summary>
        /// Возвращает ортогональную матрицу.
        /// </summary>
        public double[,] Q
        {
            get
            {
                return this.q;
            }
        }
        /// <summary>
        /// Возвращает трехдиагональную матрицу.
        /// </summary>
        public double[,] T
        {
            get
            {
                return this.t;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Lanczos decomposition.
        /// </summary>
        /// <param name="a">Symmetric matrix</param>
        /// <param name="n">Matrix length</param>
        /// <param name="full">Full re-orthogonalization</param>
        private void lanczos(double[,] a, int n, bool full)
        {
            // This function uses the Lanczos algorithm with full
            // re-orthogonalization to compute k x k symmetric tridiagonal
            // matrix T that approximates mat up to rank k with respect to
            // transformation Q. That is, A = Q * T * Q'.

            // params
            int i, j, y, k = n - 1;
            double[] v = Matrice.Rand(n), z;
            double[] u = v.Div(Matrice.Norm(v));
            this.q = Matrice.Eye(n, n).SetCol(u, 0);
            this.t = new double[n, n];
            double beta, alpha;

            // do job
            for (i = 0; i < n; i++)
            {
                z = u.Dot(a);
                alpha = u.Dot(z);
                t[i, i] = alpha;

                // full re-orthogonalization (x1)
                for (j = 0; j <= i; j++)
                {
                    v = q.GetCol(j);
                    alpha = v.Dot(z);

                    for (y = 0; y < n; y++)
                        z[y] = z[y] - v[y] * alpha;
                }

                // full re-orthogonalization (x2)
                if (full)
                {
                    for (j = 0; j <= i; j++)
                    {
                        v = q.GetCol(j);
                        alpha = v.Dot(z);

                        for (y = 0; y < n; y++)
                            z[y] = z[y] - v[y] * alpha;
                    }
                }

                // tridiagonal matrix
                if (i < k)
                {
                    beta = Matrice.Norm(z);
                    u = z.Div(beta);
                    q = Matrice.SetCol(q, u, i + 1);
                    t[i, i + 1] = beta;
                    t[i + 1, i] = beta;
                }
            }

            return;
        }
        #endregion
    }
    #endregion
}