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
    /// Это представление квадратной матрицы в виде произведения трех матриц: A = P * H * P', где H - матрица Хессенберга, а P - унитарная матрица.
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
            get { return Matrice.FromJagged(matrices); }
        }
        /// <summary>
        /// Получает форму Хессенберга.
        /// </summary>
        public double[,] H
        {
            get { return Matrice.FromJagged(hessenberg); }
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
            this.matrices = Matrice.ToJagged(new double[n, n]);
            this.hessenberg = Matrice.ToJagged(A);
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
            this.eps = eps;

            if (Matrice.IsSymmetric(A))
            {
                hessenberg = Matrice.ToJagged(new double[n, n]);
                matrices = Matrice.ToJagged(A);

                tred2(); // Tridiagonalize.
                tql2();  // Diagonalize.
            }
            else
            {
                matrices = Matrice.ToJagged(new double[n, n]);
                hessenberg = Matrice.ToJagged(A);
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
            get { return Matrice.FromJagged(matrices); }
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
            get { return Matrice.FromJagged(hessenberg); }
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
    /// Данное разложение используется для выделения диагональных матриц в других разложениях (например, LDU-разложение).
    /// </remarks>
    /// </summary>
    public class Diagonal
    {
        #region Private data
        private double[][] matrix;
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

            int n = A.GetLength(0), i, j;
            this.matrix = Matrice.ToJagged(A);
            this.diag = new double[n];

            // Определение значений диагональной матрицы:
            for (i = 0; i < n; i++)
            {
                diag[i] = matrix[i][i];
            }

            double alpha;

            // Вычисление новой матрицы:
            for (j = 0; j < n; j++)
            {
                // digaonal element:
                alpha = diag[j];

                // dividing:
                if (alpha != 0)
                {
                    for (i = 0; i < n; i++)
                    {
                        matrix[i][j] /= alpha;
                    }
                }
            }
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает квадратную матрицу.
        /// </summary>
        public double[,] B
        {
            get { return Matrice.FromJagged(matrix); }
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
            get { return Matrice.FromJagged(lower); }
        }
        /// <summary>
        /// Получает верхнюю треугольную матрицу.
        /// </summary>
        public double[,] U
        {
            get { return Matrice.FromJagged(upper); }
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
            double alpha = 0, beta;
            this.upper = Matrice.ToJagged(new double[n, n]);
            this.lower = Matrice.ToJagged(new double[n, n]);

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
        /// <param name="incomplete">Неполное разложение или нет</param>
        public Cholesky(double[,] A, bool incomplete = false)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            // Cholesky decomposition:
            if (incomplete)
            {
                this.lower = A.ToJagged();
                ichol(this.lower);
            }
            else
            {
                chol(A.ToJagged());
            }
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает нижнюю треугольную матрицу L.
        /// </summary>
        public double[,] L
        {
            get { return Matrice.FromJagged(lower); }
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
            // full Cholesky decomposition
            int n = a.GetLength(0);
            this.lower = Matrice.ToJagged(new double[n, n]);
            int j, i, k;
            double alpha;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i <= j; i++)
                {
                    if (i == j)
                    {
                        alpha = 0;
                        for (k = 0; k < i; k++)
                        {
                            alpha += lower[i][k] * lower[i][k];
                        }
                        lower[i][i] = Math.Sqrt(a[i][i] - alpha);
                    }
                    else if (i < j)
                    {
                        alpha = 0;
                        for (k = 0; k < i; k++)
                        {
                            alpha += lower[i][k] * lower[j][k];
                        }

                        lower[j][i] = 1.0 / lower[i][i] * (a[j][i] - alpha);
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Реализует неполное разложение Холецкого.
        /// </summary>
        /// <param name="a">Квадратная матрица</param>
        private void ichol(double[][] a)
        {
            // incomplete Cholesky decomposition
            int n = a.GetLength(0);
            int j, i, k;

            // cholesky:
            for (k = 0; k < n; k++)
            {
                a[k][k] = Math.Sqrt(a[k][k]);

                for (i = k + 1; i < n; i++)
                {
                    if (a[i][k] != 0)
                    {
                        a[i][k] /= a[k][k];
                    }
                }

                for (j = k + 1; j < n; j++)
                {
                    for (i = j; i < n; i++)
                    {
                        if (a[i][j] != 0)
                        {
                            a[i][j] -= a[i][k] * a[j][k];
                        }
                    }
                }
            }

            // apply zeros:
            for (i = 0; i < n; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    a[i][j] = 0;
                }
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
            double[][] p = Matrice.ToJagged(a);
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
                return reversed ? Matrice.FromJagged(Vr) : Matrice.FromJagged(Ur); 
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
                return reversed ? Matrice.FromJagged(Ur) : Matrice.FromJagged(Vr);
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
                return V.Mul(Matrice.One(m).Div(S)).Dot(U.Transponate());
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
            this.Ur = Matrice.ToJagged(A);
            this.Sr = new double[m];
            this.Vr = Matrice.ToJagged(new double[m, m]);
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
            u = U.Dot(H); p = V.Mul(S).Dot(H);
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
    public class NNMF
    {
        #region Private data
        private double[,] Wr;  // W is m x r (weights),
        private double[,] Hr;  // H is r x n (transformed data) (transposed).
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует неотрицательную матричную факторизацию.
        /// </summary>
        /// <param name="A">Неотрицательная матрица</param>
        /// <param name="r">Размерность новых матриц</param>
        /// <param name="iterations">Количество итераций</param>
        public NNMF(double[,] A, int r, int iterations = 100)
        {
            if (!Matrice.IsPositive(A))
                throw new Exception("Матрица должна быть неотрицательной");

            // properties:
            int n = A.GetLength(0);
            int m = A.GetLength(1);

            if (n < m)
                throw new Exception("Высота матрицы должна быть больше ширины");
            
            // none-negative matrix factorization: 
            nnmf(Matrice.ToJagged(A), n, m, r, iterations);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Получает левую матрицу.
        /// </summary>
        public double[,] W
        {
            get { return Wr; }
        }
        /// <summary>
        /// Получает правую матрицу.
        /// </summary>
        public double[,] H
        {
            get { return Hr; }
        }
        #endregion

        #region Private data
        /// <summary>
        /// Представляет NNMF, основанный на мультипликативном методе.
        /// </summary>
        /// <param name="X">Матрица</param>
        /// <param name="m">Высота матрицы</param>
        /// <param name="n">Ширина матрицы</param>
        /// <param name="r">Новая размерность</param>
        /// <param name="iterations">Количество итераций</param>
        private void nnmf(double[][] X, int n, int m, int r, int iterations)
        {
            // chose W and H randomly, W with unit norm:
            Wr = Matrice.Rand(m, r);
            Hr = Matrice.Rand(r, n);
            double[,] newW, newH, Z = new double[r, r];
            double eps = 10e-9, s, d;
            int i, j, l;

            // Iterative algorithm:
            for (int t = 0; t < iterations; t++)
            {
                newW = new double[m, r];
                newH = new double[r, n];

                // Update H using the multiplicative
                // H = H .* (W' * A) ./ (W' * W * H + eps) 
                for (i = 0; i < r; i++)
                {
                    for (j = i; j < r; j++)
                    {
                        s = 0.0;
                        for (l = 0; l < m; l++)
                            s += Wr[l, i] * Wr[l, j];
                        Z[i, j] = Z[j, i] = s;
                    }

                    for (j = 0; j < n; j++)
                    {
                        d = 0.0;
                        for (l = 0; l < r; l++)
                            d += Z[i, l] * Hr[l, j];

                        s = 0.0;
                        for (l = 0; l < m; l++)
                            s += Wr[l, i] * X[j][l];

                        newH[i, j] = Hr[i, j] * s / (d + eps);
                    }
                }

                // Update W using the multiplicative
                //   W = W .* (A * H') ./ (W * H * H' + eps)
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
                            d += Wr[i, l] * Z[j, l];

                        s = 0.0;
                        for (l = 0; l < n; l++)
                            s += X[l][i] * newH[j, l];

                        newW[i, j] = Wr[i, j] * s / (d + eps);
                    }
                }

                Wr = newW;
                Hr = newH;
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет QR-разложение квадратной матрицы.
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
        private double[][] qr;
        private double[] diag;
        private double[,] q;
        private double[,] r;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует QR-разложение квадратной матрицы.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        public QR(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("Матрица должна быть квадратной");

            qrdecomp(A);
        }
        #endregion

        #region Standart voids
        /// <summary>
        /// Возвращает верхнюю треугольную матрицу R.
        /// </summary>
        public double[,] R
        {
            get
            {
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
                return q;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Реализует QR-разложение.
        /// </summary>
        /// <param name="A">Квадратная матрица</param>
        private void qrdecomp(double[,] A)
        {
            int n = A.GetLength(0);
            this.qr = Matrice.ToJagged(A);
            this.diag = new double[n];
            this.q = new double[n, n];
            this.r = new double[n, n];
            double norm, s;
            int k, i, j;

            // QR-decomposition with diagonal matrix and qr-matrix
            for (k = 0; k < n; k++)
            {
                // Compute 2-norm of k-th column without under/overflow.
                norm = 0;
                for (i = k; i < n; i++)
                {
                    norm = Maths.Hypotenuse(norm, qr[i][k]);
                }

                if (norm != 0)
                {
                    // Form k-th Householder Matrice.
                    if (qr[k][k] < 0)
                    {
                        norm = -norm;
                    }
                    for (i = k; i < n; i++)
                    {
                        qr[i][k] /= norm;
                    }

                    qr[k][k] += 1;

                    // Apply transformation to remaining columns.
                    for (j = k + 1; j < n; j++)
                    {
                        s = 0;
                        for (i = k; i < n; i++)
                        {
                            s += qr[i][k] * qr[i][j];
                        }

                        s = -s / qr[k][k];
                        for (i = k; i < n; i++)
                        {
                            qr[i][j] += s * qr[i][k];
                        }
                    }
                }

                this.diag[k] = -norm;
            }

            // Finding Q-matrix:
            for (k = n - 1; k >= 0; k--)
            {
                for (i = 0; i < n; i++)
                {
                    q[i, k] = 0;
                }

                q[k, k] = 1;
                for (j = k; j < n; j++)
                {
                    if (qr[k][k] != 0)
                    {
                        s = 0;
                        for (i = k; i < n; i++)
                        {
                            s += qr[i][k] * q[i, j];
                        }

                        s = -s / qr[k][k];
                        for (i = k; i < n; i++)
                        {
                            q[i, j] += s * qr[i][k];
                        }
                    }
                }
            }

            // Finding R-matrix:
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (i < j)
                    {
                        r[i, j] = qr[i][j];
                    }
                    else if (i == j)
                    {
                        r[i, j] = diag[i];
                    }
                }
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
                return Matrice.FromJagged(matrices);
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
            double[,] eye = Matrice.Eye(n, n);
            this.matrices = Matrice.ToJagged(new double[n, n]);

            // M = I - w * w':
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    matrices[i][j] = eye[i, j] - w[i] * w[j];
                }
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
            this.matrices = Matrice.ToJagged(a);

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
            // matrix properties:
            int n = A.GetLength(0);
            int m = A.GetLength(1);

            if (n != m)
                throw new Exception("Матрица должна быть квадратной");

            // arnoldi decomposition:
            arnoldi(A, n, m);
        }
        /// <summary>
        /// Возвращает ортогональную матрица.
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
        /// Разложение Арнольди.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="n">Высота</param>
        /// <param name="m">Ширина</param>
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
    #endregion
}