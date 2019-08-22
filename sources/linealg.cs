// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Imaging;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Static components
    /// <summary>
    /// Используется для реализации стандартных алгебраических операций над матрицами и векторами.
    /// </summary>
    public static class Matrice
    {
        // Matrix voids

        #region Matrix booleans
        /// <summary>
        /// Проверяет равенство двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет равенство двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица вектором.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
        public static bool IsVector(this double[,] m)
        {
            if (m.GetLength(0) == 1 || m.GetLength(1) == 1)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли матрица квадратной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
        public static bool IsSquare(this double[,] m)
        {
            if (m.GetLength(0) == m.GetLength(1))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли матрица положительной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица симметричной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица кососимметричной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица диагональной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица вектором.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
        public static bool IsVector(this Complex[,] m)
        {
            if (m.GetLength(0) == 1 || m.GetLength(1) == 1)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли матрица квадратной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
        public static bool IsSquare(this Complex[,] m)
        {
            if (m.GetLength(0) == m.GetLength(1))
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Проверяет является ли матрица симметричной (эрмитовой).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица кососимметричной (антиэрмитовой).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет является ли матрица диагональной.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Логическое значение</returns>
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
        /// Реализует операцию инвертирования матрицы.
        /// </summary>
        /// <param name="m">Квадратная матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Invert(this double[,] m)
        {
            if (!Matrice.IsSquare(m))
            {
                throw new Exception("Матрица должна быть квадратной");
            }

            // Построение матрицы дополнения:
            int n = m.GetLength(0);
            int n2 = n * 2;
            double[][] a = Matrice.ToJagged(new double[n, n2]);
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i][j] = m[i, j];
                }
                a[i][i + n] = 1;
            }

            // Вычисление обратной матрицы через матрицу дополнения:
            const double epsilon = 1e-4; // вычислительная погрешность.
            int k, c, l, t;
            double temp, factor, div1, div2;

            for (i = 0; i < n; i++)
            {
                // first decomposition:
                for (k = i + 1; k < n; k++)
                {
                    if (Math.Abs(a[k][i]) > epsilon)
                    {
                        for (c = 0; c < n2; c++)
                        {
                            temp = a[i][c];
                            a[i][c] = a[k][c];
                            a[k][c] = temp;
                        }
                        break;
                    }
                }
                {
                    // second decomposition:
                    div1 = a[i][i];

                    for (j = 0; j < n2; j++)
                    {
                        if (j != i)
                        {
                            a[i][j] /= div1;
                        }
                    }
                    a[i][i] = 1;
                    div2 = a[i][i];

                    for (t = 0; t < n; t++)
                    {
                        if (t != i)
                        {
                            factor = a[t][i] / div2;

                            for (l = 0; l < n2; l++)
                            {
                                a[t][l] -= factor * a[i][l];
                            }
                        }
                    }
                }
            }

            // building invert matrix:
            double[,] inv = new double[n, n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    inv[i, j] = a[i][j + n];
                }
            }

            return inv;
        }
        /// <summary>
        /// Реализует операцию траспонирования матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Transponate(this double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r1, r0]; // Транспонированная матрица
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
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Зубчатый массив</returns>
        public static double[][] ToJagged(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[][] jagged = new double[ml][];
            double[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new double[mr];
                for (j = 0; j < mr; j++)
                {
                    data[j] = m[i, j];
                }
                jagged[i] = data;
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static double[,] FromJagged(this double[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            double[,] m = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    m[i, j] = jagged[i][j];
                }
            }
            return m;
        }
        /// <summary>
        /// Реализует операцию инвертирования матрицы.
        /// </summary>
        /// <param name="m">Квадратная матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Invert(this Complex[,] m)
        {
            if (!Matrice.IsSquare(m))
            {
                throw new Exception("Матрица должна быть квадратной");
            }

            // Построение матрицы дополнения:
            int n = m.GetLength(0);
            int n2 = n * 2;
            Complex[][] a = Matrice.ToJagged(new Complex[n, n2]);
            int i, j;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a[i][j] = m[i, j];
                }
                a[i][i + n] = 1;
            }

            // Вычисление обратной матрицы через матрицу дополнения:
            const double epsilon = 1e-4; // вычислительная погрешность.
            int k, c, l, t;
            Complex temp, factor, div1, div2;

            for (i = 0; i < n; i++)
            {
                // first decomposition:
                for (k = i + 1; k < n; k++)
                {
                    if (Maths.Abs(a[k][i]) > epsilon)
                    {
                        for (c = 0; c < n2; c++)
                        {
                            temp = a[i][c];
                            a[i][c] = a[k][c];
                            a[k][c] = temp;
                        }
                        break;
                    }
                }
                {
                    // second decomposition:
                    div1 = a[i][i];

                    for (j = 0; j < n2; j++)
                    {
                        if (j != i)
                        {
                            a[i][j] /= div1;
                        }
                    }
                    a[i][i] = 1;
                    div2 = a[i][i];

                    for (t = 0; t < n; t++)
                    {
                        if (t != i)
                        {
                            factor = a[t][i] / div2;

                            for (l = 0; l < n2; l++)
                            {
                                a[t][l] -= factor * a[i][l];
                            }
                        }
                    }
                }
            }

            // building invert matrix:
            Complex[,] inv = new Complex[n, n];

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    inv[i, j] = a[i][j + n];
                }
            }

            return inv;
        }
        /// <summary>
        /// Реализует операцию траспонирования матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Transponate(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r1, r0]; // Транспонированная матрица
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
        /// Возвращает комплексно-сопряженную матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Conjugate(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r0, r1]; // Транспонированная матрица
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
        /// Реализует операцию эрмитово-сопряжения матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Hermitian(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] H = new Complex[r1, r0]; // Транспонированная матрица
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
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex[][] ToJagged(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[][] jagged = new Complex[ml][];
            Complex[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new Complex[mr];
                for (j = 0; j < mr; j++)
                {
                    data[j] = m[i, j];
                }
                jagged[i] = data;
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static Complex[,] FromJagged(this Complex[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            Complex[,] m = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    m[i, j] = jagged[i][j];
                }
            }
            return m;
        }
        #endregion

        #region Matrix properties
        /// <summary>
        /// Возвращает значение ранга матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static int Rank(this double[,] m)
        {
            int height = m.GetLength(0), width = m.GetLength(1);
            int rank = width, row, col, i;
            double mult;

            for (row = 0; row < rank; row++)
            {
                // Before we visit current row 'row', we make
                // sure that mat[row][0],....mat[row][row-1]
                // are 0.

                // Diagonal element is not zero
                if (m[row, row] != 0)
                {
                    for (col = 0; col < height; col++)
                    {
                        if (col != row)
                        {
                            // This makes all entries of current
                            // column as 0 except entry 'mat[row][row]'
                            mult = (double)m[col, row] / m[row, row];
                            for (i = 0; i < rank; i++)
                            {
                                m[col, i] -= mult * m[row, i];
                            }
                        }
                    }
                }

                // Diagonal element is already zero. Two cases
                // arise:
                // 1) If there is a row below it with non-zero
                //    entry, then swap this row with that row
                //    and process that row
                // 2) If all elements in current column below
                //    mat[r][row] are 0, then remvoe this column
                //    by swapping it with last column and
                //    reducing number of columns by 1.
                else
                {
                    bool reduce = true;

                    /* Find the non-zero element in current
                        column  */
                    for (i = row + 1; i < height; i++)
                    {
                        // Swap the row with non-zero element
                        // with this row.
                        if (m[i, row] != 0)
                        {
                            LinealgOptions.Swap(m, row, i, rank);
                            reduce = false;
                            break;
                        }
                    }

                    // If we did not find any row with non-zero
                    // element in current columnm, then all
                    // values in this column are 0.
                    if (reduce)
                    {
                        // Reduce number of columns
                        rank--;

                        // Copy the last column here
                        for (i = 0; i < height; i++)
                        {
                            m[i, row] = m[i, rank];
                        }
                    }

                    // Process this row again
                    row--;
                }
            }
            return rank;
        }
        /// <summary>
        /// Возвращает значение следа квадратной матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Число двойночй точности с плавающей запятой</returns>
        public static double Trace(this double[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("Матрица должна быть квадратной");

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
        /// Возвращает значение определителя матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Число двойночй точности с плавающей запятой</returns>
        public static double Det(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new Exception("Матрица должна быть квадратной");

            unsafe
            {
                fixed (double* pm = &m[0, 0])
                    return LinealgOptions.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Возвращает P-норму матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="p">Параметр p</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает норму матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double Norm(this double[,] m)
        {
            return Matrice.Norm(m, 2);
        }
        /// <summary>
        /// Выделяет целую часть матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Fix(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[,] H = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = Maths.Fix(m[i, j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Реализует построение квадратной матрицы перестановки.
        /// </summary>
        /// <param name="m">Квадратная матрица</param>
        /// <returns>Квадратная матрица</returns>
        public static double[,] Permutation(this double[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("Матрица должна быть квадратной");

            int i, j, r, n = m.GetLength(0);
            double[] temp; double diagonal;
            double[][] perm = Matrice.ToJagged(Matrice.Eye(n, n));

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
            return Matrice.FromJagged(perm);
        }
        /// <summary>
        /// Возвращает значение ранга матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static int Rank(this Complex[,] m)
        {
            int height = m.GetLength(0), width = m.GetLength(1);
            int rank = width, row, col, i;
            Complex mult;

            for (row = 0; row < rank; row++)
            {
                // Before we visit current row 'row', we make
                // sure that mat[row][0],....mat[row][row-1]
                // are 0.

                // Diagonal element is not zero
                if (m[row, row] != 0)
                {
                    for (col = 0; col < height; col++)
                    {
                        if (col != row)
                        {
                            // This makes all entries of current
                            // column as 0 except entry 'mat[row][row]'
                            mult = (Complex)m[col, row] / m[row, row];
                            for (i = 0; i < rank; i++)
                            {
                                m[col, i] -= mult * m[row, i];
                            }
                        }
                    }
                }

                // Diagonal element is already zero. Two cases
                // arise:
                // 1) If there is a row below it with non-zero
                //    entry, then swap this row with that row
                //    and process that row
                // 2) If all elements in current column below
                //    mat[r][row] are 0, then remvoe this column
                //    by swapping it with last column and
                //    reducing number of columns by 1.
                else
                {
                    bool reduce = true;

                    /* Find the non-zero element in current
                        column  */
                    for (i = row + 1; i < height; i++)
                    {
                        // Swap the row with non-zero element
                        // with this row.
                        if (m[i, row] != 0)
                        {
                            LinealgOptions.Swap(m, row, i, rank);
                            reduce = false;
                            break;
                        }
                    }

                    // If we did not find any row with non-zero
                    // element in current columnm, then all
                    // values in this column are 0.
                    if (reduce)
                    {
                        // Reduce number of columns
                        rank--;

                        // Copy the last column here
                        for (i = 0; i < height; i++)
                        {
                            m[i, row] = m[i, rank];
                        }
                    }

                    // Process this row again
                    row--;
                }
            }
            return rank;
        }
        /// <summary>
        /// Возвращает значение следа квадратной матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Комплексное число</returns>
        public static Complex Trace(this Complex[,] m)
        {
            if (!Matrice.IsSquare(m))
                throw new Exception("Матрица должна быть квадратной");

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
        /// Возвращает значение определителя матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Число двойночй точности с плавающей запятой</returns>
        public static Complex Det(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);

            if (mr != ml)
                throw new Exception("Матрица должна быть квадратной");

            unsafe
            {
                fixed (Complex* pm = &m[0, 0])
                    return LinealgOptions.Determinant(pm, mr);
            }
        }
        /// <summary>
        /// Возвращает P-норму двумерной матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="p">Параметр p</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает норму двумерной матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double Norm(this Complex[,] m)
        {
            return Matrice.Norm(m, 2);
        }
        /// <summary>
        /// Выделяет целую часть двумерной матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Fix(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[,] H = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    H[i, j] = Maths.Fix(m[i, j]);
                }
            }

            return H;
        }
        #endregion

        #region Matrix kronecker product
        /// <summary>
        /// Возвращает матричное произведение Кронекера.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует матричное произведение Кронекера.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует матричное произведение Кронекера.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует матричное произведение Кронекера.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Add(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает сумму двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Add(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает сумму двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Add(this Complex[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает сумму двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Add(this double[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает сумму матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает сумму числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Sub(this double[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает разность двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Sub(this Complex[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает разность двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Sub(this Complex[,] m, double[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает разность двух матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Sub(this double[,] m, Complex[,] n)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            int nl = n.GetLength(0), nr = n.GetLength(1);

            if (ml != nl || mr != nr)
                throw new Exception("Матрицы должны быть одинаковых размеров");

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
        /// Возвращает разность матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность матрицы и числа.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность числа и матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Реализует умножение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует умножение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует умножение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Реализует умножение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Умножает все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит матрицу на матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Делит матрицу на матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Делит матрицу на матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Делит матрицу на матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Делит все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит все элементы матрицы на число.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит число на элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит число на элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит число на элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Делит число на элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возводит все элементы матрицы в степень.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="pow">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возводит все элементы матрицы в степень.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="pow">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возводит все элементы матрицы в степень.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="pow">Число</param>
        /// <returns>Матрица</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Логарифмирует все элементы матрицы по основанию.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Экспонирует все элементы матрицы по основанию.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Логарифмирует все элементы матрицы по основанию.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="a">Число</param>
        /// <returns>Матрица</returns>
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
        /// Экспонирует все элементы матрицы по основанию.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Инвертирует все элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Инвертирует все элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает комплексную матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой принадлежат интервалу [0, 255].
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой принадлежат интервалу [0, 1].
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Берет модуль для всех элементов матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Берет модуль для всех элементов матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Берет угол для всех элементов матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
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
        /// Берет действительную часть для всех элементов матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Real(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Re;
                }
            }

            return H;
        }
        /// <summary>
        /// Берет мнимую часть для всех элементов матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Imaginary(this Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] H = new double[r0, r1];
            int i, j;

            for (i = 0; i < r0; i++)
            {
                for (j = 0; j < r1; j++)
                {
                    H[i, j] = m[i, j].Im;
                }
            }

            return H;
        }
        #endregion

        #region Matrix statistics
        /// <summary>
        /// Сортирует матрицу.
        /// </summary>
        /// <param name="m">Матрица</param>
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
        /// Возвращает вектор сумм матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Sum(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] kernel = new double[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] += m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает вектор сумм матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Sum(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[] kernel = new Complex[mr];
            int i, j;

            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] += m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает вектор произведений матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] *= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает вектор произведений матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] *= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает вектор частных матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] /= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает вектор частных матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
                    if (Maths.Singular(m[j, i])) continue;
                    kernel[i] /= m[j, i];
                }
            }
            return kernel;
        }
        /// <summary>
        /// Возвращает максимальный вектор матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Max(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[mr];
            int i, j;

            // Первичное заполнение вектора
            // максимумов:
            for (i = 0; i < mr; i++)
            {
                v[i] = double.MinValue;
            }
            // Поиск максимального вектора:
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
        /// Возвращает минимальный вектор матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Min(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[mr];
            int i, j;

            // Первичное заполнение вектора
            // минимумов:
            for (i = 0; i < mr; i++)
            {
                v[i] = double.MaxValue;
            }
            // Поиск минимального вектора:
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
        /// Возвращает вектор матрицы, соответствующий указанному пороговому значению.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="threshold">Пороговое значение</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Morph(this double[,] m, int threshold)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[] v = new double[ml];
            double[] u = new double[mr];
            int i, j;

            // Поиск минимального вектора:
            for (i = 0; i < mr; i++)
            {
                for (j = 0; j < ml; j++)
                {
                    v[j] = m[j, i];
                }

                // Выделение по порогу:
                Array.Sort(v);
                u[i] = v[threshold];
            }
            return u;
        }
        /// <summary>
        /// Возвращает вектор математического ожидания матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Mean(this double[,] m)
        {
            return Matrice.Div(Matrice.Sum(m), m.GetLength(0));
        }
        /// <summary>
        /// Возвращает вектор математического ожидания матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Mean(this Complex[,] m)
        {
            return Matrice.Div(Matrice.Sum(m), m.GetLength(0));
        }
        /// <summary>
        /// Возвращает вектор дисперсий матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор дисперсий матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор дисперсий матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор дисперсий матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор среднеквадратических отклонений матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] StnDev(this double[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5);
        }
        /// <summary>
        /// Возвращает вектор среднеквадратических отклонений матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] StnDev(this Complex[,] m)
        {
            return Matrice.Pow(Matrice.Var(m), 0.5);
        }
        /// <summary>
        /// Возвращает вектор среднеквадратических отклонений матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] StnDev(this double[,] m, double[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5);
        }
        /// <summary>
        /// Возвращает вектор среднеквадратических отклонений матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] StnDev(this Complex[,] m, Complex[,] n)
        {
            return Matrice.Pow(Matrice.Var(m, n), 0.5);
        }
        /// <summary>
        /// Возвращает матрицу ковариации.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Cov(this double[,] m)
        {
            // Получаем вектора мат. ожиданий μ:
            double[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            double[,] H = (double[,])m.Clone();
            int i, j;

            // Нахождение матрицы E = A - μ:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    H[i, j] -= v[j];
                }
            }
            // Вычисление: cov(A) = 1 / (n  - 1) * E' * E,
            // где ' - знак эрмитова-сопряжения.
            return H.Transponate().Dot(H).Div(height - 1);
        }
        /// <summary>
        /// Возвращает матрицу ковариации.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Cov(this Complex[,] m)
        {
            // Получаем вектора мат. ожиданий μ:
            Complex[] v = Matrice.Mean(m);
            int width = m.GetLength(1), height = m.GetLength(0);
            Complex[,] H = (Complex[,])m.Clone();
            int i, j;

            // Нахождение матрицы E = A - μ:
            for (i = 0; i < height; i++)
            {
                for (j = 0; j < width; j++)
                {
                    H[i, j] -= v[j];
                }
            }
            // Вычисление: cov(A) = 1 / (n  - 1) * E' * E,
            // где ' - знак эрмитова-сопряжения.
            return H.Hermitian().Dot(H).Div(height - 1);
        }
        /// <summary>
        /// Возвращает вектор энтропии матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует скалярное произведение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static double[,] Dot(this double[,] m, double[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Реализует скалярное произведение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Dot(this Complex[,] m, Complex[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Реализует скалярное произведение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Dot(this Complex[,] m, double[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        /// <summary>
        /// Реализует скалярное произведение матриц.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Dot(this double[,] m, Complex[,] n)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Mul(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n)));
        }
        #endregion

        #region Matrix convolutions
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
        public static double[,] Conv(this double[,] m, double[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Conv(this Complex[,] m, Complex[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Conv(this Complex[,] m, double[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Conv(this double[,] m, Complex[,] n, bool normalize = true)
        {
            return LinealgOptions.FromJagged(LinealgOptions.Conv(LinealgOptions.ToJagged(m), LinealgOptions.ToJagged(n), normalize));
        }
        #endregion

        #region Matrix convolutions (separable)
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу-результат морфологического сужения.
        /// </summary>
        /// <param name="m">Двумерный массив</param>
        /// <param name="r0">Радиус обработки по высоте</param>
        /// <param name="r1">Радиус обработки по ширине</param>
        public static double[,] Min(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MinVertical(
                                             LinealgOptions.MinHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологического расширения.
        /// </summary>
        /// <param name="m">Двумерный массив</param>
        /// <param name="r0">Радиус обработки по высоте</param>
        /// <param name="r1">Радиус обработки по ширине</param>
        public static double[,] Max(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MaxVertical(
                                             LinealgOptions.MaxHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологии.
        /// </summary>
        /// <param name="m">Двумерный массив</param>
        /// <param name="r0">Радиус обработки по высоте</param>
        /// <param name="r1">Радиус обработки по ширине</param>
        /// <param name="threshold">Пороговое значение</param>
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
        /// Возвращает матрицу-результат локального усреднения.
        /// </summary>
        /// <param name="m">Двумерный массив</param>
        /// <param name="r0">Радиус обработки по высоте</param>
        /// <param name="r1">Радиус обработки по ширине</param>
        public static double[,] Mean(this double[,] m, int r0, int r1)
        {
            // both processing
            float[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MeanVertical(
                                             LinealgOptions.MeanHorizontal(mm, r1), r0));
        }
        /// <summary>
        /// Возвращает матрицу-результат локального усреднения.
        /// </summary>
        /// <param name="m">Двумерный массив</param>
        /// <param name="r0">Радиус обработки по высоте</param>
        /// <param name="r1">Радиус обработки по ширине</param>
        public static Complex[,] Mean(this Complex[,] m, int r0, int r1)
        {
            // both processing
            LinealgOptions.Complex32[][] mm = LinealgOptions.ToJagged(m);
            return LinealgOptions.FromJagged(LinealgOptions.MeanVertical(
                                             LinealgOptions.MeanHorizontal(mm, r1), r0));
        }
        #endregion

        // Matrix <-> Bitmap voids

        #region Bitmap matrix voids
        /// <summary>
        /// Преобразовывает точечный рисунок в матрицу значений усредненного канала.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Прямоугольная матрица</returns>
        public static double[,] FromBitmap(this Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[,] rgb = FromBitmap(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в матрицу значений усредненного канала.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Прямоугольная матрица</returns>
        public unsafe static double[,] FromBitmap(this BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] rgb = new double[height, width];
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    rgb[j, i] = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0;
                }
            });

            return rgb;
        }
        /// <summary>
        /// Преобразовывает прямоугольную матрицу значений каналов в монохромный точечный рисунок.
        /// </summary>
        /// <param name="m">Прямоугольная матрица</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap ToBitmap(this double[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0);
                    p[k + 3] = 255;
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        #endregion

        // Vector voids

        #region Vector booleans
        /// <summary>
        /// Проверяет являются ли векторы одинаковыми.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет являются ли векторы одинаковыми.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
        public static bool IsEquals(this Complex[] a, double[] b)
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
        /// Проверяет являются ли векторы коллинеарными.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет являются ли векторы коллинеарными.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет являются ли векторы коллинеарными.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет являются ли векторы коллинеарными.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Проверяет являются ли векторы коллинеарными.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Логическое значение</returns>
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
        /// Возвращает P-норму вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="p">Параметр p</param>
        /// <returns>Норма</returns>
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
        /// Возвращает норму вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <returns>Норма</returns>
        public static double Norm(this double[] a)
        {
            return Norm(a, 2);
        }
        /// <summary>
        /// Возвращает P-норму вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="p">Параметр p</param>
        /// <returns>Норма</returns>
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
        /// Возвращает норму вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <returns>Норма</returns>
        public static double Norm(this Complex[] a)
        {
            return Norm(a, 2);
        }
        #endregion

        #region Vector angle, projection, cosine
        /// <summary>
        /// Возвращает угол между двумя векторами.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Angle(this double[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает угол между двумя векторами.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Angle(this Complex[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает угол между двумя векторами.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Angle(this double[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает угол между двумя векторами.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Angle(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(a) / Matrice.Norm(b);
        }

        /// <summary>
        /// Возвращает проекцию двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Proj(this double[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает проекцию двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Proj(this Complex[] a, double[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает проекцию двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Proj(this double[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }
        /// <summary>
        /// Возвращает проекцию двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Proj(this Complex[] a, Complex[] b)
        {
            return Matrice.Dot(a, b) / Matrice.Norm(b);
        }

        /// <summary>
        /// Возвращает направляющие косинусы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает направляющие косинусы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает сумму вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность двух векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность вектора и числа.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность числа и вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность числа и вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность числа и вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность числа и вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[,] Mul(this double[,] m, double[] v)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            double[,] temp = new double[r0, r1];
            double alpha;
            int i, j;

            // Вычисление новой матрицы:
            for (j = 0; j < r1; j++)
            {
                alpha = v[j];
                for (i = 0; i < r0; i++)
                {
                    temp[i, j] = m[i, j] * alpha;
                }
            }

            return temp;
        }
        /// <summary>
        /// Реализует умножение матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Mul(this Complex[,] m, double[] v)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
            int i, j;

            // Вычисление новой матрицы:
            for (j = 0; j < r1; j++)
            {
                alpha = v[j];
                for (i = 0; i < r0; i++)
                {
                    temp[i, j] = m[i, j] * alpha;
                }
            }

            return temp;
        }
        /// <summary>
        /// Реализует умножение матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Mul(this double[,] m, Complex[] v)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
            int i, j;

            // Вычисление новой матрицы:
            for (j = 0; j < r1; j++)
            {
                alpha = v[j];
                for (i = 0; i < r0; i++)
                {
                    temp[i, j] = m[i, j] * alpha;
                }
            }

            return temp;
        }
        /// <summary>
        /// Реализует умножение матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Mul(this Complex[,] m, Complex[] v)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            Complex[,] temp = new Complex[r0, r1];
            Complex alpha;
            int i, j;

            // Вычисление новой матрицы:
            for (j = 0; j < r1; j++)
            {
                alpha = v[j];
                for (i = 0; i < r0; i++)
                {
                    temp[i, j] = m[i, j] * alpha;
                }
            }

            return temp;
        }

        /// <summary>
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует умножение вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[,] Div(this double[,] m, double[] v)
        {
            int c = m.GetLength(1);
            int r = m.GetLength(0);
            double[,] H = new double[r, c];
            double alpha;
            int i, j;

            for (j = 0; j < c; j++)
            {
                // digaonal element:
                alpha = v[j];

                // dividing:
                if (alpha != 0)
                {
                    for (i = 0; i < r; i++)
                    {
                        H[i, j] = m[i, j] / alpha;
                    }
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует деление матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Div(this Complex[,] m, double[] v)
        {
            int c = m.GetLength(1);
            int r = m.GetLength(0);
            Complex[,] H = new Complex[r, c];
            Complex alpha;
            int i, j;

            for (j = 0; j < c; j++)
            {
                // digaonal element:
                alpha = v[j];

                // dividing:
                if (alpha != 0)
                {
                    for (i = 0; i < r; i++)
                    {
                        H[i, j] = m[i, j] / alpha;
                    }
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует деление матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Div(this double[,] m, Complex[] v)
        {
            int c = m.GetLength(1);
            int r = m.GetLength(0);
            Complex[,] H = new Complex[r, c];
            Complex alpha;
            int i, j;

            for (j = 0; j < c; j++)
            {
                // digaonal element:
                alpha = v[j];

                // dividing:
                if (alpha != 0)
                {
                    for (i = 0; i < r; i++)
                    {
                        H[i, j] = m[i, j] / alpha;
                    }
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует деление матрицы на вектор вида: A * diag(v).
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Div(this Complex[,] m, Complex[] v)
        {
            int c = m.GetLength(1);
            int r = m.GetLength(0);
            Complex[,] H = new Complex[r, c];
            Complex alpha;
            int i, j;

            for (j = 0; j < c; j++)
            {
                // digaonal element:
                alpha = v[j];

                // dividing:
                if (alpha != 0)
                {
                    for (i = 0; i < r; i++)
                    {
                        H[i, j] = m[i, j] / alpha;
                    }
                }
            }
            return H;
        }

        /// <summary>
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует поэлементное произведение векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление вектора на число.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление числа на вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление числа на вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление числа на вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует деление числа на вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит элементы вектора в степень.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="power">Степень</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит элементы вектора в степень.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="power">Степень</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит элементы вектора в степень.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="power">Степень</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возводит число поэлементно в степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Логарифмирует элементы вектора по основанию.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Логарифмирует элементы вектора по основанию.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="a">Число</param>
        /// <returns>Одномерный массив</returns>
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
        /// Экспонирует элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Экспонирует элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор, значения которого принадлежат интервалу [0, 1].
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор, значения которого принадлежат интервалу [0, 255].
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает модуль элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Инвертирует все элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Инвертирует все элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает комплексный вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает модуль элементов комплексного вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает угол элементов комплексного вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает действительную часть элементов комплексного вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Real(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Re;
            }
            return H;
        }
        /// <summary>
        /// Возвращает мнимую часть элементов комплексного вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Imaginary(this Complex[] v)
        {
            int length = v.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = v[i].Im;
            }
            return H;
        }
        /// <summary>
        /// Возвращает комплексно-сопряженный вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает значение ковариации вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение ковариации вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
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
        /// Возвращает общее значение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Sum(this double[] v)
        {
            double total = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total += v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает общее значение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
        public static Complex Sum(this Complex[] v)
        {
            Complex total = 0;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total += v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает общее произведение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Mul(this double[] v)
        {
            double total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total *= v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает общее произведение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Mul(this Complex[] v)
        {
            Complex total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total *= v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает общее частное вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Div(this double[] v)
        {
            double total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total /= v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает общее частное вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static Complex Div(this Complex[] v)
        {
            Complex total = 1;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (Maths.Singular(v[i])) continue;
                total /= v[i];
            }
            return total;
        }
        /// <summary>
        /// Возвращает значение энтропии вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает среднее значение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Mean(this double[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Возвращает среднее значение вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
        public static Complex Mean(this Complex[] v)
        {
            return Matrice.Sum(v) / v.Length;
        }
        /// <summary>
        /// Возвращает значение дисперсии.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение дисперсии.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
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
        /// Возвращает значение дисперсии.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает значение дисперсии.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
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
        /// Возвращает значение среднеквадратического отклонения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double StnDev(this double[] v)
        {
            return Math.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Возвращает значение среднеквадратического отклонения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комлпексное число</returns>
        public static Complex StnDev(this Complex[] v)
        {
            return Maths.Sqrt(Matrice.Var(v));
        }
        /// <summary>
        /// Возвращает значение среднеквадратического отклонения.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double StnDev(this double[] x, double[] y)
        {
            return Math.Sqrt(Matrice.Var(x, y));
        }
        /// <summary>
        /// Возвращает значение среднеквадратического отклонения.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Комлпексное число</returns>
        public static Complex StnDev(this Complex[] x, Complex[] y)
        {
            return Maths.Sqrt(Matrice.Var(x, y));
        }
        /// <summary>
        /// Возвращает значение моды вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число одиночной точности с плавающей запятой</returns>
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
        /// Возвращает значение моды вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Комплексное число</returns>
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
        /// Возвращает значения минимума и максимума вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Пара целых чисел, представляющая отрезок</returns>
        public static RangeDouble Extremum(this double[] v)
        {
            double max = int.MinValue, min = int.MaxValue;
            int length = v.Length;

            for (int i = 0; i < length; i++)
            {
                if (v[i] > max)
                    max = v[i];
                if (v[i] < min)
                    min = v[i];
            }
            return new RangeDouble(min, max);
        }
        /// <summary>
        /// Получает значение минимального элемента сигнала.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число одиночной точности с плавающей запятой</returns>
        public static double Min(this double[] v)
        {
            int index = 0, length = v.Length;
            double minimum = v[index];

            for (int i = 1; i < length; i++)
            {
                if (v[i] < minimum)
                {
                    minimum = v[i];
                    index = i;
                }
            }
            return minimum;
        }
        /// <summary>
        /// Получает значение максимального элемента вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Число одиночной точности с плавающей запятой</returns>
        public static double Max(this double[] v)
        {
            int index = 0, length = v.Length;
            double maximum = v[index];

            for (int i = 1; i < length; i++)
            {
                if (v[i] > maximum)
                {
                    maximum = v[i];
                    index = i;
                }
            }
            return maximum;
        }
        /// <summary>
        /// Получает значение элемента вектора, соответствующего пороговому значению.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="threshold">Пороговое значение</param>
        /// <returns>Число одиночной точности с плавающей запятой</returns>
        public static double Morph(this double[] v, int threshold)
        {
            double[] u = (double[])v.Clone();
            Array.Sort(u);
            return u[threshold];
        }
        /// <summary>
        /// Сортирует исходный вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
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
        /// Сортирует исходный вектор.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Определяет диапазон сортировки. Начальная точка</param>
        /// <param name="l">Определяет диапазон сортировки. Конечная точка</param>
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

        // Vector special

        #region Vector dot
        /// <summary>
        /// Реализует скалярное произведение векторов вида: a * b'.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Реализует скалярное произведение векторов вида: a * b'.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Реализует скалярное произведение векторов вида: a * b'.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Реализует скалярное произведение векторов вида: a * b'.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Реализует скалярное умножение вектора на матрицу вида: 
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Dot(this double[] v, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            if (v.Length != r1)
                throw new Exception("Размерность вектора должна быть равна ширине матрицы");

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
        /// Реализует скалярное умножение вектора на матрицу.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Dot(this double[] v, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            if (v.Length != r1)
                throw new Exception("Размерность вектора должна быть равна ширине матрицы");

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
        /// Реализует скалярное умножение вектора на матрицу.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Dot(this Complex[] v, Complex[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            if (v.Length != r1)
                throw new Exception("Размерность вектора должна быть равна ширине матрицы");

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
        /// Реализует скалярное умножение вектора на матрицу.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="m">Матрица</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Dot(this Complex[] v, double[,] m)
        {
            int r0 = m.GetLength(0), r1 = m.GetLength(1);
            if (v.Length != r1)
                throw new Exception("Размерность вектора должна быть равна ширине матрицы");

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
        /// Реализует скалярное произведение векторов вида: a' * b, 
        /// где ' - знак транспонирования.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static double[,] Dotp(this double[] a, double[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            double[,] H = new double[l0, l1];
            int i, j;

            for (j = 0; j < l0; j++)
            {
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += a[j] * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует скалярное произведение векторов вида: a' * b, 
        /// где ' - знак транспонирования.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Dotp(this Complex[] a, Complex[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            int i, j;

            for (j = 0; j < l0; j++)
            {
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += a[j] * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует скалярное произведение векторов вида: a' * b, 
        /// где ' - знак транспонирования.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Dotp(this double[] a, Complex[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            int i, j;

            for (j = 0; j < l0; j++)
            {
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += a[j] * b[i];
                }
            }
            return H;
        }
        /// <summary>
        /// Реализует скалярное произведение векторов вида: a' * b, 
        /// где ' - знак транспонирования.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[,] Dotp(this Complex[] a, double[] b)
        {
            int l0 = a.Length, l1 = b.Length;
            Complex[,] H = new Complex[l0, l1];
            int i, j;

            for (j = 0; j < l0; j++)
            {
                for (i = 0; i < l1; i++)
                {
                    H[j, i] += a[j] * b[i];
                }
            }
            return H;
        }
        #endregion

        #region Vector convolutions
        /// <summary>
        /// Возвращает вектор-результат одномерной свертки.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="u">Одномерный массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор-результат одномерной свертки.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="u">Одномерный массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор-результат одномерной свертки.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="u">Одномерный массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор-результат одномерной свертки.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="u">Одномерный массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует одномерный фильтр морфологического расширения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Радиус обработки</param>
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
        /// Реализует одномерный фильтр морфологического сужения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Радиус обработки</param>
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
        /// Возвращает вектор-результат одномерной морфологии.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Радиус</param>
        /// <param name="threshold">Пороговое значение</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Morph(this double[] v, int r, int threshold)
        {
            int n = v.Length;
            double[] u = new double[r];
            double[] uv = new double[n];
            int r2 = r / 2;
            int i, j, k, c, p;

            // Операция морфологии:
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

                // Сортировка массива:
                Array.Sort(u);
                // Выборка:
                uv[i] = u[threshold];
            }

            return uv;
        }
        #endregion

        #region Vector mean
        /// <summary>
        /// Реализует одномерный фильтр локального усреднения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Радиус</param>
        public static double[] Mean(this double[] v, int r)
        {
            // Определение параметров фильтра:
            int l = v.Length;
            double[] output = new double[l];
            int r2 = r >> 1;
            int dl = l - r2;
            double s = 0;
            int i;

            // Вычисление глобальной суммы [0, h):
            for (i = 0; i < r; i++)
            {
                s += v[i];
            }
            // Вычисление фильтра на отрезке [0, v):
            for (i = 0; i < r2; i++)
            {
                output[i] = s / r;
            }
            // Вычисление фильтра на отрезке [v, l-v):
            for (i = r2; i < dl; i++)
            {
                s = s - v[i - r2] + v[i + r2];
                output[i] = s / r;
            }
            // Вычисление фильтра на отрезке [l-v, l):
            for (i = dl; i < l; i++)
            {
                s = s - v[i - r2] + v[i];
                output[i] = s / r;
            }

            return output;
        }
        /// <summary>
        /// Реализует одномерный фильтр локального усреднения.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="r">Радиус</param>
        public static Complex[] Mean(this Complex[] v, int r)
        {
            // Определение параметров фильтра:
            int l = v.Length;
            Complex[] output = new Complex[l];
            int r2 = r >> 1;
            int dl = l - r2;
            Complex s = 0;
            int i;

            // Вычисление глобальной суммы [0, h):
            for (i = 0; i < r; i++)
            {
                s += v[i];
            }
            // Вычисление фильтра на отрезке [0, v):
            for (i = 0; i < r2; i++)
            {
                output[i] = s / r;
            }
            // Вычисление фильтра на отрезке [v, l-v):
            for (i = r2; i < dl; i++)
            {
                s = s - v[i - r2] + v[i + r2];
                output[i] = s / r;
            }
            // Вычисление фильтра на отрезке [l-v, l):
            for (i = dl; i < l; i++)
            {
                s = s - v[i - r2] + v[i];
                output[i] = s / r;
            }

            return output;
        }
        #endregion

        // MATLAB voids

        #region MATLAB voids
        #region Get/set rows and columns
        /// <summary>
        /// Возвращает вектор столбца матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="r">Номер столбца</param>
        /// <returns>Одномерный массив</returns>
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
        /// Задает вектор столбца матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Одномерный массив</param>
        /// <param name="r">Номер столбца</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает вектор строки матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="r">Номер строки</param>
        /// <returns>Одномерный массив</returns>
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
        /// Задает вектор строки матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Одномерный массив</param>
        /// <param name="r">Номер строки</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает вектор столбца матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="r">Номер столбца</param>
        /// <returns>Одномерный массив</returns>
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
        /// Задает вектор столбца матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Одномерный массив</param>
        /// <param name="r">Номер столбца</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает вектор строки матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="r">Номер строки</param>
        /// <returns>Одномерный массив</returns>
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
        /// Задает вектор строки матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="n">Одномерный массив</param>
        /// <param name="r">Номер строки</param>
        /// <returns>Матрица</returns>
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
        /// Релизует блочную перегруппировку матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="m">Количество позиций, на которое происходит сдвиг по высоте</param>
        /// <param name="l">Количество позиций, на которое происходит сдвиг по ширине</param>
        /// <returns>Матрица</returns>
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
        /// Релизует блочную перегруппировку матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="m">Количество позиций, на которое происходит сдвиг по высоте</param>
        /// <param name="l">Количество позиций, на которое происходит сдвиг по ширине</param>
        /// <returns>Матрица</returns>
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
        /// Релизует сдвиг элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="l">Количество позиций, на которое происходит сдвиг</param>
        /// <returns>Одномерный массив</returns>
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
        /// Релизует сдвиг элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="l">Количество позиций, на которое происходит сдвиг</param>
        /// <returns>Одномерный массив</returns>
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
        /// Отображает элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <returns>Матрица</returns>
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
        /// Отображает элементы матрицы.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="direction">Направление обработки</param>
        /// <returns>Матрица</returns>
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
        /// Отображает элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Отображает элементы вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует слияние векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует слияние векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует слияние векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует слияние векторов.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="b">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает заданную часть вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="start">Начальная позиция</param>
        /// <param name="length">Длина вектора</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Cut(this double[] a, int start, int length)
        {
            int na = a.Length, i;
            double[] v = new double[length];

            for (i = 0; i < length; i++)
                v[i] = a[Maths.Mod(start + i, na)];

            return v;
        }
        /// <summary>
        /// Возвращает заданную часть вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="start">Начальная позиция</param>
        /// <param name="length">Длина вектора</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Cut(this Complex[] a, int start, int length)
        {
            int na = a.Length, i;
            Complex[] v = new Complex[length];

            for (i = 0; i < length; i++)
                v[i] = a[Maths.Mod(start + i, na)];

            return v;
        }
        /// <summary>
        /// Обрезает матрицу до заданных размеров.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="y">Начальная позиция по высоте</param>
        /// <param name="x">Начальная позиция по ширине</param>
        /// <param name="height">Высота</param>
        /// <param name="width">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Обрезает матрицу до заданных размеров.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="y">Начальная позиция по высоте</param>
        /// <param name="x">Начальная позиция по ширине</param>
        /// <param name="height">Высота</param>
        /// <param name="width">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, сформированную из вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="height">Высота матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, сформированную из вектора.
        /// </summary>
        /// <param name="a">Одномерный массив</param>
        /// <param name="height">Высота матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает вектор, сформированный из матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="length">Размерность вектора</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор, сформированный из матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="length">Размерность вектора</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует приведение вектора к диагональной матрице.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует приведение вектора к диагональной матрице.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает вектор, элементы которого лежат на диагонали матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает вектор, элементы которого лежат на диагонали матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="i">Первый ряд или столбец</param>
        /// <param name="j">Второй ряд или столбец</param>
        /// <param name="direction">Направление обработки</param>
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
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="i">Первый ряд или столбец</param>
        /// <param name="j">Второй ряд или столбец</param>
        /// <param name="direction">Направление обработки</param>
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
        /// Реализует перестановку элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="i">Номер первого элемента</param>
        /// <param name="j">Номер второго элемента</param>
        public static void Swap(this double[] v, int i, int j)
        {
            // get elements:
            double e1 = v[i], e2 = v[j];
            // swapping vector elements:
            v[j] = e1; v[i] = e2;
            return;
        }
        /// <summary>
        /// Реализует перестановку элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="i">Номер первого элемента</param>
        /// <param name="j">Номер второго элемента</param>
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
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="i">Первый ряд или столбец</param>
        /// <param name="length">Длина</param>
        /// <param name="direction">Направление обработки</param>
        /// <returns>Матрица</returns>
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
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="i">Первый ряд или столбец</param>
        /// <param name="length">Длина</param>
        /// <param name="direction">Направление обработки</param>
        /// <returns>Матрица</returns>
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
        /// Реализует удаление элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="i">Номер элемента</param>
        /// <param name="length">Длина</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует удаление элементов вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="i">Номер элемента</param>
        /// <param name="length">Длина</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует операцию взятия минора матрицы.
        /// </summary>
        /// <param name="m">Квадратная матрица</param>
        /// <param name="n">Номер строки и столбца</param>
        /// <returns>Квадратная матрица</returns>
        public static double[,] Minor(this double[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new Exception("Матрица должна быть квадратной");
            if (n >= height || n < 0) throw new Exception("Номер строки и столбца указан не верно");

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
        /// Реализует операцию взятия минора матрицы.
        /// </summary>
        /// <param name="m">Квадратная матрица</param>
        /// <param name="n">Номер строки и столбца</param>
        /// <returns>Квадратная матрица</returns>
        public static Complex[,] Minor(this Complex[,] m, int n)
        {
            // matrix sizes:
            int height = m.GetLength(0), width = m.GetLength(1);

            // errors:
            if (height != width) throw new Exception("Матрица должна быть квадратной");
            if (n >= height || n < 0) throw new Exception("Номер строки и столбца указан не верно");

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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="n">Порядок</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="n">Порядок</param>
        /// <param name="direction">Направление обработки</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="n">Порядок</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает разность элементов массива.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="n">Порядок</param>
        /// <param name="reverse">Обратное направление обработки или нет</param>
        /// <returns>Одномерный массив</returns>
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
        /// Расширяет вектор до заданной длины.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="length">Длина</param>
        /// <returns>Одномерный массив</returns>
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
        /// Расширяет матрицу до заданных размеров.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="height">Высота</param>
        /// <param name="width">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Расширяет вектор до заданной длины.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="length">Длина</param>
        /// <returns>Одномерный массив</returns>
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
        /// Расширяет матрицу до заданных размеров.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="height">Высота</param>
        /// <param name="width">Ширина</param>
        /// <returns>Матрица</returns>
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
        #endregion

        // Extra voids

        #region Compute methods
        /// <summary>
        /// Возвращает массив значений аргумента, реализованный в заданном отрезке с заданным шагом.
        /// </summary>
        /// <param name="min">Начало отрезка</param>
        /// <param name="max">Конец отрезка</param>
        /// <param name="step">Шаг</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Compute(double min, double max, double step)
        {
            // ******************************
            // MATLAB vector computing void
            // designed by Asiryan Valeriy
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
        /// Возвращает массив значений функции.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает массив значений функции.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает матрицу, значения которой определяются по двум массивам аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="x">Массив значений первого аргумента</param>
        /// <param name="y">Массив значений второго аргумента</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой определяются по двум массивам аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="x">Массив значений первого аргумента</param>
        /// <param name="y">Массив значений второго аргумента</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой определяются по двум массивам аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="x">Массив значений первого аргумента</param>
        /// <param name="y">Массив значений второго аргумента</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой определяются по двум массивам аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="x">Массив значений первого аргумента</param>
        /// <param name="y">Массив значений второго аргумента</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой определяются по массиву аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает матрицу, значения которой определяются по массиву аргумента и делегату непрерывной функции.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <param name="function">Делегат непрерывной функции</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение вектора единиц.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
        public static double[] One(int n)
        {
            double[] v = new double[n];

            for (int i = 0; i < n; i++)
            {
                v[i] = 1;
            }

            return v;
        }
        #endregion

        #region Matrix products
        /// <summary>
        /// Возвращает вектор Хаусхолдера.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует построение сопровождающей матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Вандерморта.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение неполной матрицы Ганкеля.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Ганкеля.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Тёплица.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Коши.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение цикличной матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение симметричной матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение сопровождающей матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Вандерморта.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение неплоной матрицы Ганкеля.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Ганкеля.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Тёплица.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Коши.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="y">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение цикличной матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение симметричной матрицы.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение единичной матрицы.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы единиц.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение обменной матрицы.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Лемера.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Редхеффера.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
        public static double[,] Redheffer(int n)
        {
            double[,] H = new double[n, n];
            int i, j;

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    // Главная диагональ:
                    if (i == j)
                    {
                        H[i - 1, j - 1] = 1;
                    }
                    //  i делится на j без остатка:
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
        /// Реализует построение матрицы Гильберта.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение цикличной матрицы.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение симметричной матрицы.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы НОД-ов.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы Стирлинга первого или второго рода.
        /// </summary>
        /// <param name="n">Размерность матрицы</param>
        /// <param name="second">Второго рода или нет</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение магического квадрата.
        /// </summary>
        /// <param name="n">Размер матрицы (нечетное число)</param>
        /// <returns>Матрица</returns>
        public static double[,] Magic(int n)
        {
            if (Maths.Mod(n, 2) != 1)
                throw new Exception("Размерность матрица должна быть нечетным числом");

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
        /// Генератор случайных чисел.
        /// </summary>
        private static Random rnd = new Random();

        #region Rand and randc voids
        /// <summary>
        /// Реализует построение вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует построение комплексного вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует построение матрицы случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение матрицы случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Randi(int n)
        {
            return Randi(n, 1, n + 1);
        }
        /// <summary>
        /// Реализует построение комплексного вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Randic(int n)
        {
            return Randic(n, 1, n + 1);
        }
        /// <summary>
        /// Реализует построение вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует построение комплексного вектора случайных дробных чисел, значения которой распределены по равномерному закону.
        /// </summary>
        /// <param name="n">Размерность</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует построение матрицы случайных целых чисел.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
        public static double[,] Randi(int m, int l)
        {
            return Randi(m, l, 1, l + 1);
        }
        /// <summary>
        /// Реализует построение матрицы случайных целых чисел.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Матрица</returns>
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
        /// Реализует комплексной построение матрицы случайных целых чисел.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
        public static Complex[,] Randic(int m, int l)
        {
            return Randic(m, l, 1, l + 1);
        }
        /// <summary>
        /// Реализует комплексной построение матрицы случайных целых чисел.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="a">Нижний предел</param>
        /// <param name="b">Верхний предел</param>
        /// <returns>Матрица</returns>
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
        /// Переводит исходную строку в матрицу вещественных чисел.
        /// <remarks>
        /// Пример входной строки: "[1, 2, 3; 4, 5, 6; 7, 8, 9]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="s">Исходная строка</param>
        /// <returns>Матрица</returns>
        public static double[,] Parse(this double[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length, n = nums.Length;
            double[,] H = new double[r, n];
            int i, j;

            // collecting rows:
            for (i = 0; i < r; i++)
            {
                nums = rows[i].Split('|');
                n = nums.Length;

                for (j = 0; j < n; j++)
                {
                    H[i, j] = double.Parse(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Пробует перевести исходную строку в матрицу вещественных чисел.
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <param name="result">Матрица вещественных чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string s, ref double[,] result)
        {
            double[,] zero = new double[0, 0];
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
        /// Переводит исходную строку в матрицу комплексных чисел.
        /// </summary>
        /// <remarks>
        /// Пример входной строки: "[1 + 2i, 2 + 4i; 3 + 6i, 4 + 8i]";
        /// </remarks>
        /// <param name="a">Матрица</param>
        /// <param name="s">Исходная строка</param>
        /// <returns>Матрица комплексных чисел</returns>
        public static Complex[,] Parse(this Complex[,] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length, n = nums.Length;
            Complex[,] H = new Complex[r, n];
            int i, j;

            // collecting rows:
            for (i = 0; i < r; i++)
            {
                nums = rows[i].Split('|');
                n = nums.Length;

                for (j = 0; j < n; j++)
                {
                    H[i, j] = StringOptions.Compar(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Пробует перевести исходную строку в матрицу комплексных чисел.
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <param name="result">Матрица комплексных чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string s, ref Complex[,] result)
        {
            Complex[,] zero = new Complex[0, 0];
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
        /// Переводит исходную строку в матрицу вещественных чисел.
        /// <remarks>
        /// Пример входной строки: "[1, 2, 3, 4]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="s">Исходная строка</param>
        /// <returns>Матрица</returns>
        public static double[] Parse(this double[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
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
                throw new Exception("Входная строка имела неверный формат");
            }
        }
        /// <summary>
        /// Пробует перевести исходную строку в матрицу вещественных чисел.
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <param name="result">Матрица вещественных чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string s, ref double[] result)
        {
            double[] zero = new double[0];
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
        /// Переводит исходную строку в матрицу вещественных чисел.
        /// <remarks>
        /// Пример входной строки: "[1 + 2i, 2 + 0.3i, 3 + i, 4 - 11i]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Матрица</param>
        /// <param name="s">Исходная строка</param>
        /// <returns>Матрица</returns>
        public static Complex[] Parse(this Complex[] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums = rows[0].Split('|');
            int r = rows.Length;

            // vector?
            if (r < 2)
            {
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
                throw new Exception("Входная строка имела неверный формат");
            }
        }
        /// <summary>
        /// Пробует перевести исходную строку в матрицу вещественных чисел.
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <param name="result">Матрица вещественных чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string s, ref Complex[] result)
        {
            Complex[] zero = new Complex[0];
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
    }
    #endregion

    #region Linealg options
    /// <summary>
    /// Определяет класс оптимизаций матричных операций.
    /// </summary>
    internal class LinealgOptions
    {
        #region Private data
        /// <summary>
        /// Высота матрицы.
        /// </summary>
        private static int height;
        /// <summary>
        /// Ширина матрицы.
        /// </summary>
        private static int width;
        /// <summary>
        /// Длина матрицы.
        /// </summary>
        private static int length;
        /// <summary>
        /// Текст исключения при умножении.
        /// </summary>
        private static string exception = "Длина строки матрицы A должна быть равна длине столбца матрицы B";
        #endregion

        #region Determinant
        /// <summary>
        /// Итеративный расчет детерминанта.
        /// </summary>
        /// <param name="element">"Элемент</param>
        /// <param name="n">Радиус матрицы</param>
        /// <returns>Число двойной точности</returns>
        public unsafe static double Determinant(double* element, int n)
        {
            double* mtx_u_ii, mtx_ii_j;
            double* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
            double val, det = 1;
            int d = 0;
            // rmX указывает на (i,i) элемент на каждом шаге и называется ведущим
            for (double* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
            {
                // Ищем максимальный элемент в столбце(под ведущим)
                {
                    //Ищем максимальный элемент и его позицию
                    val = Math.Abs(*(mtx_ii_j = element));
                    for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                    {
                        if (val < Math.Abs(*mtx_u_ii))
                            val = Math.Abs(*(mtx_ii_j = mtx_u_ii));
                    }
                    //Если максимальный элемент = 0 -> матрица вырожденная
                    if (Math.Abs(val - 0) < double.Epsilon) return double.NaN;
                    //Если ведущий элемент не является максимальным - делаем перестановку строк и меняем знак определителя
                    if (mtx_ii_j != element)
                    {
                        det = -det;
                        for (mtx_u_ii = element; mtx_u_ii < mtx_ii_end; mtx_ii_j++, mtx_u_ii++)
                        {
                            val = *mtx_u_ii;
                            *mtx_u_ii = *mtx_ii_j;
                            *mtx_ii_j = val;
                        }
                    }
                }
                //Обнуляем элементы под ведущим
                for (mtx_u_ii = element + n, mtx_u_ii_j = mtx_end + n; mtx_u_ii < mtx_u_ii_j; mtx_u_ii += d)
                {
                    val = *(mtx_u_ii++) / *element;
                    for (mtx_ii_j = element + 1; mtx_ii_j < mtx_ii_end; mtx_u_ii++, mtx_ii_j++)
                        *mtx_u_ii -= *mtx_ii_j * val;
                }
                det *= *element;
            }
            return det * *element;
        }
        /// <summary>
        /// Итеративный расчет детерминанта.
        /// </summary>
        /// <param name="element">"Элемент</param>
        /// <param name="n">Радиус матрицы</param>
        /// <returns>Комплексное число</returns>
        public unsafe static Complex Determinant(Complex* element, int n)
        {
            Complex* mtx_u_ii, mtx_ii_j;
            Complex* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
            Complex val, det = (Complex)1;
            int d = 0;
            // rmX указывает на (i,i) элемент на каждом шаге и называется ведущим
            for (Complex* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
            {
                // Ищем максимальный элемент в столбце(под ведущим) 
                {
                    //Ищем максимальный элемент и его позицию
                    val = (Complex)Maths.Abs(*(mtx_ii_j = element));
                    for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                    {
                        if (val.Abs < (Maths.Abs(*mtx_u_ii)))
                            val = (Complex)Maths.Abs(*(mtx_ii_j = mtx_u_ii));
                    }
                    //Если максимальный элемент = 0 -> матрица вырожденная
                    if (Maths.Abs(val - 0) < double.Epsilon) return (Complex)double.NaN;
                    //Если ведущий элемент не является максимальным - делаем перестановку строк и меняем знак определителя
                    if (mtx_ii_j != element)
                    {
                        det = -det;
                        for (mtx_u_ii = element; mtx_u_ii < mtx_ii_end; mtx_ii_j++, mtx_u_ii++)
                        {
                            val = *mtx_u_ii;
                            *mtx_u_ii = *mtx_ii_j;
                            *mtx_ii_j = val;
                        }
                    }
                }
                //Обнуляем элементы под ведущим
                for (mtx_u_ii = element + n, mtx_u_ii_j = mtx_end + n; mtx_u_ii < mtx_u_ii_j; mtx_u_ii += d)
                {
                    val = *(mtx_u_ii++) / *element;
                    for (mtx_ii_j = element + 1; mtx_ii_j < mtx_ii_end; mtx_u_ii++, mtx_ii_j++)
                        *mtx_u_ii -= *mtx_ii_j * val;
                }
                det *= *element;
            }
            return det * *element;
        }
        #endregion

        #region Multiplication
        /// <summary>
        /// Реализует умножение матриц, представленных в виде зубчатых массивов.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] Mul(float[][] A, float[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            height = A.GetLength(0); width = B[0].GetLength(0); length = B.GetLength(0);
            float[][] C = LinealgOptions.ToJagged(new double[height, width]);

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A[i], B, C[i]);
            });

            return C;
        }
        /// <summary>
        /// Реализует умножение матриц, представленных в виде зубчатых массивов.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Mul(Complex32[][] A, Complex32[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            height = A.GetLength(0); width = B[0].GetLength(0); length = B.GetLength(0);
            Complex32[][] C = LinealgOptions.ToJagged(new Complex[height, width]);

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A[i], B, C[i]);
            });

            return C;
        }
        /// <summary>
        /// Реализует умножение матриц, представленных в виде зубчатых массивов.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Mul(Complex32[][] A, float[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            height = A.GetLength(0); width = B[0].GetLength(0); length = B.GetLength(0);
            Complex32[][] C = LinealgOptions.ToJagged(new Complex[height, width]);

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A[i], B, C[i]);
            });

            return C;
        }
        /// <summary>
        /// Реализует умножение матриц, представленных в виде зубчатых массивов.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Mul(float[][] A, Complex32[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            height = A.GetLength(0); width = B[0].GetLength(0); length = B.GetLength(0);
            Complex32[][] C = LinealgOptions.ToJagged(new Complex[height, width]);

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A[i], B, C[i]);
            });

            return C;
        }
        #endregion

        #region Convolution
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] Conv(float[][] A, float[][] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r0 = B.GetLength(0), r1 = B[0].GetLength(0);
            int r0p = r0 / 2, r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0; div = 0;
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != 0)
                                {
                                    s += A[ir][jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0;
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != 0)
                                {
                                    s += A[ir][jr] * k;
                                }
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Conv(Complex32[][] A, Complex32[][] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0), r1 = B[0].GetLength(0);
            int r0p = r0 / 2, r1p = r1 / 2;


            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != new Complex32(0, 0))
                                {
                                    s += A[ir][jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != new Complex32(0, 0))
                                {
                                    s += A[ir][jr] * k;
                                }
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Conv(float[][] A, Complex32[][] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0), r1 = B[0].GetLength(0);
            int r0p = r0 / 2, r1p = r1 / 2;


            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != new Complex32(0, 0))
                                {
                                    s += A[ir][jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != new Complex32(0, 0))
                                {
                                    s += A[ir][jr] * k;
                                }
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат двумерной дискретной свертки.
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] Conv(Complex32[][] A, float[][] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0), r1 = B[0].GetLength(0);
            int r0p = r0 / 2, r1p = r1 / 2;


            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k;
                    Complex32 s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != 0f)
                                {
                                    s += A[ir][jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k;
                    Complex32 s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i][j];

                                if (k != 0f)
                                {
                                    s += A[ir][jr] * k;
                                }
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        #endregion

        #region Convolution (separable)
        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] ConvHorizontal(float[][] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0; div = 0;
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != 0)
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0;
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != 0)
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] ConvVertical(float[][] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0; div = 0;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != 0)
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = 0;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != 0)
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }

        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvHorizontal(float[][] A, Complex32[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvVertical(float[][] A, Complex32[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }

        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvHorizontal(Complex32[][] A, Complex32[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvVertical(Complex32[][] A, Complex32[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }

        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvHorizontal(Complex32[][] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = new Complex32((float)B[j], 0);

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = new Complex32((float)B[j], 0);

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат одномерной дискретной свертки (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="B">Зубчатый массив</param>
        /// <param name="normalize">Нормализованная свертка или нет</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ConvVertical(Complex32[][] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0); div = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = new Complex32((float)B[i], 0);

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex32(0, 0))
                        {
                            s /= div;
                        }

                        H[y][x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex32 k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex32(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = new Complex32((float)B[i], 0);

                            if (k != new Complex32(0, 0))
                            {
                                s += A[ir][jr] * k;
                            }
                        }

                        H[y][x] = s;
                    }
                });
            }

            return H;
        }
        #endregion

        #region Morphology (separable)
        /// <summary>
        /// Возвращает матрицу-результат морфологии (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r1">Размер</param>
        /// <param name="threshold">Пороговое значение</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MorphHorizontal(float[][] A, int r1, int threshold)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r1p = r1 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float[] s;
                int j, x;
                int xr, yr = y;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = new float[r1];
                    xr = x - r1p;
                    ir = yr;

                    for (j = 0; j < r1; j++)
                    {
                        jr = xr + j;
                        if (jr < 0) continue; if (jr >= width) break;

                        s[j] = A[ir][jr];
                    }

                    Array.Sort(s);

                    H[y][x] = s[threshold];
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологии (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r0">Размер</param>
        /// <param name="threshold">Пороговое значение</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MorphVertical(float[][] A, int r0, int threshold)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r0p = r0 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float[] s;
                int i, x;
                int yr = y - r0p;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = new float[r0];

                    for (i = 0; i < r0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = x;

                        s[i] = A[ir][jr];
                    }

                    Array.Sort(s);

                    H[y][x] = s[threshold];
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологического сужения (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r1">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MinHorizontal(float[][] A, int r1)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r1p = r1 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float s;
                int j, x;
                int xr, yr = y;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = float.MaxValue;
                    xr = x - r1p;
                    ir = yr;

                    for (j = 0; j < r1; j++)
                    {
                        jr = xr + j;
                        if (jr < 0) continue; if (jr >= width) break;

                        s = (A[ir][jr] < s) ? A[ir][jr] : s;
                    }

                    H[y][x] = s;
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологического сужения (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r0">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MinVertical(float[][] A, int r0)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r0p = r0 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float s;
                int i, x;
                int yr = y - r0p;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = float.MaxValue;

                    for (i = 0; i < r0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = x;

                        s = (A[ir][jr] < s) ? A[ir][jr] : s;
                    }

                    H[y][x] = s;
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологического расширения (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r1">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MaxHorizontal(float[][] A, int r1)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r1p = r1 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float s;
                int j, x;
                int xr, yr = y;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = float.MinValue;
                    xr = x - r1p;
                    ir = yr;

                    for (j = 0; j < r1; j++)
                    {
                        jr = xr + j;
                        if (jr < 0) continue; if (jr >= width) break;

                        s = (A[ir][jr] > s) ? A[ir][jr] : s;
                    }

                    H[y][x] = s;
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат морфологического расширения (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r0">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MaxVertical(float[][] A, int r0)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int r0p = r0 / 2;

            // morphology:
            Parallel.For(0, height, y =>
            {
                float s;
                int i, x;
                int yr = y - r0p;
                int ir, jr;

                for (x = 0; x < width; x++)
                {
                    s = float.MinValue;

                    for (i = 0; i < r0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = x;

                        s = (A[ir][jr] > s) ? A[ir][jr] : s;
                    }

                    H[y][x] = s;
                }
            });

            return H;
        }
        #endregion

        #region Mean (separable)
        /// <summary>
        /// Возвращает матрицу-результат локального усреднения (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r1">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MeanHorizontal(float[][] A, int r1)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int h = r1 >= width ? width - 1 : r1;
            int v = h >> 1;
            int dl = width - v;
            
            Parallel.For(0, height, y =>
            {
                float s = 0;
                int x;

                // Вычисление глобальной суммы [0, h):
                for (x = 0; x < h; x++)
                {
                    s += A[y][x];
                }
                // Вычисление фильтра на отрезке [0, v):
                for (x = 0; x < v; x++)
                {
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [v, l-v):
                for (x = v; x < dl; x++)
                {
                    s = s - A[y][x - v] + A[y][x + v];
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [l-v, l):
                for (x = dl; x < width; x++)
                {
                    s = s - A[y][x - v] + A[y][x];
                    H[y][x] = s / h;
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат локального усреднения (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r0">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] MeanVertical(float[][] A, int r0)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            float[][] H = LinealgOptions.ToJagged(new double[height, width]);
            int h = r0 >= height ? height - 1 : r0;
            int v = h >> 1;
            int dl = height - v;

            Parallel.For(0, width, x =>
            {
                float s = 0;
                int y;
                
                // Вычисление глобальной суммы [0, h):
                for (y = 0; y < h; y++)
                {
                    s += A[y][x];
                }
                // Вычисление фильтра на отрезке [0, v):
                for (y = 0; y < v; y++)
                {
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [v, l-v):
                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v][x] + A[y + v][x];
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [l-v, l):
                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v][x] + A[y][x];
                    H[y][x] = s / h;
                }

            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат локального усреднения (по горизонтали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r1">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] MeanHorizontal(Complex32[][] A, int r1)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int h = r1 >= width ? width - 1 : r1;
            int v = h >> 1;
            int dl = width - v;

            Parallel.For(0, height, y =>
            {
                Complex32 s = new Complex32(0, 0);
                int x;

                // Вычисление глобальной суммы [0, h):
                for (x = 0; x < h; x++)
                {
                    s += A[y][x];
                }
                // Вычисление фильтра на отрезке [0, v):
                for (x = 0; x < v; x++)
                {
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [v, l-v):
                for (x = v; x < dl; x++)
                {
                    s = s - A[y][x - v] + A[y][x + v];
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [l-v, l):
                for (x = dl; x < width; x++)
                {
                    s = s - A[y][x - v] + A[y][x];
                    H[y][x] = s / h;
                }
            });

            return H;
        }
        /// <summary>
        /// Возвращает матрицу-результат локального усреднения (по вертикали).
        /// </summary>
        /// <param name="A">Зубчатый массив</param>
        /// <param name="r0">Размер</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] MeanVertical(Complex32[][] A, int r0)
        {
            int height = A.GetLength(0), width = A[0].GetLength(0);
            Complex32[][] H = LinealgOptions.ToJagged(new Complex[height, width]);
            int h = r0 >= height ? height - 1 : r0;
            int v = h >> 1;
            int dl = height - v;

            Parallel.For(0, width, x =>
            {
                Complex32 s = new Complex32(0, 0);
                int y;

                // Вычисление глобальной суммы [0, h):
                for (y = 0; y < h; y++)
                {
                    s += A[y][x];
                }
                // Вычисление фильтра на отрезке [0, v):
                for (y = 0; y < v; y++)
                {
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [v, l-v):
                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v][x] + A[y + v][x];
                    H[y][x] = s / h;
                }
                // Вычисление фильтра на отрезке [l-v, l):
                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v][x] + A[y][x];
                    H[y][x] = s / h;
                }

            });

            return H;
        }
        #endregion

        #region Jagged Array Options
        /// <summary>
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Зубчатый массив</returns>
        public static float[][] ToJagged(double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[][] jagged = new float[ml][];
            float[] dummy;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                dummy = new float[mr];
                for (j = 0; j < mr; j++)
                {
                    dummy[j] = (float)m[i, j];
                }
                jagged[i] = dummy;
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static double[,] FromJagged(float[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            double[,] m = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    m[i, j] = jagged[i][j];
                }
            }
            return m;
        }
        /// <summary>
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Матрица</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[][] ToJagged(Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[][] jagged = new Complex32[ml][];
            Complex32[] dummy;
            Complex mij;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                dummy = new Complex32[mr];
                for (j = 0; j < mr; j++)
                {
                    mij = m[i, j];
                    dummy[j] = new Complex32((float)mij.Re, (float)mij.Im);
                }
                jagged[i] = dummy;
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static Complex[,] FromJagged(Complex32[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            Complex[,] m = new Complex[ml, mr];
            Complex32 jaggedij;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    jaggedij = jagged[i][j];
                    m[i, j] = new Complex(jaggedij.Re, jaggedij.Im);
                }
            }
            return m;
        }

        /// <summary>
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Вектор</param>
        /// <returns>Зубчатый массив</returns>
        public static float[] ToJagged(double[] m)
        {
            int n = m.GetLength(0);
            float[] jagged = new float[n];
            int i;

            for (i = 0; i < n; i++)
            {
                jagged[i] = (float)m[i];
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static double[] FromJagged(float[] jagged)
        {
            int n = jagged.GetLength(0);
            double[] m = new double[n];
            int i;

            for (i = 0; i < n; i++)
            {
                m[i] = jagged[i];
            }
            return m;
        }
        /// <summary>
        /// Возвращает зубчатый массив.
        /// </summary>
        /// <param name="m">Вектор</param>
        /// <returns>Зубчатый массив</returns>
        public static Complex32[] ToJagged(Complex[] m)
        {
            int n = m.GetLength(0);
            Complex32[] jagged = new Complex32[n];
            Complex mi;
            int i;

            for (i = 0; i < n; i++)
            {
                mi = m[i];
                jagged[i] = new Complex32((float)mi.Re, (float)mi.Im);
            }
            return jagged;
        }
        /// <summary>
        /// Возвращает матрицу.
        /// </summary>
        /// <param name="jagged">Зубчатый массив</param>
        /// <returns>Матрица</returns>
        public static Complex[] FromJagged(Complex32[] jagged)
        {
            int n = jagged.GetLength(0);
            Complex[] m = new Complex[n];
            Complex32 mi;
            int i;

            for (i = 0; i < n; i++)
            {
                mi = jagged[i];
                m[i] = new Complex(mi.Re, mi.Im);
            }
            return m;
        }
        #endregion

        #region Complex32
        /// <summary>
        /// Определяет комплексное число, вещественные части которой представляются числами одиночной точности с плавающей запятой.
        /// </summary>
        public struct Complex32
        {
            #region Private data
            /// <summary>
            /// Действительная часть числа.
            /// </summary>
            public float Re;
            /// <summary>
            /// Мнимая часть числа.
            /// </summary>
            public float Im;
            #endregion

            #region Struct Components
            /// <summary>
            /// Инициализирует комплексное число.
            /// </summary>
            /// <param name="re">Действительная часть комплексного числа</param>
            /// <param name="im">Мнимая часть комплексного числа</param>
            public Complex32(float re, float im)
            {
                this.Re = re; this.Im = im;
            }
            #endregion

            #region Overrides
            /// <summary>
            /// Возвращает хэш-код для данного объекта.
            /// </summary>
            /// <returns>Целое число со знаком</returns>
            public override int GetHashCode()
            {
                return Re.GetHashCode() ^ Im.GetHashCode();
            }
            /// <summary>
            /// Возвращает значение, указывающее, равен ли данный экземляр заданному значению типа Complex.
            /// </summary>
            /// <param name="obj">Объект</param>
            /// <returns>Логическое значение</returns>
            public override bool Equals(object obj)
            {
                return (obj is Complex32) ? (this == (Complex32)obj) : false;
            }
            #endregion

            #region Bools
            /// <summary>
            /// Проверяет равны ли два объекта типа Comlex между собой.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Логическое значение</returns>
            public static bool operator ==(Complex32 a, Complex32 b)
            {
                return ((a.Re == b.Re) && (a.Im == b.Im));
            }
            /// <summary>
            /// Проверяет не равны ли два объекта типа Complex между собой.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Логическое значение</returns>
            public static bool operator !=(Complex32 a, Complex32 b)
            {
                return !(a == b);
            }
            #endregion

            #region Operations
            /// <summary>
            /// Сумма двух комплексных чисел
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator +(Complex32 a, Complex32 b)
            {
                return new Complex32(a.Re + b.Re, a.Im + b.Im);
            }
            /// <summary>
            /// Сумма комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator +(Complex32 a, float b)
            {
                return new Complex32(a.Re + b, a.Im);
            }
            /// <summary>
            /// Сумма комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator +(float b, Complex32 a)
            {
                return new Complex32(a.Re + b, a.Im);
            }
            /// <summary>
            /// Разность двух комплексных чисел
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator -(Complex32 a, Complex32 b)
            {
                return new Complex32(a.Re - b.Re, a.Im - b.Im);
            }
            /// <summary>
            /// Разность комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator -(Complex32 a, float b)
            {
                return new Complex32(a.Re - b, a.Im);
            }
            /// <summary>
            /// Разность комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator -(float b, Complex32 a)
            {
                return new Complex32(-a.Re + b, a.Im);
            }
            /// <summary>
            /// Произведение двух комплексных чисел
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator *(Complex32 a, Complex32 b)
            {
                float aRe = a.Re, aIm = a.Im;
                float bRe = b.Re, bIm = b.Im;

                return new Complex32(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
            }
            /// <summary>
            /// Произведение комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator *(Complex32 a, float b)
            {
                return new Complex32(a.Re * b, a.Im * b);
            }
            /// <summary>
            /// Произведение комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator *(float b, Complex32 a)
            {
                return new Complex32(a.Re * b, a.Im * b);
            }

            /// <summary>
            /// Частное двух комплексных чисел
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator /(Complex32 a, Complex32 b)
            {
                float aRe = a.Re, aIm = a.Im;
                float bRe = b.Re, bIm = b.Im;
                float abs = bRe * bRe + bIm * bIm;
                float inv = 1 / abs;

                return new Complex32((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
            }
            /// <summary>
            /// Частное комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Комплексное число</param>
            /// <param name="b">Число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator /(Complex32 a, float b)
            {
                return new Complex32(a.Re / b, a.Im / b);
            }
            /// <summary>
            /// Частное комплексного числа и действительного числа.
            /// </summary>
            /// <param name="a">Числа</param>
            /// <param name="b">Комплексное число</param>
            /// <returns>Комплексное число</returns>
            public static Complex32 operator /(float a, Complex32 b)
            {
                if (b.Im == 0)
                {
                    return new Complex32(a / b.Re, 0);
                }
                else if (b.Re == 0)
                {
                    return new Complex32(0, a / b.Im);
                }
                return new Complex32(a / b.Re, a / b.Im);
            }
            #endregion
        }
        #endregion

        #region Swap voids
        /// <summary>
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="row1">Первый ряд</param>
        /// <param name="row2">Второй ряд</param>
        /// <param name="col">Размерность столбца</param>
        public static void Swap(double[,] A, int row1, int row2, int col)
        {
            double temp;
            for (int i = 0; i < col; i++)
            {
                temp = A[row1, i];
                A[row1, i] = A[row2, i];
                A[row2, i] = temp;
            }
        }
        /// <summary>
        /// Реализует перестановку векторов матрицы.
        /// </summary>
        /// <param name="A">Матрица</param>
        /// <param name="row1">Первый ряд</param>
        /// <param name="row2">Второй ряд</param>
        /// <param name="col">Размерность столбца</param>
        public static void Swap(Complex[,] A, int row1, int row2, int col)
        {
            Complex temp;
            for (int i = 0; i < col; i++)
            {
                temp = A[row1, i];
                A[row1, i] = A[row2, i];
                A[row2, i] = temp;
            }
        }
        #endregion

        #region Modified Whittle's Multiply Optimizations
        /// <summary>
        /// Реализует умножние матриц с использованием модифицированной оптимизации Уиттла.
        /// </summary>
        /// <param name="iRowA">Строка матрицы A</param>
        /// <param name="B">Матрица B</param>
        /// <param name="iRowC">Строка матрицы C</param>
        private static void Whittle_Mul(float[] iRowA, float[][] B, float[] iRowC)
        {
            float[] kRowB;
            float ikA;
            int k, j;

            for (k = 0; k < length; k++)
            {
                kRowB = B[k]; ikA = iRowA[k];
                for (j = 0; j < width; j++)
                {
                    iRowC[j] += ikA * kRowB[j];
                }
            }
            return;
        }
        /// <summary>
        /// Реализует умножние матриц с использованием модифицированной оптимизации Уиттла.
        /// </summary>
        /// <param name="iRowA">Строка матрицы A</param>
        /// <param name="B">Матрица B</param>
        /// <param name="iRowC">Строка матрицы C</param>
        private static void Whittle_Mul(Complex32[] iRowA, Complex32[][] B, Complex32[] iRowC)
        {
            Complex32[] kRowB;
            Complex32 ikA;
            int k, j;

            for (k = 0; k < length; k++)
            {
                kRowB = B[k]; ikA = iRowA[k];

                for (j = 0; j < width; j++)
                {
                    iRowC[j] += ikA * kRowB[j];
                }
            }
        }
        /// <summary>
        /// Реализует умножние матриц с использованием модифицированной оптимизации Уиттла.
        /// </summary>
        /// <param name="iRowA">Строка матрицы A</param>
        /// <param name="B">Матрица B</param>
        /// <param name="iRowC">Строка матрицы C</param>
        private static void Whittle_Mul(Complex32[] iRowA, float[][] B, Complex32[] iRowC)
        {
            float[] kRowB;
            Complex32 ikA;
            int k, j;

            for (k = 0; k < length; k++)
            {
                kRowB = B[k]; ikA = iRowA[k];

                for (j = 0; j < width; j++)
                {
                    iRowC[j] += ikA * kRowB[j];
                }
            }
        }
        /// <summary>
        /// Реализует умножние матриц с использованием модифицированной оптимизации Уиттла.
        /// </summary>
        /// <param name="iRowA">Строка матрицы A</param>
        /// <param name="B">Матрица B</param>
        /// <param name="iRowC">Строка матрицы C</param>
        private static void Whittle_Mul(float[] iRowA, Complex32[][] B, Complex32[] iRowC)
        {
            Complex32[] kRowB;
            float ikA;
            int k, j;

            for (k = 0; k < length; k++)
            {
                kRowB = B[k]; ikA = iRowA[k];

                for (j = 0; j < width; j++)
                {
                    iRowC[j] += ikA * kRowB[j];
                }
            }
        }
        #endregion
    }
    #endregion

    #region Enums
    /// <summary>
    /// Определяет направление обработки.
    /// </summary>
    public enum Direction
    {
        #region Direction
        /// <summary>
        /// Горизонтальное направление.
        /// </summary>
        Horizontal,
        /// <summary>
        /// Вертикальное направление.
        /// </summary>
        Vertical,
        /// <summary>
        /// Двустороннее направление.
        /// </summary>
        Both,
        #endregion
    }
    #endregion
}