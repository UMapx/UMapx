using System;
using System.Threading.Tasks;

namespace UMapx.Core
{
    /// <summary>
    /// Defines the class of optimizations of matrix operations.
    /// </summary>
    internal class LinealgOptions
    {
        #region Private data
        /// <summary>
        /// 
        /// </summary>
        private static string exception = "The length of the row of matrix A must be equal to the length of the column of matrix B";
        #endregion

        #region Determinant
        /// <summary>
        /// Iterative calculation of the determinant.
        /// </summary>
        /// <param name="element">"Element</param>
        /// <param name="n">Radius</param>
        /// <returns>Double precision floating point number</returns>
        public unsafe static double Determinant(double* element, int n)
        {
            double* mtx_u_ii, mtx_ii_j;
            double* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
            double val, det = 1;
            int d = 0;

            for (double* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
            {
                {
                    val = Math.Abs(*(mtx_ii_j = element));
                    for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                    {
                        if (val < Math.Abs(*mtx_u_ii))
                            val = Math.Abs(*(mtx_ii_j = mtx_u_ii));
                    }

                    if (Math.Abs(val - 0) < double.Epsilon) return double.NaN;

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
        /// Iterative calculation of the determinant.
        /// </summary>
        /// <param name="element">"Element</param>
        /// <param name="n">Radius</param>
        /// <returns>Complex number</returns>
        public unsafe static Complex Determinant(Complex* element, int n)
        {
            Complex* mtx_u_ii, mtx_ii_j;
            Complex* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
            Complex val, det = (Complex)1;
            int d = 0;

            for (Complex* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
            {
                {

                    val = (Complex)Maths.Abs(*(mtx_ii_j = element));
                    for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                    {
                        if (val.Abs < (Maths.Abs(*mtx_u_ii)))
                            val = (Complex)Maths.Abs(*(mtx_ii_j = mtx_u_ii));
                    }

                    if (Maths.Abs(val - 0) < double.Epsilon) return (Complex)double.NaN;

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
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static float[][] Mul(float[][] A, float[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            float[][] C = new float[height][];

            Parallel.For(0, height, i =>
            {
                C[i] = new float[width];
                LinealgOptions.Whittle_Mul(A[i], B, C[i], length, width);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex32[][] Mul(Complex32[][] A, Complex32[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex32[][] C = new Complex32[height][];

            Parallel.For(0, height, i =>
            {
                C[i] = new Complex32[width];
                LinealgOptions.Whittle_Mul(A[i], B, C[i], length, width);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex32[][] Mul(Complex32[][] A, float[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex32[][] C = new Complex32[height][];

            Parallel.For(0, height, i =>
            {
                C[i] = new Complex32[width];
                LinealgOptions.Whittle_Mul(A[i], B, C[i], length, width);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex32[][] Mul(float[][] A, Complex32[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex32[][] C = new Complex32[height][];

            Parallel.For(0, height, i =>
            {
                C[i] = new Complex32[width];
                LinealgOptions.Whittle_Mul(A[i], B, C[i], length, width);
            });

            return C;
        }
        #endregion

        #region Convolution
        /// <summary>
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete convolution of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <param name="normalize">Normalized convolution or not</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r1">Size</param>
        /// <param name="threshold">Threshold</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r0">Size</param>
        /// <param name="threshold">Threshold</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology minimum of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r1">Size</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology minimum of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r0">Size</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology maximum of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r1">Size</param>
        /// <returns>Jagged array</returns>
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
        /// Implements discrete morphology maximum of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r0">Size</param>
        /// <returns>Jagged array</returns>
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
        /// Implements local average of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r1">Size</param>
        /// <returns>Jagged array</returns>
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

                for (x = 0; x < h; x++)
                {
                    s += A[y][x];
                }

                for (x = 0; x < v; x++)
                {
                    H[y][x] = s / h;
                }

                for (x = v; x < dl; x++)
                {
                    s = s - A[y][x - v] + A[y][x + v];
                    H[y][x] = s / h;
                }

                for (x = dl; x < width; x++)
                {
                    s = s - A[y][x - v] + A[y][x];
                    H[y][x] = s / h;
                }
            });

            return H;
        }
        /// <summary>
        /// Implements local average of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r0">Size</param>
        /// <returns>Jagged array</returns>
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

                for (y = 0; y < h; y++)
                {
                    s += A[y][x];
                }

                for (y = 0; y < v; y++)
                {
                    H[y][x] = s / h;
                }

                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v][x] + A[y + v][x];
                    H[y][x] = s / h;
                }

                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v][x] + A[y][x];
                    H[y][x] = s / h;
                }

            });

            return H;
        }
        /// <summary>
        /// Implements local average of matrices (horizontal).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r1">Size</param>
        /// <returns>Jagged array</returns>
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


                for (x = 0; x < h; x++)
                {
                    s += A[y][x];
                }

                for (x = 0; x < v; x++)
                {
                    H[y][x] = s / h;
                }

                for (x = v; x < dl; x++)
                {
                    s = s - A[y][x - v] + A[y][x + v];
                    H[y][x] = s / h;
                }

                for (x = dl; x < width; x++)
                {
                    s = s - A[y][x - v] + A[y][x];
                    H[y][x] = s / h;
                }
            });

            return H;
        }
        /// <summary>
        /// Implements local average of matrices (vertical).
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="r0">Size</param>
        /// <returns>Jagged array</returns>
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


                for (y = 0; y < h; y++)
                {
                    s += A[y][x];
                }

                for (y = 0; y < v; y++)
                {
                    H[y][x] = s / h;
                }

                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v][x] + A[y + v][x];
                    H[y][x] = s / h;
                }

                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v][x] + A[y][x];
                    H[y][x] = s / h;
                }

            });

            return H;
        }
        #endregion

        #region Jagged array
        /// <summary>
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Jagged array</returns>
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
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
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
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Jagged array</returns>
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
                    dummy[j] = new Complex32((float)mij.Real, (float)mij.Imag);
                }
                jagged[i] = dummy;
            }
            return jagged;
        }
        /// <summary>
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
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
                    m[i, j] = new Complex(jaggedij.Real, jaggedij.Imag);
                }
            }
            return m;
        }

        /// <summary>
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Array</param>
        /// <returns>Jagged array</returns>
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
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
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
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Array</param>
        /// <returns>Jagged array</returns>
        public static Complex32[] ToJagged(Complex[] m)
        {
            int n = m.GetLength(0);
            Complex32[] jagged = new Complex32[n];
            Complex mi;
            int i;

            for (i = 0; i < n; i++)
            {
                mi = m[i];
                jagged[i] = new Complex32((float)mi.Real, (float)mi.Imag);
            }
            return jagged;
        }
        /// <summary>
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
        public static Complex[] FromJagged(Complex32[] jagged)
        {
            int n = jagged.GetLength(0);
            Complex[] m = new Complex[n];
            Complex32 mi;
            int i;

            for (i = 0; i < n; i++)
            {
                mi = jagged[i];
                m[i] = new Complex(mi.Real, mi.Imag);
            }
            return m;
        }
        #endregion

        #region Complex32
        /// <summary>
        /// Complex (32 bit).
        /// </summary>
        public struct Complex32
        {
            #region Private data
            /// <summary>
            /// Real part.
            /// </summary>
            public float Real;
            /// <summary>
            /// Imaginary part.
            /// </summary>
            public float Imag;
            #endregion

            #region Struct Components
            /// <summary>
            /// Complex (32 bit).
            /// </summary>
            /// <param name="re">Real part of the complex number</param>
            /// <param name="im">Imaginary part of a complex number</param>
            public Complex32(float re, float im)
            {
                this.Real = re; this.Imag = im;
            }
            #endregion

            #region Overrides
            /// <summary>
            /// Returns the hash code for this object.
            /// </summary>
            /// <returns>Integer number</returns>
            public override int GetHashCode()
            {
                return Real.GetHashCode() ^ Imag.GetHashCode();
            }
            /// <summary>
            /// Gets a value indicating whether this instance is equal to the given value of type Complex32.
            /// </summary>
            /// <param name="obj">Object</param>
            /// <returns>Boolean</returns>
            public override bool Equals(object obj)
            {
                return (obj is Complex32) ? (this == (Complex32)obj) : false;
            }
            #endregion

            #region Bools
            /// <summary>
            /// Checks if two objects of type Complex are equal to each other.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Boolean</returns>
            public static bool operator ==(Complex32 a, Complex32 b)
            {
                return ((a.Real == b.Real) && (a.Imag == b.Imag));
            }
            /// <summary>
            /// Checks if two objects of type Complex are not equal to each other.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Boolean</returns>
            public static bool operator !=(Complex32 a, Complex32 b)
            {
                return !(a == b);
            }
            #endregion

            #region Operators
            /// <summary>
            /// The sum of two complex numbers.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator +(Complex32 a, Complex32 b)
            {
                return new Complex32(a.Real + b.Real, a.Imag + b.Imag);
            }
            /// <summary>
            /// The sum of a complex number and a real number.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator +(Complex32 a, float b)
            {
                return new Complex32(a.Real + b, a.Imag);
            }
            /// <summary>
            /// The sum of a complex number and a real number.
            /// </summary>
            /// <param name="a">Number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator +(float a, Complex32 b)
            {
                return new Complex32(b.Real + a, b.Imag);
            }


            /// <summary>
            /// The difference of two complex numbers.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator -(Complex32 a, Complex32 b)
            {
                return new Complex32(a.Real - b.Real, a.Imag - b.Imag);
            }
            /// <summary>
            /// The difference between a complex number and a real number.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator -(Complex32 a, float b)
            {
                return new Complex32(a.Real - b, a.Imag);
            }
            /// <summary>
            /// The difference between a complex number and a real number.
            /// </summary>
            /// <param name="a">Number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator -(float a, Complex32 b)
            {
                return new Complex32(a - b.Real, b.Imag);
            }
            /// <summary>
            /// Inverts complex number.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator -(Complex32 a)
            {
                return new Complex32(-a.Real, -a.Imag);
            }


            /// <summary>
            /// Multiplies one complex number by another.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator *(Complex32 a, Complex32 b)
            {
                float aRe = a.Real, aIm = a.Imag;
                float bRe = b.Real, bIm = b.Imag;

                return new Complex32(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
            }
            /// <summary>
            /// Multiplies real number by complex number.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator *(float a, Complex32 b)
            {
                return new Complex32(b.Real * a, b.Imag * a);
            }
            /// <summary>
            /// Multiplies complex number by real number.
            /// </summary>
            /// <param name="a">Number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator *(Complex32 a, float b)
            {
                return new Complex32(a.Real * b, a.Imag * b);
            }


            /// <summary>
            /// Divides one complex number by another.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator /(Complex32 a, Complex32 b)
            {
                float aRe = a.Real, aIm = a.Imag;
                float bRe = b.Real, bIm = b.Imag;
                float abs = bRe * bRe + bIm * bIm;
                float inv = 1 / abs;

                return new Complex32((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
            }
            /// <summary>
            /// Divides complex number by real number.
            /// </summary>
            /// <param name="a">Complex number</param>
            /// <param name="b">Number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator /(Complex32 a, float b)
            {
                return new Complex32(a.Real / b, a.Imag / b);
            }
            /// <summary>
            /// Divides real number by complex number.
            /// </summary>
            /// <param name="a">Number</param>
            /// <param name="b">Complex number</param>
            /// <returns>Complex number</returns>
            public static Complex32 operator /(float a, Complex32 b)
            {
                if (b.Imag == 0)
                {
                    return new Complex32(a / b.Real, 0);
                }
                else if (b.Real == 0)
                {
                    return new Complex32(0, a / b.Imag);
                }
                return new Complex32(a / b.Real, a / b.Imag);
            }
            #endregion
        }
        #endregion

        #region Modified Whittle multiply optimizations
        /// <summary>
        /// Implements matrix multiplication using modified Whittle optimization.
        /// </summary>
        /// <param name="iRowA">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="iRowC">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        private static void Whittle_Mul(float[] iRowA, float[][] B, float[] iRowC, int length, int width)
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
        /// Implements matrix multiplication using modified Whittle optimization.
        /// </summary>
        /// <param name="iRowA">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="iRowC">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        private static void Whittle_Mul(Complex32[] iRowA, Complex32[][] B, Complex32[] iRowC, int length, int width)
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
        /// Implements matrix multiplication using modified Whittle optimization.
        /// </summary>
        /// <param name="iRowA">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="iRowC">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        private static void Whittle_Mul(Complex32[] iRowA, float[][] B, Complex32[] iRowC, int length, int width)
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
        /// Implements matrix multiplication using modified Whittle optimization.
        /// </summary>
        /// <param name="iRowA">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="iRowC">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        private static void Whittle_Mul(float[] iRowA, Complex32[][] B, Complex32[] iRowC, int length, int width)
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
}
