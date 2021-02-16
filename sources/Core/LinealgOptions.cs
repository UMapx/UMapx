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
        /// <returns>float precision floating point number</returns>
        public unsafe static float Determinant(float* element, int n)
        {
            float* mtx_u_ii, mtx_ii_j;
            float* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
            float val, det = 1;
            int d = 0;

            for (float* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
            {
                {
                    val = Math.Abs(*(mtx_ii_j = element));
                    for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                    {
                        if (val < Math.Abs(*mtx_u_ii))
                            val = Math.Abs(*(mtx_ii_j = mtx_u_ii));
                    }

                    if (Math.Abs(val - 0) < float.Epsilon) return float.NaN;

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

                    if (Maths.Abs(val - 0) < float.Epsilon) return (Complex)float.NaN;

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

        #region Inversion
        /// <summary>
        /// Implements the matrix inversion operation.
        /// </summary>
        /// <param name="working">Square matrix</param>
        /// <returns>Matrix</returns>
        public static float[][] Invert(float[][] working)
        {
            // There are faster ways to do this, but for simplicity
            // and to get something working quickly, I'll just write a 
            // simple Gaussian-elimination matrix inverter here, and look 
            // at speeding things up later.
            // https://gswce.net/?p=585

            float epsilon = 1e-32f;

            // This routine destroys the matrix it's working on, so I'll first 
            // make a copy of it, as well as setting up the output matrix invA 
            // as a unit matrix of the appropriate size:
            int dimension = working.GetLength(0);
            float[][] inverse = Jagged.Zero(dimension, dimension);
            // C# will set the initial values to zero, so to create a unit
            // matrix, I just need to fill in the diagonal elements:
            for (int loop = 0; loop < dimension; loop++) inverse[loop][loop] = 1.0f;

            // OK, first convert working to upper triangular form:
            for (int loop = 0; loop < dimension; loop++) // for each row
            {
                int currentRow = loop;

                // First step is pivoting: make sure the biggest element
                // remaining in any column is on the next row.  First, find
                // the biggest element remaining in the current column:
                float biggestSoFar = 0.0f; int biggestRow = currentRow;
                for (int x = currentRow; x < dimension; x++)
                {
                    float sizeOfThis = (float)Maths.Abs(working[x][currentRow]);
                    if (sizeOfThis > biggestSoFar)
                    {
                        biggestSoFar = sizeOfThis;
                        biggestRow = x;
                    }
                }

                // and if this is not at the top, swop the rows of working
                // and inverse around until it is:
                if (biggestRow != currentRow)
                {
                    float temp;
                    for (int lop = currentRow; lop < dimension; lop++)
                    {
                        temp = working[currentRow][lop];
                        working[currentRow][lop] = working[biggestRow][lop];
                        working[biggestRow][lop] = temp;
                    }
                    for (int lop = 0; lop < dimension; lop++)
                    {
                        temp = inverse[currentRow][lop];
                        inverse[currentRow][lop] = inverse[biggestRow][lop];
                        inverse[biggestRow][lop] = temp;
                    }
                }

                // Then, go down the matrix subtracting as necessary
                // to get rid of the lower-triangular elements:
                for (int lop = currentRow + 1; lop < dimension; lop++)
                {
                    // Matrix might be ill-conditioned.  I should check:
                    if (working[currentRow][currentRow] == 0)
                    {
                        working[currentRow][currentRow] = epsilon;
                    }
                    float factor = working[lop][currentRow] / working[currentRow][currentRow];

                    // If the matrix is fairly sparse (quite common for this
                    // application), it might make sense to check that the 
                    // lower elements are not already zero before doing all
                    // the scaling and replacing:
                    if (factor != 0.0)
                    {
                        // Only have to do from current row on in working, but due
                        // to pivoting, might have to do the entire row in inverse:
                        for (int lp = currentRow; lp < dimension; lp++)
                            working[lop][lp] -= factor * working[currentRow][lp];
                        for (int lp = 0; lp < dimension; lp++)
                            inverse[lop][lp] -= factor * inverse[currentRow][lp];
                    }
                }
                // That's it for this row, now on to the next one...
            }

            // Now with the working matrix in upper-triangular form, continue the same
            // process amongst the upper-triangular elements to convert working into
            // diagonal form:
            for (int loop = dimension - 1; loop >= 0; loop--) // for each row
            {
                int currentRow = loop;

                // Matrix might be ill-conditioned.  I should check:
                if (working[currentRow][currentRow] == 0)
                {
                    working[currentRow][currentRow] = epsilon;
                }

                // Then, go up the matrix subtracting as necessary to get 
                // rid of the remaining upper-triangular elements:
                for (int lop = currentRow - 1; lop >= 0; lop--)
                {
                    float factor = working[lop][currentRow] / working[currentRow][currentRow];

                    // There's only one element in working to change (the other elements
                    // in the row of working are all zero), and that will always be set
                    // to zero; but you might have to do the entire row in inverse:
                    working[lop][currentRow] = 0.0f;

                    if (factor != 0.0f)
                    {
                        for (int lp = 0; lp < dimension; lp++)
                        {
                            inverse[lop][lp] -= factor * inverse[currentRow][lp];
                        }
                    }
                }
                // That's it for this row, now on to the next one...
            }

            // Should now have working as a diagonal matrix.  Final thing is 
            // to scale all the rows:
            for (int loop = 0; loop < dimension; loop++)
            {
                float scale = working[loop][loop];
                for (int lop = 0; lop < dimension; lop++) inverse[loop][lop] /= scale;
            }

            // That's it.  inverse should now be the inverse of the original matrix.
            return inverse;
        }
        /// <summary>
        /// Implements the matrix inversion operation.
        /// </summary>
        /// <param name="working">Square matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Invert(Complex[][] working)
        {
            // There are faster ways to do this, but for simplicity
            // and to get something working quickly, I'll just write a 
            // simple Gaussian-elimination matrix inverter here, and look 
            // at speeding things up later.
            // https://gswce.net/?p=585

            float epsilon = 1e-32f;

            // This routine destroys the matrix it's working on, so I'll first 
            // make a copy of it, as well as setting up the output matrix invA 
            // as a unit matrix of the appropriate size:
            int dimension = working.GetLength(0);
            Complex[][] inverse = Jagged.ToComplex(Jagged.Zero(dimension, dimension));
            // C# will set the initial values to zero, so to create a unit
            // matrix, I just need to fill in the diagonal elements:
            for (int loop = 0; loop < dimension; loop++) inverse[loop][loop] = new Complex(1.0f, 0);

            // OK, first convert working to upper triangular form:
            for (int loop = 0; loop < dimension; loop++) // for each row
            {
                int currentRow = loop;

                // First step is pivoting: make sure the biggest element
                // remaining in any column is on the next row.  First, find
                // the biggest element remaining in the current column:
                float biggestSoFar = 0.0f; int biggestRow = currentRow;
                for (int x = currentRow; x < dimension; x++)
                {
                    float sizeOfThis = working[x][currentRow].Abs;
                    if (sizeOfThis > biggestSoFar)
                    {
                        biggestSoFar = sizeOfThis;
                        biggestRow = x;
                    }
                }

                // and if this is not at the top, swop the rows of working
                // and inverse around until it is:
                if (biggestRow != currentRow)
                {
                    Complex temp;
                    for (int lop = currentRow; lop < dimension; lop++)
                    {
                        temp = working[currentRow][lop];
                        working[currentRow][lop] = working[biggestRow][lop];
                        working[biggestRow][lop] = temp;
                    }
                    for (int lop = 0; lop < dimension; lop++)
                    {
                        temp = inverse[currentRow][lop];
                        inverse[currentRow][lop] = inverse[biggestRow][lop];
                        inverse[biggestRow][lop] = temp;
                    }
                }

                // Then, go down the matrix subtracting as necessary
                // to get rid of the lower-triangular elements:
                for (int lop = currentRow + 1; lop < dimension; lop++)
                {
                    // Matrix might be ill-conditioned.  I should check:
                    if (working[currentRow][currentRow] == new Complex(0, 0))
                    {
                        working[currentRow][currentRow] = new Complex(epsilon, 0);
                    }
                    Complex factor = working[lop][currentRow] / working[currentRow][currentRow];

                    // If the matrix is fairly sparse (quite common for this
                    // application), it might make sense to check that the 
                    // lower elements are not already zero before doing all
                    // the scaling and replacing:
                    if (factor != new Complex(0, 0))
                    {
                        // Only have to do from current row on in working, but due
                        // to pivoting, might have to do the entire row in inverse:
                        for (int lp = currentRow; lp < dimension; lp++)
                            working[lop][lp] -= factor * working[currentRow][lp];
                        for (int lp = 0; lp < dimension; lp++)
                            inverse[lop][lp] -= factor * inverse[currentRow][lp];
                    }
                }
                // That's it for this row, now on to the next one...
            }

            // Now with the working matrix in upper-triangular form, continue the same
            // process amongst the upper-triangular elements to convert working into
            // diagonal form:
            for (int loop = dimension - 1; loop >= 0; loop--) // for each row
            {
                int currentRow = loop;

                // Matrix might be ill-conditioned.  I should check:
                if (working[currentRow][currentRow] == new Complex(0, 0))
                {
                    working[currentRow][currentRow] = new Complex(epsilon, 0);
                }

                // Then, go up the matrix subtracting as necessary to get 
                // rid of the remaining upper-triangular elements:
                for (int lop = currentRow - 1; lop >= 0; lop--)
                {
                    Complex factor = working[lop][currentRow] / working[currentRow][currentRow];

                    // There's only one element in working to change (the other elements
                    // in the row of working are all zero), and that will always be set
                    // to zero; but you might have to do the entire row in inverse:
                    working[lop][currentRow] = new Complex(0, 0);

                    if (factor != new Complex(0, 0))
                    {
                        for (int lp = 0; lp < dimension; lp++)
                        {
                            inverse[lop][lp] -= factor * inverse[currentRow][lp];
                        }
                    }
                }
                // That's it for this row, now on to the next one...
            }

            // Should now have working as a diagonal matrix.  Final thing is 
            // to scale all the rows:
            for (int loop = 0; loop < dimension; loop++)
            {
                Complex scale = working[loop][loop];
                for (int lop = 0; lop < dimension; lop++) inverse[loop][lop] /= scale;
            }

            // That's it.  inverse should now be the inverse of the original matrix.
            return inverse;
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
            float[][] C = Jagged.Zero(height, width);

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A, B, C, length, width, i);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex[][] Mul(Complex[][] A, Complex[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex[][] C = Jagged.Zero(height, width).ToComplex();

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A, B, C, length, width, i);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex[][] Mul(Complex[][] A, float[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex[][] C = Jagged.Zero(height, width).ToComplex();

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A, B, C, length, width, i);
            });

            return C;
        }
        /// <summary>
        /// Implements the multiplication of matrices presented in the form of jagged arrays.
        /// </summary>
        /// <param name="A">Jagged array</param>
        /// <param name="B">Jagged array</param>
        /// <returns>Jagged array</returns>
        public static Complex[][] Mul(float[][] A, Complex[][] B)
        {
            if (A[0].GetLength(0) != B.GetLength(0))
                throw new Exception(exception);

            int height = A.GetLength(0);
            int width = B[0].GetLength(0);
            int length = B.GetLength(0);
            Complex[][] C = Jagged.Zero(height, width).ToComplex();

            Parallel.For(0, height, i =>
            {
                LinealgOptions.Whittle_Mul(A, B, C, length, width, i);
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
        public static float[,] Conv(float[,] A, float[,] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
            int r0 = B.GetLength(0), r1 = B.GetLength(1);
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

                                k = B[i,j];

                                if (k != 0)
                                {
                                    s += A[ir,jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y,x] = s;
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

                                k = B[i,j];

                                if (k != 0)
                                {
                                    s += A[ir,jr] * k;
                                }
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] Conv(Complex[,] A, Complex[,] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0), r1 = B.GetLength(1);
            int r0p = r0 / 2, r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != new Complex(0, 0))
                                {
                                    s += A[ir,jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != new Complex(0, 0))
                                {
                                    s += A[ir,jr] * k;
                                }
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] Conv(float[,] A, Complex[,] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0), r1 = B.GetLength(1);
            int r0p = r0 / 2, r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != new Complex(0, 0))
                                {
                                    s += A[ir,jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != new Complex(0, 0))
                                {
                                    s += A[ir,jr] * k;
                                }
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] Conv(Complex[,] A, float[,] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0), r1 = B.GetLength(1);
            int r0p = r0 / 2, r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k;
                    Complex s, div;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != 0f)
                                {
                                    s += A[ir,jr] * k;
                                    div += k;
                                }
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    float k;
                    Complex s;
                    int i, j, x;
                    int xr, yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[i,j];

                                if (k != 0f)
                                {
                                    s += A[ir,jr] * k;
                                }
                            }
                        }

                        H[y,x] = s;
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
        public static float[,] ConvHorizontal(float[,] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y,x] = s;
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
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static float[,] ConvVertical(float[,] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != 0)
                        {
                            s /= div;
                        }

                        H[y,x] = s;
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
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvHorizontal(float[,] A, Complex[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvVertical(float[,] A, Complex[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvHorizontal(Complex[,] A, Complex[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = B[j];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvVertical(Complex[,] A, Complex[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = B[i];

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvHorizontal(Complex[,] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r1 = B.GetLength(0);
            int r1p = r1 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = new Complex(B[j], 0);

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int j, x;
                    int xr, yr = y;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);
                        xr = x - r1p;
                        ir = yr;

                        for (j = 0; j < r1; j++)
                        {
                            jr = xr + j;
                            if (jr < 0) continue; if (jr >= width) break;

                            k = new Complex((float)B[j], 0);

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static Complex[,] ConvVertical(Complex[,] A, float[] B, bool normalize = true)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int r0 = B.GetLength(0);
            int r0p = r0 / 2;

            if (normalize)
            {
                // normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s, div;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0); div = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = new Complex(B[i], 0);

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                                div += k;
                            }
                        }

                        if (div != new Complex(0, 0))
                        {
                            s /= div;
                        }

                        H[y,x] = s;
                    }
                });
            }
            else
            {
                // non-normalize convolution:
                Parallel.For(0, height, y =>
                {
                    Complex k, s;
                    int i, x;
                    int yr = y - r0p;
                    int ir, jr;

                    for (x = 0; x < width; x++)
                    {
                        s = new Complex(0, 0);

                        for (i = 0; i < r0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            jr = x;

                            k = new Complex(B[i], 0);

                            if (k != new Complex(0, 0))
                            {
                                s += A[ir,jr] * k;
                            }
                        }

                        H[y,x] = s;
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
        public static float[,] MorphHorizontal(float[,] A, int r1, int threshold)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s[j] = A[ir,jr];
                    }

                    Array.Sort(s);

                    H[y,x] = s[threshold];
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
        public static float[,] MorphVertical(float[,] A, int r0, int threshold)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s[i] = A[ir,jr];
                    }

                    Array.Sort(s);

                    H[y,x] = s[threshold];
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
        public static float[,] MinHorizontal(float[,] A, int r1)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s = (A[ir,jr] < s) ? A[ir,jr] : s;
                    }

                    H[y,x] = s;
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
        public static float[,] MinVertical(float[,] A, int r0)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s = (A[ir,jr] < s) ? A[ir,jr] : s;
                    }

                    H[y,x] = s;
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
        public static float[,] MaxHorizontal(float[,] A, int r1)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s = (A[ir,jr] > s) ? A[ir,jr] : s;
                    }

                    H[y,x] = s;
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
        public static float[,] MaxVertical(float[,] A, int r0)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
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

                        s = (A[ir,jr] > s) ? A[ir,jr] : s;
                    }

                    H[y,x] = s;
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
        public static float[,] MeanHorizontal(float[,] A, int r1)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
            int h = r1 >= width ? width - 1 : r1;
            int v = h >> 1;
            int dl = width - v;

            Parallel.For(0, height, y =>
            {
                float s = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    s += A[y, x];
                }

                for (x = 0; x < v; x++)
                {
                    H[y,x] = s / h;
                }

                for (x = v; x < dl; x++)
                {
                    s = s - A[y, x - v] + A[y, x + v];
                    H[y, x] = s / h;
                }

                for (x = dl; x < width; x++)
                {
                    s = s - A[y, x - v] + A[y, x];
                    H[y, x] = s / h;
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
        public static float[,] MeanVertical(float[,] A, int r0)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            float[,] H = new float[height, width];
            int h = r0 >= height ? height - 1 : r0;
            int v = h >> 1;
            int dl = height - v;

            Parallel.For(0, width, x =>
            {
                float s = 0;
                int y;

                for (y = 0; y < h; y++)
                {
                    s += A[y, x];
                }

                for (y = 0; y < v; y++)
                {
                    H[y, x] = s / h;
                }

                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v, x] + A[y + v, x];
                    H[y, x] = s / h;
                }

                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v, x] + A[y, x];
                    H[y, x] = s / h;
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
        public static Complex[,] MeanHorizontal(Complex[,] A, int r1)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int h = r1 >= width ? width - 1 : r1;
            int v = h >> 1;
            int dl = width - v;

            Parallel.For(0, height, y =>
            {
                Complex s = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    s += A[y, x];
                }

                for (x = 0; x < v; x++)
                {
                    H[y, x] = s / h;
                }

                for (x = v; x < dl; x++)
                {
                    s = s - A[y, x - v] + A[y, x + v];
                    H[y, x] = s / h;
                }

                for (x = dl; x < width; x++)
                {
                    s = s - A[y, x - v] + A[y, x];
                    H[y, x] = s / h;
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
        public static Complex[,] MeanVertical(Complex[,] A, int r0)
        {
            int height = A.GetLength(0), width = A.GetLength(1);
            Complex[,] H = new Complex[height, width];
            int h = r0 >= height ? height - 1 : r0;
            int v = h >> 1;
            int dl = height - v;

            Parallel.For(0, width, x =>
            {
                Complex s = 0;
                int y;

                for (y = 0; y < h; y++)
                {
                    s += A[y, x];
                }

                for (y = 0; y < v; y++)
                {
                    H[y, x] = s / h;
                }

                for (y = v; y < dl; y++)
                {
                    s = s - A[y - v, x] + A[y + v, x];
                    H[y, x] = s / h;
                }

                for (y = dl; y < height; y++)
                {
                    s = s - A[y - v, x] + A[y, x];
                    H[y, x] = s / h;
                }

            });

            return H;
        }
        #endregion

        #region Modified Whittle multiply optimizations
        /// <summary>
        /// Implements matrix multiplication using modified Whittle optimization.
        /// </summary>
        /// <param name="A">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="C">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        /// <param name="i">Index</param>
        private static void Whittle_Mul(float[][] A, float[][] B, float[][] C, int length, int width, int i)
        {
            float[] iRowA = A[i];
            float[] iRowC = C[i];
            int k, j;

            for (k = 0; k < length; k++)
            {
                float[] kRowB = B[k];
                float ikA = iRowA[k];

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
        /// <param name="A">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="C">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        /// <param name="i">Index</param>
        private static void Whittle_Mul(Complex[][] A, Complex[][] B, Complex[][] C, int length, int width, int i)
        {
            Complex[] iRowA = A[i];
            Complex[] iRowC = C[i];
            int k, j;

            for (k = 0; k < length; k++)
            {
                Complex[] kRowB = B[k];
                Complex ikA = iRowA[k];

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
        /// <param name="A">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="C">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        /// <param name="i">Index</param>
        private static void Whittle_Mul(Complex[][] A, float[][] B, Complex[][] C, int length, int width, int i)
        {
            Complex[] iRowA = A[i];
            Complex[] iRowC = C[i];
            int k, j;

            for (k = 0; k < length; k++)
            {
                float[] kRowB = B[k];
                Complex ikA = iRowA[k];

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
        /// <param name="A">Row of A</param>
        /// <param name="B">Matrix B</param>
        /// <param name="C">Row of C</param>
        /// <param name="length">Length</param>
        /// <param name="width">Width</param>
        /// <param name="i">Index</param>
        private static void Whittle_Mul(float[][] A, Complex[][] B, Complex[][] C, int length, int width, int i)
        {
            float[] iRowA = A[i];
            Complex[] iRowC = C[i];
            int k, j;

            for (k = 0; k < length; k++)
            {
                Complex[] kRowB = B[k];
                Complex ikA = iRowA[k];

                for (j = 0; j < width; j++)
                {
                    iRowC[j] += ikA * kRowB[j];
                }
            }
            return;
        }
        #endregion
    }
}
