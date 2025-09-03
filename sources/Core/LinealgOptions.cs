using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace UMapx.Core
{
    /// <summary>
    /// Defines the class of optimizations of matrix operations.
    /// </summary>
    internal static class LinealgOptions
    {
        #region Operation
        /// <summary>
        /// Defines matrix operation class. 
        /// </summary>
        public static class MatrixOperation
        {
            #region Private data
            /// <summary>
            /// Exception message.
            /// </summary>
            private static readonly string exception = "The length of the row of matrix A must be equal to the length of the column of matrix B";
            #endregion

            #region Determinant
            /// <summary>
            /// Iterative calculation of the determinant.
            /// </summary>
            /// <param name="element">"Element</param>
            /// <param name="n">Radius</param>
            /// <returns>Value</returns>
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
            public unsafe static ComplexF Determinant(ComplexF* element, int n)
            {
                ComplexF* mtx_u_ii, mtx_ii_j;
                ComplexF* mtx_end = element + n * (n - 1), mtx_u_ii_j = null;
                ComplexF val, det = (ComplexF)1;
                int d = 0;

                for (ComplexF* mtx_ii_end = element + n; element < mtx_end; element += n + 1, mtx_ii_end += n, d++)
                {
                    {

                        val = (ComplexF)MathF.Abs(*(mtx_ii_j = element));
                        for (mtx_u_ii = element + n; mtx_u_ii < mtx_end; mtx_u_ii += n)
                        {
                            if (val.Abs < (MathF.Abs(*mtx_u_ii)))
                                val = (ComplexF)MathF.Abs(*(mtx_ii_j = mtx_u_ii));
                        }

                        if (MathF.Abs(val - 0) < float.Epsilon) return (ComplexF)float.NaN;

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
                        float sizeOfThis = (float)MathF.Abs(working[x][currentRow]);
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
            public static ComplexF[][] Invert(ComplexF[][] working)
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
                ComplexF[][] inverse = Jagged.ToComplex(Jagged.Zero(dimension, dimension));
                // C# will set the initial values to zero, so to create a unit
                // matrix, I just need to fill in the diagonal elements:
                for (int loop = 0; loop < dimension; loop++) inverse[loop][loop] = new ComplexF(1.0f, 0);

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
                        ComplexF temp;
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
                        if (working[currentRow][currentRow] == new ComplexF(0, 0))
                        {
                            working[currentRow][currentRow] = new ComplexF(epsilon, 0);
                        }
                        ComplexF factor = working[lop][currentRow] / working[currentRow][currentRow];

                        // If the matrix is fairly sparse (quite common for this
                        // application), it might make sense to check that the 
                        // lower elements are not already zero before doing all
                        // the scaling and replacing:
                        if (factor != new ComplexF(0, 0))
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
                    if (working[currentRow][currentRow] == new ComplexF(0, 0))
                    {
                        working[currentRow][currentRow] = new ComplexF(epsilon, 0);
                    }

                    // Then, go up the matrix subtracting as necessary to get 
                    // rid of the remaining upper-triangular elements:
                    for (int lop = currentRow - 1; lop >= 0; lop--)
                    {
                        ComplexF factor = working[lop][currentRow] / working[currentRow][currentRow];

                        // There's only one element in working to change (the other elements
                        // in the row of working are all zero), and that will always be set
                        // to zero; but you might have to do the entire row in inverse:
                        working[lop][currentRow] = new ComplexF(0, 0);

                        if (factor != new ComplexF(0, 0))
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
                    ComplexF scale = working[loop][loop];
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
                    throw new ArgumentException(exception);

                int height = A.GetLength(0);
                int width = B[0].GetLength(0);
                int length = B.GetLength(0);
                float[][] C = Jagged.Zero(height, width);

                Parallel.For(0, height, i =>
                {
                    LinealgOptions.MatrixOperation.Whittle_Mul(A, B, C, length, width, i);
                });

                return C;
            }
            /// <summary>
            /// Implements the multiplication of matrices presented in the form of jagged arrays.
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="B">Jagged array</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[][] Mul(ComplexF[][] A, ComplexF[][] B)
            {
                if (A[0].GetLength(0) != B.GetLength(0))
                    throw new ArgumentException(exception);

                int height = A.GetLength(0);
                int width = B[0].GetLength(0);
                int length = B.GetLength(0);
                ComplexF[][] C = Jagged.Zero(height, width).ToComplex();

                Parallel.For(0, height, i =>
                {
                    LinealgOptions.MatrixOperation.Whittle_Mul(A, B, C, length, width, i);
                });

                return C;
            }
            /// <summary>
            /// Implements the multiplication of matrices presented in the form of jagged arrays.
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="B">Jagged array</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[][] Mul(ComplexF[][] A, float[][] B)
            {
                if (A[0].GetLength(0) != B.GetLength(0))
                    throw new ArgumentException(exception);

                int height = A.GetLength(0);
                int width = B[0].GetLength(0);
                int length = B.GetLength(0);
                ComplexF[][] C = Jagged.Zero(height, width).ToComplex();

                Parallel.For(0, height, i =>
                {
                    LinealgOptions.MatrixOperation.Whittle_Mul(A, B, C, length, width, i);
                });

                return C;
            }
            /// <summary>
            /// Implements the multiplication of matrices presented in the form of jagged arrays.
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="B">Jagged array</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[][] Mul(float[][] A, ComplexF[][] B)
            {
                if (A[0].GetLength(0) != B.GetLength(0))
                    throw new ArgumentException(exception);

                int height = A.GetLength(0);
                int width = B[0].GetLength(0);
                int length = B.GetLength(0);
                ComplexF[][] C = Jagged.Zero(height, width).ToComplex();

                Parallel.For(0, height, i =>
                {
                    LinealgOptions.MatrixOperation.Whittle_Mul(A, B, C, length, width, i);
                });

                return C;
            }

            #region Modified Whittle matrix multiplication
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
            private static void Whittle_Mul(ComplexF[][] A, ComplexF[][] B, ComplexF[][] C, int length, int width, int i)
            {
                ComplexF[] iRowA = A[i];
                ComplexF[] iRowC = C[i];
                int k, j;

                for (k = 0; k < length; k++)
                {
                    ComplexF[] kRowB = B[k];
                    ComplexF ikA = iRowA[k];

                    for (j = 0; j < width; j++)
                    {
                        iRowC[j] += ikA * kRowB[j];
                    }
                }
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
            private static void Whittle_Mul(ComplexF[][] A, float[][] B, ComplexF[][] C, int length, int width, int i)
            {
                ComplexF[] iRowA = A[i];
                ComplexF[] iRowC = C[i];
                int k, j;

                for (k = 0; k < length; k++)
                {
                    float[] kRowB = B[k];
                    ComplexF ikA = iRowA[k];

                    for (j = 0; j < width; j++)
                    {
                        iRowC[j] += ikA * kRowB[j];
                    }
                }
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
            private static void Whittle_Mul(float[][] A, ComplexF[][] B, ComplexF[][] C, int length, int width, int i)
            {
                float[] iRowA = A[i];
                ComplexF[] iRowC = C[i];
                int k, j;

                for (k = 0; k < length; k++)
                {
                    ComplexF[] kRowB = B[k];
                    ComplexF ikA = iRowA[k];

                    for (j = 0; j < width; j++)
                    {
                        iRowC[j] += ikA * kRowB[j];
                    }
                }
            }
            #endregion

            #endregion

            #region Copy
            /// <summary>
            /// Copies matrix.
            /// </summary>
            /// <param name="src">Source</param>
            /// <param name="dst">Destination</param>
            /// <param name="r0">R0</param>
            /// <param name="c0">C0</param>
            public static void Copy(float[,] src, float[,] dst, int r0, int c0)
            {
                int rows = src.GetLength(0), cols = src.GetLength(1);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        dst[r0 + i, c0 + j] = src[i, j];
            }
            /// <summary>
            /// Copies matrix.
            /// </summary>
            /// <param name="src">Source</param>
            /// <param name="dst">Destination</param>
            /// <param name="r0">R0</param>
            /// <param name="c0">C0</param>
            public static void Copy(ComplexF[,] src, ComplexF[,] dst, int r0, int c0)
            {
                int rows = src.GetLength(0), cols = src.GetLength(1);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        dst[r0 + i, c0 + j] = src[i, j];
            }
            /// <summary>
            /// Copies matrix.
            /// </summary>
            /// <param name="src">Source</param>
            /// <param name="dst">Destination</param>
            /// <param name="r0">R0</param>
            /// <param name="c0">C0</param>
            public static void Copy(float[,] src, ComplexF[,] dst, int r0, int c0)
            {
                int rows = src.GetLength(0), cols = src.GetLength(1);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        dst[r0 + i, c0 + j] = new ComplexF(src[i, j], 0f);
            }

            #endregion
        }
        #endregion

        #region Convolution
        /// <summary>
        /// Defines a convolution filter class.
        /// </summary>
        public static class ConvolutionFilter
        {
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

                                    k = B[i, j];

                                    if (k != 0)
                                    {
                                        s += A[ir, jr] * k;
                                        div += k;
                                    }
                                }
                            }

                            if (div != 0)
                            {
                                s /= div;
                            }

                            H[y, x] = s;
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

                                    k = B[i, j];

                                    if (k != 0)
                                    {
                                        s += A[ir, jr] * k;
                                    }
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] Conv(ComplexF[,] A, ComplexF[,] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0), r1 = B.GetLength(1);
                int r0p = r0 / 2, r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != new ComplexF(0, 0))
                                    {
                                        s += A[ir, jr] * k;
                                        div += k;
                                    }
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != new ComplexF(0, 0))
                                    {
                                        s += A[ir, jr] * k;
                                    }
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] Conv(float[,] A, ComplexF[,] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0), r1 = B.GetLength(1);
                int r0p = r0 / 2, r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != new ComplexF(0, 0))
                                    {
                                        s += A[ir, jr] * k;
                                        div += k;
                                    }
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != new ComplexF(0, 0))
                                    {
                                        s += A[ir, jr] * k;
                                    }
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] Conv(ComplexF[,] A, float[,] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0), r1 = B.GetLength(1);
                int r0p = r0 / 2, r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        float k;
                        ComplexF s, div;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != 0f)
                                    {
                                        s += A[ir, jr] * k;
                                        div += k;
                                    }
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        float k;
                        ComplexF s;
                        int i, j, x;
                        int xr, yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                for (j = 0; j < r1; j++)
                                {
                                    jr = xr + j;
                                    if (jr < 0) continue; if (jr >= width) break;

                                    k = B[i, j];

                                    if (k != 0f)
                                    {
                                        s += A[ir, jr] * k;
                                    }
                                }
                            }

                            H[y, x] = s;
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
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != 0)
                            {
                                s /= div;
                            }

                            H[y, x] = s;
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
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != 0)
                            {
                                s /= div;
                            }

                            H[y, x] = s;
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
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvHorizontal(float[,] A, ComplexF[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r1 = B.GetLength(0);
                int r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[j];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[j];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvVertical(float[,] A, ComplexF[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0);
                int r0p = r0 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = B[i];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = B[i];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvHorizontal(ComplexF[,] A, ComplexF[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r1 = B.GetLength(0);
                int r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[j];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = B[j];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvVertical(ComplexF[,] A, ComplexF[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0);
                int r0p = r0 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = B[i];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = B[i];

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvHorizontal(ComplexF[,] A, float[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r1 = B.GetLength(0);
                int r1p = r1 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = new ComplexF(B[j], 0);

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int j, x;
                        int xr, yr = y;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);
                            xr = x - r1p;
                            ir = yr;

                            for (j = 0; j < r1; j++)
                            {
                                jr = xr + j;
                                if (jr < 0) continue; if (jr >= width) break;

                                k = new ComplexF((float)B[j], 0);

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
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
            public static ComplexF[,] ConvVertical(ComplexF[,] A, float[] B, bool normalize = true)
            {
                int height = A.GetLength(0), width = A.GetLength(1);
                ComplexF[,] H = new ComplexF[height, width];
                int r0 = B.GetLength(0);
                int r0p = r0 / 2;

                if (normalize)
                {
                    // normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s, div;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0); div = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = new ComplexF(B[i], 0);

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                    div += k;
                                }
                            }

                            if (div != new ComplexF(0, 0))
                            {
                                s /= div;
                            }

                            H[y, x] = s;
                        }
                    });
                }
                else
                {
                    // non-normalize convolution:
                    Parallel.For(0, height, y =>
                    {
                        ComplexF k, s;
                        int i, x;
                        int yr = y - r0p;
                        int ir, jr;

                        for (x = 0; x < width; x++)
                        {
                            s = new ComplexF(0, 0);

                            for (i = 0; i < r0; i++)
                            {
                                ir = yr + i;
                                if (ir < 0) continue; if (ir >= height) break;

                                jr = x;

                                k = new ComplexF(B[i], 0);

                                if (k != new ComplexF(0, 0))
                                {
                                    s += A[ir, jr] * k;
                                }
                            }

                            H[y, x] = s;
                        }
                    });
                }

                return H;
            }
        }
        #endregion

        #region Mean (separable)

        /// <summary>
        /// Defines a fast mean filter class.
        /// </summary>
        public static class MeanFilter
        {
            /// <summary>
            ///  Implements local average of vector.
            /// </summary>
            /// <param name="v">Array</param>
            /// <param name="r">Radius</param>
            public static float[] Mean(float[] v, int r)
            {
                int l = v.Length;

                if (l < 2 || r < 2)
                    return v;

                float[] output = new float[l];
                int h = r >= l ? l - 1 : r;
                int w = r >> 1;
                int dl = l - w;
                float s = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    s += v[x];
                }

                for (x = 0; x < w; x++)
                {
                    output[x] = s / h;
                }

                for (x = w; x < dl; x++)
                {
                    s = s - v[x - w] + v[x + w];
                    output[x] = s / h;
                }

                for (x = dl; x < l; x++)
                {
                    s = s - v[x - w] + v[x];
                    output[x] = s / h;
                }

                return output;
            }
            /// <summary>
            ///  Implements local average of vector.
            /// </summary>
            /// <param name="v">Array</param>
            /// <param name="r">Radius</param>
            public static ComplexF[] Mean(ComplexF[] v, int r)
            {
                int l = v.Length;

                if (l < 2 || r < 2)
                    return v;

                ComplexF[] output = new ComplexF[l];
                int h = r >= l ? l - 1 : r;
                int w = r >> 1;
                int dl = l - w;
                ComplexF s = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    s += v[x];
                }

                for (x = 0; x < w; x++)
                {
                    output[x] = s / h;
                }

                for (x = w; x < dl; x++)
                {
                    s = s - v[x - w] + v[x + w];
                    output[x] = s / h;
                }

                for (x = dl; x < l; x++)
                {
                    s = s - v[x - w] + v[x];
                    output[x] = s / h;
                }

                return output;
            }

            /// <summary>
            /// Implements local average of matrix (horizontal).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="r1">Size</param>
            /// <returns>Jagged array</returns>
            public static float[,] MeanHorizontal(float[,] A, int r1)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (width < 2 || r1 < 2)
                    return A;

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
            /// Implements local average of matrix (vertical).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="r0">Size</param>
            /// <returns>Jagged array</returns>
            public static float[,] MeanVertical(float[,] A, int r0)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (height < 2 || r0 < 2)
                    return A;

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
            /// Implements local average of matrix (horizontal).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="r1">Size</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[,] MeanHorizontal(ComplexF[,] A, int r1)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (width < 2 || r1 < 2)
                    return A;

                ComplexF[,] H = new ComplexF[height, width];
                int h = r1 >= width ? width - 1 : r1;
                int v = h >> 1;
                int dl = width - v;

                Parallel.For(0, height, y =>
                {
                    ComplexF s = 0;
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
            /// Implements local average of matrix (vertical).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="r0">Size</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[,] MeanVertical(ComplexF[,] A, int r0)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (height < 2 || r0 < 2)
                    return A;

                ComplexF[,] H = new ComplexF[height, width];
                int h = r0 >= height ? height - 1 : r0;
                int v = h >> 1;
                int dl = height - v;

                Parallel.For(0, width, x =>
                {
                    ComplexF s = 0;
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
            /// Implements weighted local average of vector.
            /// </summary>
            /// <param name="values">Array of values</param>
            /// <param name="weights">Array of weights (same length as values)</param>
            /// <param name="r">Radius</param>
            /// <returns>Weighted blurred array</returns>
            public static float[] MeanWeighted(float[] values, float[] weights, int r)
            {
                int l = values.Length;

                if (l < 2 || r < 2)
                    return values;

                float[] output = new float[l];
                int h = r >= l ? l - 1 : r;
                int w = r >> 1;
                int dl = l - w;

                float sumVal = 0;
                float sumW = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    sumVal += values[x] * weights[x];
                    sumW += weights[x];
                }

                for (x = 0; x < w; x++)
                {
                    output[x] = sumVal / (sumW + 1e-8f);
                }

                for (x = w; x < dl; x++)
                {
                    int xAdd = x + w;
                    int xSub = x - w - 1;

                    if (xAdd < l)
                    {
                        sumVal += values[xAdd] * weights[xAdd];
                        sumW += weights[xAdd];
                    }
                    if (xSub >= 0)
                    {
                        sumVal -= values[xSub] * weights[xSub];
                        sumW -= weights[xSub];
                    }

                    output[x] = sumVal / (sumW + 1e-8f);
                }

                for (x = dl; x < l; x++)
                {
                    int xSub = x - w - 1;

                    if (xSub >= 0)
                    {
                        sumVal -= values[xSub] * weights[xSub];
                        sumW -= weights[xSub];
                    }

                    output[x] = sumVal / (sumW + 1e-8f);
                }

                return output;
            }
            /// <summary>
            /// Implements weighted local average of vector.
            /// </summary>
            /// <param name="values">Array of values</param>
            /// <param name="weights">Array of weights (same length as values)</param>
            /// <param name="r">Radius</param>
            /// <returns>Weighted blurred array</returns>
            public static ComplexF[] MeanWeighted(ComplexF[] values, ComplexF[] weights, int r)
            {
                int l = values.Length;

                if (l < 2 || r < 2)
                    return values;

                ComplexF[] output = new ComplexF[l];
                int h = r >= l ? l - 1 : r;
                int w = r >> 1;
                int dl = l - w;

                ComplexF sumVal = 0;
                ComplexF sumW = 0;
                int x;

                for (x = 0; x < h; x++)
                {
                    sumVal += values[x] * weights[x];
                    sumW += weights[x];
                }

                for (x = 0; x < w; x++)
                {
                    output[x] = sumVal / (sumW + 1e-8f);
                }

                for (x = w; x < dl; x++)
                {
                    int xAdd = x + w;
                    int xSub = x - w - 1;

                    if (xAdd < l)
                    {
                        sumVal += values[xAdd] * weights[xAdd];
                        sumW += weights[xAdd];
                    }
                    if (xSub >= 0)
                    {
                        sumVal -= values[xSub] * weights[xSub];
                        sumW -= weights[xSub];
                    }

                    output[x] = sumVal / (sumW + 1e-8f);
                }

                for (x = dl; x < l; x++)
                {
                    int xSub = x - w - 1;

                    if (xSub >= 0)
                    {
                        sumVal -= values[xSub] * weights[xSub];
                        sumW -= weights[xSub];
                    }

                    output[x] = sumVal / (sumW + 1e-8f);
                }

                return output;
            }

            /// <summary>
            /// Implements local weighted average of matrice (horizontal).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="weights">Weights</param>
            /// <param name="r1">Size</param>
            /// <returns>Jagged array</returns>
            /// <returns></returns>
            public static float[,] MeanHorizontalWeighted(float[,] A, float[,] weights, int r1)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (width < 2 || r1 < 2)
                    return A;

                float[,] result = new float[height, width];
                int h = r1 >= width ? width - 1 : r1;
                int v = h >> 1;
                int dl = width - v;

                Parallel.For(0, height, y =>
                {
                    float sumVal = 0;
                    float sumW = 0;
                    int x;

                    for (x = 0; x < h; x++)
                    {
                        sumVal += A[y, x] * weights[y, x];
                        sumW += weights[y, x];
                    }

                    for (x = 0; x < v; x++)
                    {
                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (x = v; x < dl; x++)
                    {
                        int xAdd = x + v;
                        int xSub = x - v - 1;

                        if (xAdd < width)
                        {
                            sumVal += A[y, xAdd] * weights[y, xAdd];
                            sumW += weights[y, xAdd];
                        }

                        if (xSub >= 0)
                        {
                            sumVal -= A[y, xSub] * weights[y, xSub];
                            sumW -= weights[y, xSub];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (x = dl; x < width; x++)
                    {
                        int xSub = x - v - 1;

                        if (xSub >= 0)
                        {
                            sumVal -= A[y, xSub] * weights[y, xSub];
                            sumW -= weights[y, xSub];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }
                });

                return result;
            }
            /// <summary>
            /// Implements local weighted average of matrice (vertical).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="weights">Weights</param>
            /// <param name="r0">Size</param>
            /// <returns>Jagged array</returns>
            public static float[,] MeanVerticalWeighted(float[,] A, float[,] weights, int r0)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (height < 2 || r0 < 2)
                    return A;

                float[,] result = new float[height, width];
                int h = r0 >= height ? height - 1 : r0;
                int v = h >> 1;
                int dl = height - v;

                Parallel.For(0, width, x =>
                {
                    float sumVal = 0;
                    float sumW = 0;
                    int y;

                    for (y = 0; y < h; y++)
                    {
                        sumVal += A[y, x] * weights[y, x];
                        sumW += weights[y, x];
                    }

                    for (y = 0; y < v; y++)
                    {
                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (y = v; y < dl; y++)
                    {
                        int yAdd = y + v;
                        int ySub = y - v - 1;

                        if (yAdd < height)
                        {
                            sumVal += A[yAdd, x] * weights[yAdd, x];
                            sumW += weights[yAdd, x];
                        }

                        if (ySub >= 0)
                        {
                            sumVal -= A[ySub, x] * weights[ySub, x];
                            sumW -= weights[ySub, x];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (y = dl; y < height; y++)
                    {
                        int ySub = y - v - 1;

                        if (ySub >= 0)
                        {
                            sumVal -= A[ySub, x] * weights[ySub, x];
                            sumW -= weights[ySub, x];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }
                });

                return result;
            }
            /// <summary>
            /// Implements local weighted average of matrice (horizontal).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="weights">Weights</param>
            /// <param name="r1">Size</param>
            /// <returns>Jagged array</returns>
            /// <returns></returns>
            public static ComplexF[,] MeanHorizontalWeighted(ComplexF[,] A, ComplexF[,] weights, int r1)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (width < 2 || r1 < 2)
                    return A;

                ComplexF[,] result = new ComplexF[height, width];
                int h = r1 >= width ? width - 1 : r1;
                int v = h >> 1;
                int dl = width - v;

                Parallel.For(0, height, y =>
                {
                    ComplexF sumVal = 0;
                    ComplexF sumW = 0;
                    int x;

                    for (x = 0; x < h; x++)
                    {
                        sumVal += A[y, x] * weights[y, x];
                        sumW += weights[y, x];
                    }

                    for (x = 0; x < v; x++)
                    {
                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (x = v; x < dl; x++)
                    {
                        int xAdd = x + v;
                        int xSub = x - v - 1;

                        if (xAdd < width)
                        {
                            sumVal += A[y, xAdd] * weights[y, xAdd];
                            sumW += weights[y, xAdd];
                        }

                        if (xSub >= 0)
                        {
                            sumVal -= A[y, xSub] * weights[y, xSub];
                            sumW -= weights[y, xSub];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (x = dl; x < width; x++)
                    {
                        int xSub = x - v - 1;

                        if (xSub >= 0)
                        {
                            sumVal -= A[y, xSub] * weights[y, xSub];
                            sumW -= weights[y, xSub];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }
                });

                return result;
            }
            /// <summary>
            /// Implements local weighted average of matrice (vertical).
            /// </summary>
            /// <param name="A">Jagged array</param>
            /// <param name="weights">Weights</param>
            /// <param name="r0">Size</param>
            /// <returns>Jagged array</returns>
            public static ComplexF[,] MeanVerticalWeighted(ComplexF[,] A, ComplexF[,] weights, int r0)
            {
                int height = A.GetLength(0), width = A.GetLength(1);

                if (height < 2 || r0 < 2)
                    return A;

                ComplexF[,] result = new ComplexF[height, width];
                int h = r0 >= height ? height - 1 : r0;
                int v = h >> 1;
                int dl = height - v;

                Parallel.For(0, width, x =>
                {
                    ComplexF sumVal = 0;
                    ComplexF sumW = 0;
                    int y;

                    for (y = 0; y < h; y++)
                    {
                        sumVal += A[y, x] * weights[y, x];
                        sumW += weights[y, x];
                    }

                    for (y = 0; y < v; y++)
                    {
                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (y = v; y < dl; y++)
                    {
                        int yAdd = y + v;
                        int ySub = y - v - 1;

                        if (yAdd < height)
                        {
                            sumVal += A[yAdd, x] * weights[yAdd, x];
                            sumW += weights[yAdd, x];
                        }

                        if (ySub >= 0)
                        {
                            sumVal -= A[ySub, x] * weights[ySub, x];
                            sumW -= weights[ySub, x];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }

                    for (y = dl; y < height; y++)
                    {
                        int ySub = y - v - 1;

                        if (ySub >= 0)
                        {
                            sumVal -= A[ySub, x] * weights[ySub, x];
                            sumW -= weights[ySub, x];
                        }

                        result[y, x] = sumVal / (sumW + 1e-8f);
                    }
                });

                return result;
            }
        }
        #endregion

        #region Morphology

        /// <summary>
        /// Defines a morphology filter class for fast sort implementation.
        /// </summary>
        public static class MorphologySortFilter
        {
            /// <summary>
            /// Applies morphology filter to 2D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r0">Radius</param>
            /// <param name="r1">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static float[,] Apply(float[,] data, int r0, int r1, MorphologyMode mode = MorphologyMode.Median)
            {
                int height = data.GetLength(0);
                int width = data.GetLength(1);

                int hW = 2 * r0 + 1;
                int wW = 2 * r1 + 1;
                int K = hW * wW;

                int rank1 = GetFilterRank(mode, K);
                int rankIndex = Math.Max(0, Math.Min(K - 1, rank1 - 1));

                var result = new float[height, width];

                Parallel.For(0, height, i =>
                {
                    float[] win = new float[K];

                    int t = 0;
                    for (int di = -r0; di <= r0; di++)
                    {
                        int ni = MathF.Range(i + di, 0, height - 1);
                        for (int dj = -r1; dj <= r1; dj++)
                        {
                            int nj = MathF.Range(0 + dj, 0, width - 1);
                            win[t++] = data[ni, nj];
                        }
                    }

                    for (int k = 1; k < K; k++)
                        FastSort(ref win, k);

                    result[i, 0] = win[rankIndex];

                    for (int j = 1; j < width; j++)
                    {
                        int outCol = MathF.Range(j - r1 - 1, 0, width - 1);
                        int inCol = MathF.Range(j + r1, 0, width - 1);

                        for (int di = -r0; di <= r0; di++)
                        {
                            int ni = MathF.Range(i + di, 0, height - 1);

                            float outVal = data[ni, outCol];
                            float inVal = data[ni, inCol];

                            RemoveOneFromSorted(win, K, outVal);

                            win[K - 1] = inVal;
                            FastSort(ref win, K - 1);
                        }

                        result[i, j] = win[rankIndex];
                    }
                });

                return result;
            }
            /// <summary>
            /// Applies morphology filter to 1D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static float[] Apply(float[] data, int r, MorphologyMode mode = MorphologyMode.Median)
            {
                if (data == null) throw new ArgumentNullException(nameof(data));
                int n = data.Length;
                if (n == 0) return Array.Empty<float>();
                if (r < 0) throw new ArgumentOutOfRangeException(nameof(r));

                int K = 2 * r + 1;
                int rank1 = GetFilterRank(mode, K);
                int rankIndex = Math.Max(0, Math.Min(K - 1, rank1 - 1));

                var result = new float[n];
                var win = new float[K];

                int t = 0;
                for (int off = -r; off <= r; off++)
                {
                    int idx = MathF.Range(0 + off, 0, n - 1);
                    win[t++] = data[idx];
                }

                for (int k = 1; k < K; k++)
                    FastSort(ref win, k);

                result[0] = win[rankIndex];

                for (int i = 1; i < n; i++)
                {
                    int outIdx = MathF.Range(i - r - 1, 0, n - 1);
                    int inIdx = MathF.Range(i + r, 0, n - 1);

                    float outVal = data[outIdx];
                    float inVal = data[inIdx];

                    RemoveOneFromSorted(win, K, outVal);

                    win[K - 1] = inVal;
                    FastSort(ref win, K - 1);

                    result[i] = win[rankIndex];
                }

                return result;
            }
            /// <summary>
            /// Removes one occurrence of x from a sorted array a of length K, shifting the elements left. 
            /// After the call, the first K-1 elements remain sorted; the last element can be overwritten.
            /// </summary>
            /// <param name="a">Array</param>
            /// <param name="K">Dimension</param>
            /// <param name="x">Value</param>
            public static void RemoveOneFromSorted(float[] a, int K, float x)
            {
                int lo = 0, hi = K - 1, pos = K;

                while (lo <= hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (a[mid] >= x) { pos = mid; hi = mid - 1; }
                    else { lo = mid + 1; }
                }

                int idx = pos;

                if (idx >= K || a[idx] != x)
                {
                    idx = -1;
                    for (int t = 0; t < K; t++)
                    {
                        if (a[t] == x) { idx = t; break; }
                    }
                    if (idx == -1) idx = Math.Min(pos, K - 1);
                }

                for (int t = idx; t < K - 1; t++)
                    a[t] = a[t + 1];
            }
            /// <summary>
            /// O(N) sort algorithm.
            /// </summary>
            /// <param name="s">Array</param>
            /// <param name="index">Index</param>
            public static void FastSort(ref float[] s, int index)
            {
                int length = s.Length - 1;

                for (int i = index; i < length; i++)
                {
                    if (s[i] > s[i + 1])
                    {
                        var t = s[i + 1];
                        s[i + 1] = s[i];
                        s[i] = t;
                    }
                    else
                        break;
                }

                for (int i = index; i > 0; i--)
                {
                    if (s[i] < s[i - 1])
                    {
                        var t = s[i - 1];
                        s[i - 1] = s[i];
                        s[i] = t;
                    }
                    else
                        break;
                }
            }
            /// <summary>
            /// Gets filter rank.
            /// </summary>
            /// <param name="mode">Mode</param>
            /// <param name="windowSize">Window size</param>
            /// <returns>Value</returns>
            public static int GetFilterRank(MorphologyMode mode, int windowSize)
            {
                return mode switch
                {
                    MorphologyMode.Erosion => 0,
                    MorphologyMode.Median => windowSize / 2,
                    _ => windowSize - 1
                };
            }
        }

        /// <summary>
        /// Defines a morphology filter class for basic implementation.
        /// </summary>
        public static class MorphologyFilter
        {
            /// <summary>
            /// Applies morphology filter to 2D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r0">Radius</param>
            /// <param name="r1">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static float[,] Apply(float[,] data, int r0, int r1, MorphologyMode mode = MorphologyMode.Median)
            {
                int height = data.GetLength(0);
                int width = data.GetLength(1);
                int windowHeight = 2 * r0 + 1;
                int windowWidth = 2 * r1 + 1;
                int windowSize = windowHeight * windowWidth;
                int rank = GetFilterRank(mode, windowSize);
                var comparer = Comparer<float>.Default;

                float[,] result = new float[height, width];

                Parallel.For(0, height, i =>
                {
                    for (int j = 0; j < width; j++)
                    {
                        var set = new HeapSet<float>(comparer);

                        for (int di = -r0; di <= r0; di++)
                        {
                            int ni = MathF.Range(i + di, 0, height - 1);
                            for (int dj = -r1; dj <= r1; dj++)
                            {
                                int nj = MathF.Range(j + dj, 0, width - 1);
                                set.Add(data[ni, nj]);
                            }
                        }

                        set.Balance(rank + 1);

                        result[i, j] = set.GetRank();
                    }
                });

                return result;
            }
            /// <summary>
            /// Applies morphology filter to 1D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static float[] Apply(float[] data, int r, MorphologyMode mode = MorphologyMode.Median)
            {
                int N = data.Length;
                int windowSize = 2 * r + 1;
                int rank = GetFilterRank(mode, windowSize);
                var comparer = Comparer<float>.Default;

                float[] result = new float[N];

                for (int i = 0; i < N; i++)
                {
                    var set = new HeapSet<float>(comparer);

                    for (int offset = -r; offset <= r; offset++)
                    {
                        int idx = MathF.Range(i + offset, 0, N - 1);
                        set.Add(data[idx]);
                    }

                    set.Balance(rank + 1);

                    result[i] = set.GetRank();
                }

                return result;
            }
            /// <summary>
            /// Gets filter rank.
            /// </summary>
            /// <param name="mode">Mode</param>
            /// <param name="windowSize">Window size</param>
            /// <returns>Value</returns>
            public static int GetFilterRank(MorphologyMode mode, int windowSize)
            {
                return mode switch
                {
                    MorphologyMode.Erosion => 0,
                    MorphologyMode.Median => windowSize / 2,
                    _ => windowSize - 1
                };
            }
        }

        /// <summary>
        /// Defines a morphology filter class for histogram implementation.
        /// </summary>
        public static class MorphologyHistogramFilter
        {
            /// <summary>
            /// Applies morphology filter to 2D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r0">Radius</param>
            /// <param name="r1">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static byte[,] Apply(byte[,] data, int r0, int r1, MorphologyMode mode = MorphologyMode.Median)
            {
                int height = data.GetLength(0);
                int width = data.GetLength(1);
                int windowHeight = 2 * r0 + 1;
                int windowWidth = 2 * r1 + 1;
                int windowSize = windowHeight * windowWidth;
                int medianPos = GetFilterRank(mode, windowSize);
                int range = byte.MaxValue + 1;

                byte[,] result = new byte[height, width];

                Parallel.For(0, height, i =>
                {
                    for (int j = 0; j < width; j++)
                    {
                        int[] histogram = new int[range];

                        for (int di = -r0; di <= r0; di++)
                        {
                            int ni = MathF.Range(i + di, 0, height - 1);

                            for (int dj = -r1; dj <= r1; dj++)
                            {
                                int nj = MathF.Range(j + dj, 0, width - 1);
                                byte val = data[ni, nj];
                                histogram[val]++;
                            }
                        }

                        int count = 0;
                        byte median = 0;

                        for (int k = 0; k < range; k++)
                        {
                            count += histogram[k];

                            if (count >= medianPos)
                            {
                                median = (byte)k;
                                break;
                            }
                        }

                        result[i, j] = median;
                    }
                });

                return result;
            }
            /// <summary>
            /// Applies morphology filter to 1D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static byte[] Apply(byte[] data, int r, MorphologyMode mode = MorphologyMode.Median)
            {
                int length = data.Length;
                int windowSize = 2 * r + 1;
                int medianPos = GetFilterRank(mode, windowSize);
                int range = byte.MaxValue + 1;

                byte[] result = new byte[length];

                for (int i = 0; i < length; i++)
                {
                    int[] histogram = new int[range];

                    for (int offset = -r; offset <= r; offset++)
                    {
                        int idx = MathF.Range(i + offset, 0, length - 1);
                        histogram[data[idx]]++;
                    }

                    int count = 0;
                    byte median = 0;

                    for (int k = 0; k < range; k++)
                    {
                        count += histogram[k];

                        if (count >= medianPos)
                        {
                            median = (byte)k;
                            break;
                        }
                    }

                    result[i] = median;
                }

                return result;
            }
            /// <summary>
            /// Gets filter rank.
            /// </summary>
            /// <param name="mode">Mode</param>
            /// <param name="windowSize">Window size</param>
            /// <returns>Value</returns>
            public static int GetFilterRank(MorphologyMode mode, int windowSize)
            {
                return mode switch
                {
                    MorphologyMode.Erosion => 1,
                    MorphologyMode.Median => (windowSize + 1) / 2,
                    _ => windowSize,
                };
            }
        }

        /// <summary>
        /// Defines a morphology filter class for fast histogram implementation.
        /// </summary>
        public static class MorphologyHistogramFastFilter
        {
            /// <summary>
            /// Applies morphology filter to 2D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r0">Radius</param>
            /// <param name="r1">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static byte[,] Apply(byte[,] data, int r0, int r1, MorphologyMode mode = MorphologyMode.Median)
            {
                int height = data.GetLength(0);
                int width = data.GetLength(1);
                int windowHeight = 2 * r0 + 1;
                int windowWidth = 2 * r1 + 1;
                int windowSize = windowHeight * windowWidth;
                int rank = GetFilterRank(mode, windowSize);
                int range = byte.MaxValue + 1;

                byte[,] result = new byte[height, width];

                Parallel.For(0, height, i =>
                {
                    int[] histogram = new int[range];

                    for (int di = -r0; di <= r0; di++)
                    {
                        int ni = MathF.Range(i + di, 0, height - 1);

                        for (int dj = -r1; dj <= r1; dj++)
                        {
                            int nj = MathF.Range(0 + dj, 0, width - 1);
                            histogram[data[ni, nj]]++;
                        }
                    }

                    result[i, 0] = GetHistogramRank(histogram, rank);

                    for (int j = 1; j < width; j++)
                    {
                        int outCol = MathF.Range(j - r1 - 1, 0, width - 1);
                        int inCol = MathF.Range(j + r1, 0, width - 1);

                        for (int di = -r0; di <= r0; di++)
                        {
                            int ni = MathF.Range(i + di, 0, height - 1);

                            byte outVal = data[ni, outCol];
                            histogram[outVal]--;

                            byte inVal = data[ni, inCol];
                            histogram[inVal]++;
                        }

                        result[i, j] = GetHistogramRank(histogram, rank);
                    }
                });

                return result;
            }
            /// <summary>
            /// Applies morphology filter to 1D array.
            /// </summary>
            /// <param name="data">Array</param>
            /// <param name="r">Radius</param>
            /// <param name="mode">Mode</param>
            /// <returns>Array</returns>
            public static byte[] Apply(byte[] data, int r, MorphologyMode mode = MorphologyMode.Median)
            {
                int length = data.Length;
                int windowSize = 2 * r + 1;
                int rank = GetFilterRank(mode, windowSize);
                int range = byte.MaxValue + 1;

                byte[] result = new byte[length];

                int[] histogram = new int[range];

                for (int i = 0; i < windowSize && i < length; i++)
                {
                    int idx = MathF.Range(i, 0, length - 1);
                    histogram[data[idx]]++;
                }

                result[0] = GetHistogramRank(histogram, rank);

                for (int i = 1; i < length; i++)
                {
                    int outIdx = MathF.Range(i - r - 1, 0, length - 1);
                    int inIdx = MathF.Range(i + r, 0, length - 1);

                    histogram[data[outIdx]]--;
                    histogram[data[inIdx]]++;

                    result[i] = GetHistogramRank(histogram, rank);
                }

                return result;
            }
            /// <summary>
            /// Gets filter rank.
            /// </summary>
            /// <param name="mode">Mode</param>
            /// <param name="windowSize">Window size</param>
            /// <returns>Value</returns>
            public static int GetFilterRank(MorphologyMode mode, int windowSize)
            {
                return mode switch
                {
                    MorphologyMode.Erosion => 1,
                    MorphologyMode.Median => (windowSize + 1) / 2,
                    _ => windowSize,
                };
            }
            /// <summary>
            /// Gets histogram rank.
            /// </summary>
            /// <param name="histogram">Histogram</param>
            /// <param name="rank">Rank</param>
            /// <returns>Value</returns>
            public static byte GetHistogramRank(int[] histogram, int rank)
            {
                int count = 0;

                for (int i = 0; i < histogram.Length; i++)
                {
                    count += histogram[i];

                    if (count >= rank)
                    {
                        return (byte)i;
                    }
                }
                return 0; // fallback (should not happen if rank valid)
            }
        }

        #endregion
    }
}
