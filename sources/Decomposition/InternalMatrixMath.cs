using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides shared double-precision work-buffer operations for matrix decompositions.</summary>
    /// <remarks>
    /// Public entry points validate input contracts; internal primitives assume compatible dimensions
    /// unless documented otherwise. Mutating operations work only on buffers owned by the current call.
    /// </remarks>
    internal static class InternalMatrixMath
    {
        #region Constants

        public const double Roundoff = 2.2204460492503131e-16;
        public const double SingleRoundoff = 1.1920928955078125e-7;

        #endregion

        #region Real work buffers and scalar arithmetic

        /// <summary>
        /// Copies a finite matrix and scales its largest magnitude to one before numerical reduction.
        /// </summary>
        /// <param name="matrix">Real rectangular matrix whose dimensions have already been validated.</param>
        /// <param name="scale">Receives the original maximum magnitude, or one for a zero matrix.</param>
        /// <returns>A scaled jagged copy; the caller must restore the scale to the resulting matrix or spectrum.</returns>
        public static double[][] ScaledCopyJagged(float[,] matrix, out double scale)
        {
            scale = 0;
            foreach (float value in matrix)
            {
                if (float.IsNaN(value) || float.IsInfinity(value))
                    throw new ArgumentException("The matrix must contain only finite values.", nameof(matrix));
                scale = Math.Max(scale, Math.Abs((double)value));
            }
            if (scale == 0) scale = 1;
            var result = CreateJagged(matrix.GetLength(0), matrix.GetLength(1));
            for (int i = 0; i < result.Length; i++)
                for (int j = 0; j < result[i].Length; j++) result[i][j] = matrix[i, j] / scale;
            return result;
        }

        /// <summary>
        /// Computes sqrt(a*a + b*b) by a ratio to avoid squaring the larger magnitude.
        /// </summary>
        /// <param name="a">First finite real component.</param>
        /// <param name="b">Second finite real component.</param>
        /// <returns>The nonnegative Euclidean length, including zero for two zero components.</returns>
        public static double Hypotenuse(double a, double b)
        {
            a = Math.Abs(a); b = Math.Abs(b);
            if (a < b) { double temporary = a; a = b; b = temporary; }
            if (a == 0) return 0;
            double ratio = b / a;
            return a * Math.Sqrt(1 + ratio * ratio);
        }

        /// <summary>
        /// Transfers an algebraic sign to a magnitude for stable Householder and QR shifts.
        /// </summary>
        /// <param name="magnitude">Real magnitude donor.</param>
        /// <param name="sign">Sign donor; zero selects the nonnegative sign.</param>
        /// <returns>The absolute magnitude with the selected sign.</returns>
        public static double CopySign(double magnitude, double sign)
        {
            return sign < 0 ? -Math.Abs(magnitude) : Math.Abs(magnitude);
        }

        /// <summary>
        /// Applies one common positive scale to a finite matrix pencil before QZ reduction.
        /// </summary>
        /// <param name="a">Square copy of A; overwritten by its scaled values.</param>
        /// <param name="b">Square copy of B of the same order; overwritten by its scaled values.</param>
        /// <returns>The common scale to restore to alpha and beta after eigenvector calculation.</returns>
        public static double ScalePair(double[][] a, double[][] b)
        {
            double scale = 0;
            for (int i = 0; i < a.Length; i++)
                for (int j = 0; j < a.Length; j++)
                {
                    if (double.IsNaN(a[i][j]) || double.IsInfinity(a[i][j]) || double.IsNaN(b[i][j]) || double.IsInfinity(b[i][j]))
                        throw new ArgumentException("The matrices must contain only finite values.");
                    scale = Math.Max(scale, Math.Max(Math.Abs(a[i][j]), Math.Abs(b[i][j])));
                }
            if (scale == 0) return 1;
            for (int i = 0; i < a.Length; i++)
                for (int j = 0; j < a.Length; j++)
                {
                    a[i][j] /= scale;
                    b[i][j] /= scale;
                }
            return scale;
        }

        /// <summary>
        /// Complex scalar division using a numerically stable branch (Smith’s method).
        /// Computes (xr + i·xi) / (yr + i·yi) and stores the real/imag parts in <paramref name="cdivr"/> / <paramref name="cdivi"/>.
        /// </summary>
        /// <param name="xr">Real part of the numerator.</param>
        /// <param name="xi">Imag part of the numerator.</param>
        /// <param name="yr">Real part of the denominator.</param>
        /// <param name="yi">Imag part of the denominator.</param>
        /// <param name="cdivr">[out] Real part of the quotient.</param>
        /// <param name="cdivi">[out] Imag part of the quotient.</param>
        /// <remarks>
        /// Chooses the scaling branch by comparing |yr| and |yi| to avoid overflow/underflow.
        /// If both <paramref name="yr"/> and <paramref name="yi"/> are zero, the result follows IEEE-754 (Inf/NaN).
        /// </remarks>
        public static void DivideComplex(double xr, double xi, double yr, double yi, ref double cdivr, ref double cdivi)
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

        /// <summary>Allocates independent rows for a real double-precision work matrix.</summary>
        /// <param name="rows">Nonnegative row count.</param>
        /// <param name="columns">Nonnegative column count.</param>
        /// <returns>A zero-initialized rectangular jagged matrix.</returns>
        public static double[][] CreateJagged(int rows, int columns)
        {
            var result = new double[rows][];
            for (int i = 0; i < rows; i++) result[i] = new double[columns];
            return result;
        }

        /// <summary>Creates a real identity accumulator in jagged double-precision storage.</summary>
        /// <param name="size">Nonnegative matrix order.</param>
        /// <returns>An independent identity matrix.</returns>
        public static double[][] EyeJagged(int size)
        {
            var result = CreateJagged(size, size);
            for (int i = 0; i < size; i++) result[i][i] = 1;
            return result;
        }

        /// <summary>Promotes a real rectangular input without scaling or discarding small entries.</summary>
        /// <param name="matrix">Input whose dimensions have already been validated by the caller.</param>
        /// <returns>An independent jagged double-precision copy; finite-value validation is left to the caller.</returns>
        public static double[][] CopyJagged(float[,] matrix)
        {
            var result = CreateJagged(matrix.GetLength(0), matrix.GetLength(1));
            for (int i = 0; i < result.Length; i++)
                for (int j = 0; j < result[i].Length; j++) result[i][j] = matrix[i, j];
            return result;
        }

        /// <summary>Promotes single-precision row arrays for an internal real numerical kernel.</summary>
        /// <param name="matrix">Validated nonnull row arrays.</param>
        /// <returns>An independent double-precision copy with the same row lengths.</returns>
        public static double[][] CopyJagged(float[][] matrix)
        {
            var result = new double[matrix.Length][];
            for (int i = 0; i < result.Length; i++)
            {
                result[i] = new double[matrix[i].Length];
                for (int j = 0; j < result[i].Length; j++) result[i][j] = matrix[i][j];
            }
            return result;
        }

        /// <summary>Narrows a real rectangular work buffer to the public single-precision representation.</summary>
        /// <param name="matrix">Nonempty jagged matrix with equal nonnull row lengths.</param>
        /// <returns>A new rectangular matrix; values outside the float range follow IEEE-754 conversion.</returns>
        public static float[,] Real(double[][] matrix)
        {
            var result = new float[matrix.Length, matrix[0].Length];
            for (int i = 0; i < matrix.Length; i++)
                for (int j = 0; j < matrix[i].Length; j++) result[i, j] = (float)matrix[i][j];
            return result;
        }

        #endregion

        #region Complex work buffers and conversions

        /// <summary>Copies a finite, nonempty real matrix into complex double-precision storage.</summary>
        /// <param name="a">Input matrix, which is not modified.</param>
        /// <param name="square">Whether equal dimensions are required.</param>
        /// <returns>An independent work buffer.</returns>
        public static C[,] Copy(float[,] a, bool square = false)
        {
            CheckShape(a, square);
            var b = new C[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++)
                {
                    if (float.IsNaN(a[i, j]) || float.IsInfinity(a[i, j]))
                        throw new ArgumentException("The matrix must contain finite values.", nameof(a));
                    b[i, j] = a[i, j];
                }
            return b;
        }

        /// <summary>Copies a finite, nonempty complex matrix into double-precision storage.</summary>
        /// <param name="a">Input matrix, which is not modified.</param>
        /// <param name="square">Whether equal dimensions are required.</param>
        /// <returns>An independent work buffer.</returns>
        public static C[,] Copy(Complex32[,] a, bool square = false)
        {
            CheckShape(a, square);
            var b = new C[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++)
                {
                    var z = a[i, j];
                    if (float.IsNaN(z.Real) || float.IsInfinity(z.Real) || float.IsNaN(z.Imag) || float.IsInfinity(z.Imag))
                        throw new ArgumentException("The matrix must contain finite values.", nameof(a));
                    b[i, j] = new C(z.Real, z.Imag);
                }
            return b;
        }

        /// <summary>Validates matrix dimensions before allocating numerical work buffers.</summary>
        /// <param name="a">A two-dimensional array.</param>
        /// <param name="square">Whether a square matrix is required.</param>
        public static void CheckShape(Array a, bool square = false)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (a.Rank != 2 || a.GetLength(0) == 0 || a.GetLength(1) == 0)
                throw new ArgumentException("The matrix must be nonempty and two-dimensional.", nameof(a));
            if (square && a.GetLength(0) != a.GetLength(1))
                throw new ArgumentException("The matrix must be square.", nameof(a));
        }

        /// <summary>Narrows a work matrix to complex single precision.</summary>
        /// <param name="a">Double-precision values.</param>
        /// <returns>A newly allocated complex matrix.</returns>
        public static Complex32[,] Single(C[,] a)
        {
            var b = new Complex32[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++) b[i, j] = new Complex32((float)a[i, j].Real, (float)a[i, j].Imaginary);
            return b;
        }

        /// <summary>Narrows a real-valued work matrix to single precision.</summary>
        /// <param name="a">Work buffer whose imaginary components are zero.</param>
        /// <returns>A newly allocated real matrix.</returns>
        public static float[,] Real(C[,] a)
        {
            var b = new float[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < b.GetLength(0); i++)
                for (int j = 0; j < b.GetLength(1); j++) b[i, j] = (float)a[i, j].Real;
            return b;
        }

        /// <summary>Narrows complex vector components to single precision.</summary>
        /// <param name="values">Non-null double-precision vector.</param>
        /// <returns>An independent complex single-precision vector.</returns>
        public static Complex32[] Single(C[] values)
        {
            var result = new Complex32[values.Length];
            for (int i = 0; i < values.Length; i++)
                result[i] = new Complex32((float)values[i].Real, (float)values[i].Imaginary);
            return result;
        }

        /// <summary>Narrows a real vector to single precision without changing its ordering.</summary>
        /// <param name="values">Non-null double-precision vector.</param>
        /// <returns>An independent real single-precision vector.</returns>
        public static float[] Single(double[] values)
        {
            var result = new float[values.Length];
            for (int i = 0; i < values.Length; i++) result[i] = (float)values[i];
            return result;
        }

        /// <summary>Extracts and narrows the real components of a work vector.</summary>
        /// <param name="values">Non-null work vector whose imaginary components are zero.</param>
        /// <returns>An independent real single-precision vector.</returns>
        public static float[] Real(C[] values)
        {
            var result = new float[values.Length];
            for (int i = 0; i < values.Length; i++) result[i] = (float)values[i].Real;
            return result;
        }

        #endregion

        #region Matrix and vector operations

        /// <summary>
        /// Swaps two columns in a jagged matrix <paramref name="M"/> (double[rows][cols]).
        /// No operation is performed if <paramref name="c1"/> equals <paramref name="c2"/>.
        /// </summary>
        /// <param name="M">Matrix represented as an array of row arrays (double[rows][cols]).</param>
        /// <param name="c1">Index of the first column.</param>
        /// <param name="c2">Index of the second column.</param>
        public static void SwapColumns(double[][] M, int c1, int c2)
        {
            if (c1 == c2) return;
            int rows = M.Length;

            for (int r = 0; r < rows; r++)
            {
                var tmp = M[r][c1];
                M[r][c1] = M[r][c2];
                M[r][c2] = tmp;
            }
        }

        /// <summary>Creates an identity matrix for accumulating unitary transformations.</summary>
        /// <param name="n">Nonnegative order.</param>
        /// <returns>The identity of order n.</returns>
        public static C[,] Eye(int n)
        {
            var a = new C[n, n];
            for (int i = 0; i < n; i++) a[i, i] = C.One;
            return a;
        }

        /// <summary>Computes conjugate transposition, including for real-valued work buffers.</summary>
        /// <param name="a">Input matrix.</param>
        /// <returns>The conjugate transpose.</returns>
        public static C[,] Adjoint(C[,] a)
        {
            var b = new C[a.GetLength(1), a.GetLength(0)];
            for (int i = 0; i < a.GetLength(0); i++)
                for (int j = 0; j < a.GetLength(1); j++) b[j, i] = C.Conjugate(a[i, j]);
            return b;
        }

        /// <summary>Multiplies compatible work matrices with complex double accumulation.</summary>
        /// <param name="a">Left matrix.</param>
        /// <param name="b">Right matrix.</param>
        /// <returns>The matrix product.</returns>
        public static C[,] Multiply(C[,] a, C[,] b)
        {
            if (a.GetLength(1) != b.GetLength(0)) throw new ArgumentException("Incompatible matrix dimensions.");
            var c = new C[a.GetLength(0), b.GetLength(1)];
            for (int i = 0; i < c.GetLength(0); i++)
                for (int k = 0; k < a.GetLength(1); k++)
                    for (int j = 0; j < c.GetLength(1); j++) c[i, j] += a[i, k] * b[k, j];
            return c;
        }

        /// <summary>Multiplies a work matrix by a compatible column vector.</summary>
        /// <param name="a">Matrix with as many columns as there are vector entries.</param>
        /// <param name="v">Input vector, which is not modified.</param>
        /// <returns>A new vector accumulated in row-major order with complex double arithmetic.</returns>
        public static C[] Multiply(C[,] a, C[] v)
        {
            var result = new C[a.GetLength(0)];
            for (int i = 0; i < result.Length; i++)
                for (int j = 0; j < v.Length; j++) result[i] += a[i, j] * v[j];
            return result;
        }

        /// <summary>Copies the tail of a column into an independent vector.</summary>
        /// <param name="a">Source matrix.</param>
        /// <param name="column">Valid source column index.</param>
        /// <param name="row">First row, from zero through the row count.</param>
        /// <returns>Entries from the selected row to the end of the column.</returns>
        public static C[] Column(C[,] a, int column, int row = 0)
        {
            var result = new C[a.GetLength(0) - row];
            for (int i = 0; i < result.Length; i++) result[i] = a[row + i, column];
            return result;
        }

        /// <summary>Copies and conjugates a row tail for a right Householder reflection.</summary>
        /// <param name="a">Source matrix.</param>
        /// <param name="row">Valid source row index.</param>
        /// <param name="column">First column, from zero through the column count.</param>
        /// <returns>The conjugated row tail as an independent column vector.</returns>
        public static C[] ConjugateRow(C[,] a, int row, int column = 0)
        {
            var result = new C[a.GetLength(1) - column];
            for (int j = 0; j < result.Length; j++) result[j] = C.Conjugate(a[row, column + j]);
            return result;
        }

        /// <summary>Copies a vector into a contiguous part of a matrix column.</summary>
        /// <param name="a">Target matrix, updated in place.</param>
        /// <param name="column">Valid target column index.</param>
        /// <param name="values">Source entries; their length must fit in the target column.</param>
        /// <param name="row">First target row.</param>
        public static void SetColumn(C[,] a, int column, C[] values, int row = 0)
        {
            for (int i = 0; i < values.Length; i++) a[row + i, column] = values[i];
        }

        /// <summary>Swaps two complex columns without changing any other entries.</summary>
        /// <param name="a">Matrix modified in place.</param>
        /// <param name="first">First valid column index.</param>
        /// <param name="second">Second valid column index.</param>
        public static void SwapColumns(C[,] a, int first, int second)
        {
            for (int i = 0; i < a.GetLength(0); i++)
            {
                C value = a[i, first]; a[i, first] = a[i, second]; a[i, second] = value;
            }
        }

        /// <summary>Divides a finite matrix by a nonzero real scale without forming its reciprocal.</summary>
        /// <param name="a">Matrix modified in place.</param>
        /// <param name="scale">Finite nonzero divisor.</param>
        public static void Divide(C[,] a, double scale)
        {
            for (int i = 0; i < a.GetLength(0); i++)
                for (int j = 0; j < a.GetLength(1); j++) a[i, j] /= scale;
        }

        /// <summary>Divides a finite vector by a nonzero real scale without forming its reciprocal.</summary>
        /// <param name="v">Vector modified in place.</param>
        /// <param name="scale">Finite nonzero divisor, normally a vector norm.</param>
        public static void Divide(C[] v, double scale)
        {
            for (int i = 0; i < v.Length; i++) v[i] /= scale;
        }

        /// <summary>Computes a scaled Euclidean norm without squaring large or tiny entries directly.</summary>
        /// <param name="v">Vector with finite entries.</param>
        /// <returns>The nonnegative Euclidean norm.</returns>
        public static double Norm(C[] v)
        {
            double scale = 0, sum = 1;
            foreach (C value in v)
                AccumulateNorm(C.Abs(value), ref scale, ref sum);
            return scale == 0 ? 0 : scale * Math.Sqrt(sum);
        }

        /// <summary>Computes a stable column norm without allocating a temporary vector.</summary>
        /// <param name="a">Finite work matrix.</param>
        /// <param name="column">Valid column index.</param>
        /// <returns>The nonnegative Euclidean column norm, including zero for a zero column.</returns>
        public static double ColumnNorm(C[,] a, int column)
        {
            double scale = 0, sum = 1;
            for (int i = 0; i < a.GetLength(0); i++)
                AccumulateNorm(C.Abs(a[i, column]), ref scale, ref sum);
            return scale == 0 ? 0 : scale * Math.Sqrt(sum);
        }

        /// <summary>Accumulates a sum of squared magnitudes using the largest magnitude as a scale.</summary>
        /// <param name="magnitude">Finite nonnegative component magnitude.</param>
        /// <param name="scale">Largest magnitude seen so far, initialized to zero.</param>
        /// <param name="sum">Scaled sum of squares, initialized to one.</param>
        /// <remarks>Rescaling the old sum avoids squaring very large or very small unscaled components.</remarks>
        private static void AccumulateNorm(double magnitude, ref double scale, ref double sum)
        {
            if (magnitude == 0) return;
            if (scale < magnitude)
            {
                double ratio = scale / magnitude;
                sum = 1 + sum * ratio * ratio;
                scale = magnitude;
            }
            else
            {
                double ratio = magnitude / scale;
                sum += ratio * ratio;
            }
        }

        /// <summary>Returns the largest entry magnitude, or zero for a zero matrix.</summary>
        /// <param name="a">Finite work matrix.</param>
        /// <returns>A nonnegative scale.</returns>
        public static double Max(C[,] a)
        {
            double scale = 0;
            foreach (C z in a) scale = Math.Max(scale, C.Abs(z));
            return scale;
        }

        /// <summary>Checks the Hermitian condition with a relative single-precision tolerance.</summary>
        /// <param name="a">Square input matrix.</param>
        public static void RequireHermitian(C[,] a)
        {
            if (!IsHermitian(a)) throw new ArgumentException("The matrix must be Hermitian (symmetric for real inputs).");
        }

        /// <summary>Detects Hermitian structure with an optional relative tolerance.</summary>
        /// <param name="a">Finite square work matrix.</param>
        /// <param name="relativeTolerance">Nonnegative fraction of the matrix scale; zero requires exact conjugate symmetry.</param>
        /// <returns>True when conjugate symmetry holds within the tolerance, defaulting to eight single-precision rounding units.</returns>
        public static bool IsHermitian(C[,] a, double relativeTolerance = 8 * SingleRoundoff)
        {
            int n = a.GetLength(0);
            if (n != a.GetLength(1)) return false;
            double tolerance = relativeTolerance * Max(a);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                    if (C.Abs(a[i, j] - C.Conjugate(a[j, i])) > tolerance)
                        return false;
            return true;
        }

        /// <summary>Copies a rectangular leading block or a block starting at specified offsets.</summary>
        /// <param name="a">Source matrix.</param>
        /// <param name="rows">Number of rows.</param>
        /// <param name="columns">Number of columns.</param>
        /// <param name="row">First source row.</param>
        /// <param name="column">First source column.</param>
        /// <returns>An independent block.</returns>
        public static C[,] Block(C[,] a, int rows, int columns, int row = 0, int column = 0)
        {
            var b = new C[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) b[i, j] = a[i + row, j + column];
            return b;
        }

        #endregion

        #region Unitary transformations

        /// <summary>Builds a unit Householder vector mapping x onto its first coordinate.</summary>
        /// <param name="x">Finite vector modified in place; a zero vector is returned unchanged.</param>
        /// <returns>The same array, containing a normalized vector v with H = I - 2 v v^H, or zero for identity.</returns>
        /// <remarks>The target is -phase(x[0])*norm(x), with phase(0)=1, to avoid cancellation.</remarks>
        public static C[] HouseholderVector(C[] x)
        {
            double norm = Norm(x);
            if (norm == 0) return x;
            C phase = C.Abs(x[0]) == 0 ? C.One : x[0] / C.Abs(x[0]);
            Divide(x, norm);
            x[0] += phase;
            norm = Norm(x);
            Divide(x, norm);
            return x;
        }

        /// <summary>Applies I - 2 v v^H to selected rows from the left, in place.</summary>
        /// <param name="a">Work matrix to update.</param>
        /// <param name="v">Normalized reflection vector, or zero for identity.</param>
        /// <param name="row">First affected row.</param>
        /// <param name="column">First affected column.</param>
        public static void ReflectLeft(C[,] a, C[] v, int row, int column)
        {
            for (int j = column; j < a.GetLength(1); j++)
            {
                C dot = 0;
                for (int i = 0; i < v.Length; i++) dot += C.Conjugate(v[i]) * a[row + i, j];
                for (int i = 0; i < v.Length; i++) a[row + i, j] -= 2 * v[i] * dot;
            }
        }

        /// <summary>Applies I - 2 v v^H to selected columns from the right, in place.</summary>
        /// <param name="a">Work matrix to update.</param>
        /// <param name="v">Normalized reflection vector, or zero for identity.</param>
        /// <param name="column">First affected column.</param>
        /// <param name="row">First affected row.</param>
        public static void ReflectRight(C[,] a, C[] v, int column, int row)
        {
            for (int i = row; i < a.GetLength(0); i++)
            {
                C dot = 0;
                for (int j = 0; j < v.Length; j++) dot += a[i, column + j] * v[j];
                for (int j = 0; j < v.Length; j++) a[i, column + j] -= 2 * dot * C.Conjugate(v[j]);
            }
        }

        /// <summary>Constructs a complex Givens rotation annihilating the second component.</summary>
        /// <param name="f">First component.</param>
        /// <param name="g">Second component.</param>
        /// <returns>Real cosine C and complex sine S defining [C,S;-conj(S),C].</returns>
        /// <remarks>The rotation maps (f,g) to (r,0); two zero inputs select identity, and f=0 selects a positive real r.</remarks>
        public static (double C, C S) Givens(C f, C g)
        {
            double af = C.Abs(f), ag = C.Abs(g);
            if (ag == 0) return (1, C.Zero);
            if (af == 0) return (0, C.Conjugate(g) / ag);
            double scale = Math.Max(af, ag);
            double norm = scale * Math.Sqrt((af / scale) * (af / scale) + (ag / scale) * (ag / scale));
            return (af / norm, (f / af) * (C.Conjugate(g) / norm));
        }

        /// <summary>Applies the transpose of G = [c,s;-conj(s),c] to an ordered pair of columns from the right.</summary>
        /// <param name="a">Target matrix.</param>
        /// <param name="i">First column.</param>
        /// <param name="j">Second column.</param>
        /// <param name="c">Real cosine with c squared plus the squared magnitude of s equal to one.</param>
        /// <param name="s">Complex sine; conjugate it to accumulate the adjoint of a left G rotation.</param>
        public static void RotateColumns(C[,] a, int i, int j, double c, C s)
        {
            for (int k = 0; k < a.GetLength(0); k++)
            {
                C x = a[k, i], y = a[k, j];
                a[k, i] = c * x + s * y;
                a[k, j] = -C.Conjugate(s) * x + c * y;
            }
        }

        /// <summary>Applies G = [c,s;-conj(s),c] to an ordered pair of rows from the left.</summary>
        /// <param name="a">Target matrix, modified in place.</param>
        /// <param name="i">First valid row index.</param>
        /// <param name="j">Second valid row index, distinct from the first.</param>
        /// <param name="c">Real cosine with c squared plus the squared magnitude of s equal to one.</param>
        /// <param name="s">Complex sine following the convention returned by Givens.</param>
        public static void RotateRows(C[,] a, int i, int j, double c, C s)
        {
            for (int k = 0; k < a.GetLength(1); k++)
            {
                C x = a[i, k], y = a[j, k];
                a[i, k] = c * x + s * y;
                a[j, k] = -C.Conjugate(s) * x + c * y;
            }
        }

        /// <summary>Applies A = U^H A U and Q = Q U for an embedded unitary transformation U.</summary>
        /// <param name="a">Full square work matrix.</param>
        /// <param name="q">Distinct square accumulator with the same dimensions as a.</param>
        /// <param name="rotation">Unitary transformation on a contiguous active block.</param>
        /// <param name="offset">First index of the active block.</param>
        public static void ApplySimilarity(C[,] a, C[,] q, C[,] rotation, int offset)
        {
            int n = a.GetLength(0), size = rotation.GetLength(0);
            var buffer = new C[size];
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < size; i++)
                {
                    buffer[i] = 0;
                    for (int k = 0; k < size; k++) buffer[i] += C.Conjugate(rotation[k, i]) * a[offset + k, j];
                }
                for (int i = 0; i < size; i++) a[offset + i, j] = buffer[i];
            }
            foreach (var target in new[] { a, q })
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < size; j++)
                    {
                        buffer[j] = 0;
                        for (int k = 0; k < size; k++) buffer[j] += target[i, offset + k] * rotation[k, j];
                    }
                    for (int j = 0; j < size; j++) target[i, offset + j] = buffer[j];
                }
        }

        #endregion

        #region Orthogonalization

        /// <summary>Completes an orthonormal basis after exact or numerical breakdown.</summary>
        /// <param name="q">Matrix containing orthonormal columns before column k.</param>
        /// <param name="k">Column to complete.</param>
        /// <returns>A unit vector orthogonal to the preceding columns.</returns>
        public static C[] Complete(C[,] q, int k)
        {
            int n = q.GetLength(0);
            var best = new C[n];
            double largest = -1;
            for (int seed = 0; seed < n; seed++)
            {
                var v = new C[n]; v[seed] = 1;
                Orthogonalize(v, q, k);
                double norm = Norm(v);
                if (norm > largest) { largest = norm; best = v; }
            }
            if (largest <= Roundoff) throw new InvalidOperationException("Unable to complete an orthonormal basis.");
            Divide(best, largest);
            return best;
        }

        /// <summary>Applies modified Gram-Schmidt using the Hermitian inner product and optional coefficient accumulation.</summary>
        /// <param name="v">Vector modified in place.</param>
        /// <param name="q">Previously computed orthonormal columns.</param>
        /// <param name="columns">Number of columns to remove.</param>
        /// <param name="passes">Positive pass count; two reorthogonalizes against roundoff.</param>
        /// <param name="coefficients">Optional matrix receiving each projection, added to its existing entries.</param>
        /// <param name="column">Target column in the coefficient matrix, ignored when it is null.</param>
        public static void Orthogonalize(C[] v, C[,] q, int columns, int passes = 2, C[,] coefficients = null, int column = 0)
        {
            for (int pass = 0; pass < passes; pass++)
                for (int j = 0; j < columns; j++)
                {
                    C dot = 0;
                    for (int i = 0; i < v.Length; i++) dot += C.Conjugate(q[i, j]) * v[i];
                    if (coefficients != null) coefficients[j, column] += dot;
                    for (int i = 0; i < v.Length; i++) v[i] -= q[i, j] * dot;
                }
        }

        #endregion
    }
}
