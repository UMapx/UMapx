using System;
using System.Numerics;

namespace UMapx.Decomposition
{
    /// <summary>Real double-precision operations on private, row-contiguous decomposition buffers.</summary>
    /// <remarks>Only Copy and Validate accept public inputs. Other operations assume compatible private buffers.</remarks>
    internal static class InternalRealMatrixMath
    {
        public const double Roundoff = InternalMatrixMath.Roundoff;
        public const double SingleRoundoff = InternalMatrixMath.SingleRoundoff;

        /// <summary>Checks finite entries and dimensions without allocating a work matrix.</summary>
        public static void Validate(float[,] a, bool square = false)
        {
            InternalMatrixMath.CheckShape(a, square);
            foreach (float value in a)
                if (float.IsNaN(value) || float.IsInfinity(value))
                    throw new ArgumentException("The matrix must contain finite values.", nameof(a));
        }

        /// <summary>Validates and promotes a real input, without changing its scale.</summary>
        public static double[][] Copy(float[,] a, bool square = false)
        {
            Validate(a, square);
            return InternalMatrixMath.CopyJagged(a);
        }

        public static double[][] Create(int rows, int columns) => InternalMatrixMath.CreateJagged(rows, columns);
        public static double[][] Eye(int n) => InternalMatrixMath.EyeJagged(n);
        public static float[,] Real(double[][] a) => InternalMatrixMath.Real(a);
        public static float[] Real(double[] a) => InternalMatrixMath.Single(a);

        /// <summary>Narrows a leading block directly, avoiding an intermediate block copy.</summary>
        public static float[,] Real(double[][] a, int rows, int columns)
        {
            var result = new float[rows, columns];
            for (int i = 0; i < rows; i++) for (int j = 0; j < columns; j++) result[i, j] = (float)a[i][j];
            return result;
        }

        public static double[][] Transpose(double[][] a)
        {
            var result = Create(a[0].Length, a.Length);
            for (int i = 0; i < a.Length; i++) for (int j = 0; j < a[i].Length; j++) result[j][i] = a[i][j];
            return result;
        }

        public static float[,] RealTranspose(double[][] a)
        {
            var result = new float[a[0].Length, a.Length];
            for (int i = 0; i < a.Length; i++) for (int j = 0; j < a[i].Length; j++) result[j, i] = (float)a[i][j];
            return result;
        }

        public static double[][] Block(double[][] a, int rows, int columns, int row = 0, int column = 0)
        {
            var result = Create(rows, columns);
            for (int i = 0; i < rows; i++) Array.Copy(a[row + i], column, result[i], 0, columns);
            return result;
        }

        public static double[] Column(double[][] a, int column, int row = 0)
        {
            var result = new double[a.Length - row];
            for (int i = 0; i < result.Length; i++) result[i] = a[row + i][column];
            return result;
        }

        public static double[] Row(double[][] a, int row, int column = 0)
        {
            var result = new double[a[row].Length - column];
            Array.Copy(a[row], column, result, 0, result.Length);
            return result;
        }

        public static void SetColumn(double[][] a, int column, double[] values, int row = 0)
        {
            for (int i = 0; i < values.Length; i++) a[row + i][column] = values[i];
        }

        public static double Max(double[][] a)
        {
            double scale = 0;
            foreach (var row in a) foreach (double value in row) scale = Math.Max(scale, Math.Abs(value));
            return scale;
        }

        public static void Divide(double[][] a, double scale)
        {
            foreach (var row in a) Divide(row, scale);
        }

        public static void Divide(double[] a, double scale)
        {
            for (int i = 0; i < a.Length; i++) a[i] /= scale;
        }

        /// <summary>Computes a scaled norm, preserving tiny values beside very large entries.</summary>
        public static double Norm(double[] a)
        {
            double scale = 0;
            foreach (double value in a) scale = Math.Max(scale, Math.Abs(value));
            if (scale == 0) return 0;
            double sum = 0;
            foreach (double value in a) { double v = value / scale; sum += v * v; }
            return scale * Math.Sqrt(sum);
        }

        public static double ColumnNorm(double[][] a, int column)
        {
            double scale = 0;
            for (int i = 0; i < a.Length; i++) scale = Math.Max(scale, Math.Abs(a[i][column]));
            if (scale == 0) return 0;
            double sum = 0;
            for (int i = 0; i < a.Length; i++) { double value = a[i][column] / scale; sum += value * value; }
            return scale * Math.Sqrt(sum);
        }

        /// <summary>Accumulates a real inner product with SIMD, retaining double precision.</summary>
        public static double Dot(double[] a, double[] b, int length, int offsetA = 0, int offsetB = 0)
        {
            int i = 0; double sum = 0;
            if (Vector.IsHardwareAccelerated && length >= Vector<double>.Count)
            {
                var total = Vector<double>.Zero;
                for (; i <= length - Vector<double>.Count; i += Vector<double>.Count)
                    total += new Vector<double>(a, offsetA + i) * new Vector<double>(b, offsetB + i);
                sum = Vector.Dot(total, Vector<double>.One);
            }
            for (; i < length; i++) sum += a[offsetA + i] * b[offsetB + i];
            return sum;
        }

        /// <summary>Computes a diagonally weighted inner product for LDL and UDL updates.</summary>
        public static double WeightedDot(double[] a, double[] b, double[] weights, int start, int length)
        {
            int i = start, end = start + length; double sum = 0;
            if (Vector.IsHardwareAccelerated && length >= Vector<double>.Count)
            {
                var total = Vector<double>.Zero;
                for (; i <= end - Vector<double>.Count; i += Vector<double>.Count)
                    total += new Vector<double>(a, i) * new Vector<double>(weights, i) * new Vector<double>(b, i);
                sum = Vector.Dot(total, Vector<double>.One);
            }
            for (; i < end; i++) sum += a[i] * weights[i] * b[i];
            return sum;
        }

        /// <summary>Checks symmetry using the same tolerance as the complex Hermitian path.</summary>
        public static void RequireSymmetric(double[][] a)
        {
            double tolerance = 8 * SingleRoundoff * Max(a);
            for (int i = 0; i < a.Length; i++) for (int j = 0; j < i; j++)
                if (Math.Abs(a[i][j] - a[j][i]) > tolerance)
                    throw new ArgumentException("The matrix must be symmetric.");
        }

        /// <summary>Normalizes a Householder vector in place; a zero vector denotes identity.</summary>
        public static double[] HouseholderVector(double[] x)
        {
            double norm = Norm(x);
            if (norm == 0) return x;
            Divide(x, norm);
            x[0] += x[0] < 0 ? -1 : 1;
            Divide(x, Norm(x));
            return x;
        }

        /// <summary>Applies I-2vv^T from the left, traversing contiguous rows in both passes.</summary>
        public static void ReflectLeft(double[][] a, double[] v, int row, int column, double[] scratch = null)
        {
            int columns = a[0].Length;
            var dots = scratch ?? new double[columns];
            Array.Clear(dots, column, columns - column);
            for (int i = 0; i < v.Length; i++)
            {
                double value = v[i]; var source = a[row + i];
                for (int j = column; j < columns; j++) dots[j] += value * source[j];
            }
            for (int i = 0; i < v.Length; i++)
            {
                double value = 2 * v[i]; var target = a[row + i];
                for (int j = column; j < columns; j++) target[j] -= value * dots[j];
            }
        }

        public static void ReflectRight(double[][] a, double[] v, int column, int row)
        {
            for (int i = row; i < a.Length; i++)
            {
                var target = a[i]; double dot = 2 * Dot(target, v, v.Length, column);
                int j = 0;
                if (Vector.IsHardwareAccelerated)
                {
                    var factor = new Vector<double>(dot);
                    for (; j <= v.Length - Vector<double>.Count; j += Vector<double>.Count)
                        (new Vector<double>(target, column + j) - factor * new Vector<double>(v, j)).CopyTo(target, column + j);
                }
                for (; j < v.Length; j++) target[column + j] -= dot * v[j];
            }
        }

        /// <summary>Accumulates a two-row QZ reflection in contiguous rows of the transposed left factor.</summary>
        public static void ReflectRows(double[][] a, int first, int second, double u2, double v1, double v2)
        {
            var x = a[first]; var y = a[second]; int j = 0;
            if (Vector.IsHardwareAccelerated)
            {
                var u = new Vector<double>(u2); var p = new Vector<double>(v1); var q = new Vector<double>(v2);
                for (; j <= x.Length - Vector<double>.Count; j += Vector<double>.Count)
                {
                    var xx = new Vector<double>(x, j); var yy = new Vector<double>(y, j);
                    var t = xx + u * yy;
                    (xx + p * t).CopyTo(x, j); (yy + q * t).CopyTo(y, j);
                }
            }
            for (; j < x.Length; j++) { double t = x[j] + u2 * y[j]; x[j] += v1 * t; y[j] += v2 * t; }
        }

        /// <summary>Accumulates a three-row QZ reflection without strided column access.</summary>
        public static void ReflectRows(double[][] a, int first, int second, int third, double u2, double u3, double v1, double v2, double v3)
        {
            var x = a[first]; var y = a[second]; var z = a[third]; int j = 0;
            if (Vector.IsHardwareAccelerated)
            {
                var u = new Vector<double>(u2); var v = new Vector<double>(u3);
                var p = new Vector<double>(v1); var q = new Vector<double>(v2); var r = new Vector<double>(v3);
                for (; j <= x.Length - Vector<double>.Count; j += Vector<double>.Count)
                {
                    var xx = new Vector<double>(x, j); var yy = new Vector<double>(y, j); var zz = new Vector<double>(z, j);
                    var t = xx + u * yy + v * zz;
                    (xx + p * t).CopyTo(x, j); (yy + q * t).CopyTo(y, j); (zz + r * t).CopyTo(z, j);
                }
            }
            for (; j < x.Length; j++)
            {
                double t = x[j] + u2 * y[j] + u3 * z[j];
                x[j] += v1 * t; y[j] += v2 * t; z[j] += v3 * t;
            }
        }

        public static double[][] Multiply(double[][] a, double[][] b)
        {
            var result = Create(a.Length, b[0].Length);
            for (int i = 0; i < a.Length; i++)
            {
                var target = result[i]; var source = a[i];
                for (int k = 0; k < source.Length; k++)
                {
                    double value = source[k]; var right = b[k];
                    for (int j = 0; j < right.Length; j++) target[j] += value * right[j];
                }
            }
            return result;
        }

        public static double[] Multiply(double[][] a, double[] v)
        {
            var result = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = Dot(a[i], v, v.Length);
            }
            return result;
        }

        public static void Orthogonalize(double[] v, double[][] q, int columns, int passes = 2, double[][] coefficients = null, int column = 0)
        {
            for (int pass = 0; pass < passes; pass++)
                for (int j = 0; j < columns; j++)
                {
                    double dot = 0;
                    for (int i = 0; i < v.Length; i++) dot += q[i][j] * v[i];
                    if (coefficients != null) coefficients[j][column] += dot;
                    for (int i = 0; i < v.Length; i++) v[i] -= q[i][j] * dot;
                }
        }

        /// <summary>Completes a basis deterministically at numerical rank breakdown.</summary>
        public static double[] Complete(double[][] q, int k)
        {
            double bestNorm = -1; double[] best = null;
            for (int seed = 0; seed < q.Length; seed++)
            {
                var v = new double[q.Length]; v[seed] = 1;
                Orthogonalize(v, q, k);
                double norm = Norm(v);
                if (norm > bestNorm) { bestNorm = norm; best = v; }
                if (norm > 0.5) break;
            }
            if (!(bestNorm > 0)) throw new InvalidOperationException("Cannot complete the orthonormal basis.");
            Divide(best, bestNorm);
            return best;
        }
    }
}
