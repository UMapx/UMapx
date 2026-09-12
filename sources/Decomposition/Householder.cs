using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides Householder reflections and symmetric or Hermitian tridiagonal reduction.</summary>
    public static class Householder
    {
        /// <summary>Reduces a symmetric matrix as A = H T H^T.</summary>
        /// <param name="matrix">Finite nonempty symmetric square matrix.</param>
        /// <returns>Orthogonal H and symmetric tridiagonal T.</returns>
        public static (float[,] H, float[,] T) Decompose(float[,] matrix)
        {
            var d = Tridiagonalize(MatrixMath.Copy(matrix, true));
            return (MatrixMath.Real(d.P), MatrixMath.Real(d.H));
        }

        /// <summary>Reduces a Hermitian matrix as A = H T H^H.</summary>
        /// <param name="matrix">Finite nonempty Hermitian square matrix.</param>
        /// <returns>Unitary H and Hermitian tridiagonal T.</returns>
        public static (Complex32[,] H, Complex32[,] T) Decompose(Complex32[,] matrix)
        {
            var d = Tridiagonalize(MatrixMath.Copy(matrix, true));
            return (MatrixMath.Single(d.P), MatrixMath.Single(d.H));
        }

        /// <summary>Constructs a reflection that maps a vector onto its first coordinate.</summary>
        /// <param name="vector">Finite nonempty vector to reduce.</param>
        /// <returns>The orthogonal reflection, or identity for a zero vector.</returns>
        public static float[,] Reflection(float[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            var a = new float[vector.Length, 1];
            for (int i = 0; i < vector.Length; i++) a[i, 0] = vector[i];
            return MatrixMath.Real(Reflect(MatrixMath.Copy(a)));
        }

        /// <summary>Constructs a reflection mapping x to -phase(x[0])*norm(x) times the first coordinate vector.</summary>
        /// <param name="vector">Finite nonempty complex vector to reduce; phase(0) is defined as one.</param>
        /// <returns>The Hermitian unitary reflection, or identity for a zero vector.</returns>
        public static Complex32[,] Reflection(Complex32[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            var a = new Complex32[vector.Length, 1];
            for (int i = 0; i < vector.Length; i++) a[i, 0] = vector[i];
            return MatrixMath.Single(Reflect(MatrixMath.Copy(a)));
        }

        /// <summary>Forms the reflection that annihilates the tail of a column vector.</summary>
        /// <param name="a">Validated single-column work matrix.</param>
        /// <returns>The reflection, with a zero vector interpreted as identity.</returns>
        private static C[,] Reflect(C[,] a)
        {
            int n = a.GetLength(0);
            var v = new C[n];
            for (int i = 0; i < n; i++) v[i] = a[i, 0];
            v = Vector(v);
            var h = MatrixMath.Eye(n);
            ApplyLeft(h, v, 0, 0);
            return h;
        }

        /// <summary>Checks Hermitian structure and removes roundoff outside the tridiagonal band.</summary>
        /// <param name="a">Private square input buffer.</param>
        /// <returns>The similarity transformation and tridiagonal matrix.</returns>
        private static (C[,] P, C[,] H) Tridiagonalize(C[,] a)
        {
            MatrixMath.RequireHermitian(a);
            var d = Hessenberg.Factor(a);
            int n = a.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                d.H[i, i] = d.H[i, i].Real;
                for (int j = i + 1; j < n; j++)
                    d.H[i, j] = j == i + 1 ? C.Conjugate(d.H[j, i]) : C.Zero;
            }
            return d;
        }

        /// <summary>Builds a unit Householder vector mapping x onto its first coordinate.</summary>
        /// <param name="x">Finite vector; zero is returned unchanged.</param>
        /// <returns>A normalized vector v with H = I - 2 v v^H.</returns>
        /// <remarks>The target is -phase(x[0])*norm(x), with phase(0)=1, to avoid cancellation.</remarks>
        internal static C[] Vector(C[] x)
        {
            double norm = MatrixMath.Norm(x);
            if (norm == 0) return x;
            C phase = C.Abs(x[0]) == 0 ? C.One : x[0] / C.Abs(x[0]);
            for (int i = 0; i < x.Length; i++) x[i] /= norm;
            x[0] += phase;
            norm = MatrixMath.Norm(x);
            for (int i = 0; i < x.Length; i++) x[i] /= norm;
            return x;
        }

        /// <summary>Applies I - 2 v v^H to selected rows from the left, in place.</summary>
        /// <param name="a">Work matrix to update.</param>
        /// <param name="v">Normalized reflection vector, or zero for identity.</param>
        /// <param name="row">First affected row.</param>
        /// <param name="column">First affected column.</param>
        internal static void ApplyLeft(C[,] a, C[] v, int row, int column)
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
        internal static void ApplyRight(C[,] a, C[] v, int column, int row)
        {
            for (int i = row; i < a.GetLength(0); i++)
            {
                C dot = 0;
                for (int j = 0; j < v.Length; j++) dot += a[i, column + j] * v[j];
                for (int j = 0; j < v.Length; j++) a[i, column + j] -= 2 * dot * C.Conjugate(v[j]);
            }
        }
    }
}
