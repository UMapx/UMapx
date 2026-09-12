using System;
using UMapx.Core;
using C = System.Numerics.Complex;

namespace UMapx.Decomposition
{
    /// <summary>Provides Householder reflections and symmetric or Hermitian tridiagonal reduction</summary>
    public static class Householder
    {
        /// <summary>Reduces a symmetric matrix as A = H T H^T</summary>
        /// <param name="matrix">Finite nonempty symmetric square matrix.</param>
        /// <returns>Orthogonal H and symmetric tridiagonal T.</returns>
        public static (float[,] H, float[,] T) Decompose(float[,] matrix)
        {
            var d = Tridiagonalize(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Real(d.P), InternalMatrixMath.Real(d.H));
        }

        /// <summary>Reduces a Hermitian matrix as A = H T H^H</summary>
        /// <param name="matrix">Finite nonempty Hermitian square matrix.</param>
        /// <returns>Unitary H and Hermitian tridiagonal T.</returns>
        public static (Complex32[,] H, Complex32[,] T) Decompose(Complex32[,] matrix)
        {
            var d = Tridiagonalize(InternalMatrixMath.Copy(matrix, true));
            return (InternalMatrixMath.Single(d.P), InternalMatrixMath.Single(d.H));
        }

        /// <summary>Constructs a reflection that maps a vector onto its first coordinate</summary>
        /// <param name="vector">Finite nonempty vector to reduce.</param>
        /// <returns>The orthogonal reflection, or identity for a zero vector.</returns>
        public static float[,] Reflection(float[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            var a = new float[vector.Length, 1];
            for (int i = 0; i < vector.Length; i++) a[i, 0] = vector[i];
            return InternalMatrixMath.Real(Reflect(InternalMatrixMath.Copy(a)));
        }

        /// <summary>Constructs a reflection mapping x to -phase(x[0])*norm(x) times the first coordinate vector</summary>
        /// <param name="vector">Finite nonempty complex vector to reduce; phase(0) is defined as one.</param>
        /// <returns>The Hermitian unitary reflection, or identity for a zero vector.</returns>
        public static Complex32[,] Reflection(Complex32[] vector)
        {
            if (vector == null) throw new ArgumentNullException(nameof(vector));
            var a = new Complex32[vector.Length, 1];
            for (int i = 0; i < vector.Length; i++) a[i, 0] = vector[i];
            return InternalMatrixMath.Single(Reflect(InternalMatrixMath.Copy(a)));
        }

        /// <summary>Forms the reflection that annihilates the tail of a column vector</summary>
        /// <param name="a">Validated single-column work matrix.</param>
        /// <returns>The reflection, with a zero vector interpreted as identity.</returns>
        private static C[,] Reflect(C[,] a)
        {
            int n = a.GetLength(0);
            var v = InternalMatrixMath.Column(a, 0);
            v = InternalMatrixMath.HouseholderVector(v);
            var h = InternalMatrixMath.Eye(n);
            InternalMatrixMath.ReflectLeft(h, v, 0, 0);
            return h;
        }

        /// <summary>Checks Hermitian structure and removes roundoff outside the tridiagonal band</summary>
        /// <param name="a">Private square input buffer.</param>
        /// <returns>The similarity transformation and tridiagonal matrix.</returns>
        private static (C[,] P, C[,] H) Tridiagonalize(C[,] a)
        {
            InternalMatrixMath.RequireHermitian(a);
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

    }
}
