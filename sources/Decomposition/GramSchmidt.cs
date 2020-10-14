using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines the Gram-Schmidt orthogonalization process.
    /// <remarks>
    /// In mathematics, in particular linear algebra and numerical analysis, the Gram-Schmidt process is a method of orthonormalizing a set of vectors
    /// in the space of internal works. This procedure is actively used for orthogonalization of bases.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    /// </remarks>
    /// </summary>
    [Serializable]
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
        /// Initializes the Gram-Schmidt orthogonalization process.
        /// </summary>
        /// <param name="A">Square matrix</param>
        public GramSchmidt(double[,] A)
        {
            if (!Matrice.IsSquare(A))
                throw new Exception("The matrix must be square");

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
        /// Gets the orthogonal matrix Q.
        /// </summary>
        public double[,] Q
        {
            get { return q; }
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
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
        /// Returns the projection of horizontal vectors.
        /// proj[e, a]' = (e * a') / (e * e') .* e
        /// </summary>
        /// <param name="e">Array</param>
        /// <param name="a">Array</param>
        /// <returns>Array</returns>
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
}
