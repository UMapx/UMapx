using UMapx.Core;

namespace UMapx.Decomposition
{
    public partial class GEVD
    {
        /// <summary>
        /// Performs the QZ reduction of matrices A and B.
        /// </summary>
        /// <param name="a">Matrix A (will be overwritten by the quasi-triangular form S).</param>
        /// <param name="b">Matrix B (will be overwritten by the upper triangular form T).</param>
        /// <param name="eps">Epsilon [0, 1]</param>
        /// <param name="z">Matrix that accumulates the right orthogonal transformations.</param>
        /// <param name="ierr">Convergence flag.</param>
        internal static void DecomposeQZ(float[][] a, float[][] b, float eps, float[][] z, ref int ierr)
        {
            int n = a.Length;
            bool matz = true;
            qzhes(n, a, b, matz, z);
            qzit(n, a, b, eps, matz, z, ref ierr);
        }
    }
}
