using System;
using UMapx.Core;

namespace UMapx.Decomposition
{
    /// <summary>
    /// Defines generalized eigenvalue decomposition.
    /// </summary>
    /// <remarks>
    /// It is the task of finding a vector of values of V such that the representation: A * V = B * V * D.
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem
    /// </remarks>
    [Serializable]
    public class GEVD
    {
        #region Private data
        private int n;
        private float[] ar;
        private float[] ai;
        private float[] beta;
        private float[][] Z;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes generalized eigenvalue decomposition.
        /// </summary>
        /// <param name="a">Matrix A</param>
        /// <param name="b">Matrix B</param>
        /// <param name="eps">Epsilon [0, 1]</param>
        public GEVD(float[,] a, float[,] b, float eps = 1e-16f)
        {
            if (a.GetLength(0) != a.GetLength(1))
                throw new ArgumentException("The matrix must be square");

            if (b.GetLength(0) != b.GetLength(1))
                throw new ArgumentException("The matrix must be square");

            if (a.GetLength(0) != b.GetLength(0) || a.GetLength(1) != b.GetLength(1))
                throw new ArgumentException("Matrices should be the same size");

            // params
            this.n = a.GetLength(0);
            this.ar = new float[n];
            this.ai = new float[n];
            this.beta = new float[n];
            this.Z = Jagged.Zero(n, n);
            var A = Jagged.ToJagged(a);
            var B = Jagged.ToJagged(b);
            bool matz = true;
            int ierr = 0;


            // reduces A to upper Hessenberg form and B to upper
            // triangular form using orthogonal transformations
            qzhes(n, A, B, matz, Z);

            // reduces the Hessenberg matrix A to quasi-triangular form
            // using orthogonal transformations while maintaining the
            // triangular form of the B matrix.
            qzit(n, A, B, Maths.Float(eps), matz, Z, ref ierr);

            // reduces the quasi-triangular matrix further, so that any
            // remaining 2-by-2 blocks correspond to pairs of complex
            // eigenvalues, and returns quantities whose ratios give the
            // generalized eigenvalues.
            qzval(n, A, B, ar, ai, beta, matz, Z);

            // computes the eigenvectors of the triangular problem and
            // transforms the results back to the original coordinate system.
            qzvec(n, A, B, ar, ai, beta, Z);
        }
        #endregion

        #region Standard voids
        /// <summary>
        /// Returns vector α.
        /// </summary>
        public Complex32[] Alpha
        {
            get
            {
                Complex32[] a = new Complex32[n];

                for (int i = 0; i < n; i++)
                    a[i] = new Complex32(ar[i], ai[i]);

                return a;
            }
        }
        /// <summary>
        /// Returns vector β.
        /// </summary>
        public float[] Beta
        {
            get { return beta; }
        }
        /// <summary>
        /// Returns a vector of eigenvalues.
        /// </summary>
        public Complex32[] Eigenvalues
        {
            get
            {
                // c = (ar + j*ai) / beta
                Complex32[] c = new Complex32[n];
                for (int i = 0; i < n; i++)
                    c[i] = new Complex32(ar[i], ai[i]) / beta[i];

                return c;
            }
        }
        /// <summary>
        /// Returns an eigenvalue matrix D.
        /// </summary>
        public float[,] D
        {
            get
            {
                float[,] x = new float[n, n];

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                        x[i, j] = 0.0f;

                    x[i, i] = ar[i] / beta[i];
                    if (ai[i] > 0)
                        x[i, i + 1] = ai[i] / beta[i];
                    else if (ai[i] < 0)
                        x[i, i - 1] = ai[i] / beta[i];
                }

                return x;
            }
        }
        /// <summary>
        /// Returns a matrix of values V.
        /// </summary>
        public float[,] V
        {
            get
            {
                return Jagged.FromJagged(Z);
            }
        }
        /// <summary>
        /// Checks whether one of the matrices is singular or not.
        /// </summary>
        public bool IsSingular
        {
            get
            {
                for (int i = 0; i < n; i++)
                    if (beta[i] == 0)
                        return true;
                return false;
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Performs the QZ reduction of matrices A and B.
        /// </summary>
        /// <param name="a">Matrix A (will be overwritten by the quasi-triangular form S)</param>
        /// <param name="b">Matrix B (will be overwritten by the upper triangular form T)</param>
        /// <param name="eps">Epsilon [0, 1]</param>
        /// <param name="z">Matrix that accumulates the right orthogonal transformations</param>
        /// <param name="ierr">Convergence flag</param>
        internal static void qzdecomp(float[][] a, float[][] b, float eps, float[][] z, ref int ierr)
        {
            int n = a.Length;
            bool matz = true;
            qzhes(n, a, b, matz, z);
            qzit(n, a, b, eps, matz, z, ref ierr);
        }
        /// <summary>
        /// First QZ stage: reduce B to upper triangular form and A to upper Hessenberg form using orthogonal transformations.
        /// </summary>
        /// <remarks>
        /// Port of EISPACK <c>QZHES</c> (Moler–Stewart, SIAM J. Numer. Anal., 1973).
        /// <para>On exit: B is upper triangular; A is upper Hessenberg. If <paramref name="matz"/> is true, Z accumulates the
        /// right transformations such that the similarity actions A ← GᵀAG, B ← GᵀBG are recorded in Z ← ZG.</para>
        /// <para>All operations are in place; input matrices are overwritten.</para>
        /// </remarks>
        /// <param name="n">Matrix order (rows = cols = n)</param>
        /// <param name="a">Input A; overwritten with its upper-Hessenberg form. Modified in place</param>
        /// <param name="b">Input B; overwritten with its upper-triangular form. Modified in place</param>
        /// <param name="matz">
        /// If true, right orthogonal transformations are accumulated into <paramref name="z"/>; otherwise Z is ignored
        /// </param>
        /// <param name="z">
        /// Right transformation accumulator Z (n×n). If <paramref name="matz"/> is true, updated as Z ← ZG; otherwise unused
        /// </param>
        private static void qzhes(int n, float[][] a, float[][] b, bool matz, float[][] z)
        {
            int i, j, k, l;
            float r, s, t;
            int l1;
            float u1, u2, v1, v2;
            int lb, nk1;
            float rho;

            // Reduce b to upper triangular form
            if (n <= 1) return;
            for (l = 0; l < n - 1; ++l)
            {
                l1 = l + 1;
                s = 0.0f;

                for (i = l1; i < n; ++i)
                    s += (System.Math.Abs(b[i][l]));

                if (s == 0.0) continue;
                s += (Math.Abs(b[l][l]));
                r = 0.0f;

                for (i = l; i < n; ++i)
                {
                    // Computing 2nd power
                    b[i][l] /= s;
                    r += b[i][l] * b[i][l];
                }

                r = Sign(Maths.Sqrt(r), b[l][l]);
                b[l][l] += r;
                rho = r * b[l][l];

                for (j = l1; j < n; ++j)
                {
                    t = 0.0f;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * b[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        b[i][j] += t * b[i][l];
                }

                for (j = 0; j < n; ++j)
                {
                    t = 0.0f;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * a[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        a[i][j] += t * b[i][l];
                }

                b[l][l] = -s * r;
                for (i = l1; i < n; ++i)
                    b[i][l] = 0.0f;
            }

            // Reduce a to upper Hessenberg form, while keeping b triangular
            if (n == 2) return;
            for (k = 0; k < n - 2; ++k)
            {
                nk1 = n - 2 - k;

                // for l=n-1 step -1 until k+1 do
                for (lb = 0; lb < nk1; ++lb)
                {
                    l = n - lb - 2;
                    l1 = l + 1;

                    // Zero a(l+1,k)
                    s = (Math.Abs(a[l][k])) + (Math.Abs(a[l1][k]));

                    if (s == 0.0) continue;
                    u1 = a[l][k] / s;
                    u2 = a[l1][k] / s;
                    r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (j = k; j < n; ++j)
                    {
                        t = a[l][j] + u2 * a[l1][j];
                        a[l][j] += t * v1;
                        a[l1][j] += t * v2;
                    }

                    a[l1][k] = 0.0f;

                    for (j = l; j < n; ++j)
                    {
                        t = b[l][j] + u2 * b[l1][j];
                        b[l][j] += t * v1;
                        b[l1][j] += t * v2;
                    }

                    // Zero b(l+1,l)
                    s = (System.Math.Abs(b[l1][l1])) + (System.Math.Abs(b[l1][l]));

                    if (s == 0.0) continue;
                    u1 = b[l1][l1] / s;
                    u2 = b[l1][l] / s;
                    r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (i = 0; i <= l1; ++i)
                    {
                        t = b[i][l1] + u2 * b[i][l];
                        b[i][l1] += t * v1;
                        b[i][l] += t * v2;
                    }

                    b[l1][l] = 0.0f;

                    for (i = 0; i < n; ++i)
                    {
                        t = a[i][l1] + u2 * a[i][l];
                        a[i][l1] += t * v1;
                        a[i][l] += t * v2;
                    }

                    if (matz)
                    {
                        for (i = 0; i < n; ++i)
                        {
                            t = z[i][l1] + u2 * z[i][l];
                            z[i][l1] += t * v1;
                            z[i][l] += t * v2;
                        }
                    }
                }
            }
        }
        /// <summary>
        /// Second QZ stage: iterative QZ step (bulge chasing) to drive A to quasi-triangular form while keeping B triangular.
        /// </summary>
        /// <remarks>
        /// Port of EISPACK <c>QZIT</c> with Ward's NASA TN D-7305 modification. Implements single- and double-shift
        /// iterations, choosing shifts from trailing 1×1 or 2×2 blocks. Converges subdiagonals of A to zero (up to tolerance),
        /// revealing 1×1 (real) and 2×2 (complex-conjugate) blocks.
        /// <para>Overwrites A, B, and (optionally) Z in place. <paramref name="ierr"/> is set if the iteration cap (≈ 30·n)
        /// is exceeded.</para>
        /// </remarks>
        /// <param name="n">Matrix order</param>
        /// <param name="a">On entry: upper-Hessenberg from <see cref="qzhes"/>; on exit: quasi-triangular S. Modified in place</param>
        /// <param name="b">On entry: upper-triangular from <see cref="qzhes"/>; on exit: upper-triangular T. Modified in place</param>
        /// <param name="eps1">
        /// Relative convergence tolerance. If zero, machine roundoff is used (via <see cref="Epsilon(float)"/>)
        /// </param>
        /// <param name="matz">If true, accumulate right transformations into <paramref name="z"/></param>
        /// <param name="z">Right orthogonal accumulator Z (updated if <paramref name="matz"/> is true)</param>
        /// <param name="ierr">
        /// Output status: 0 if all subdiagonals converged; otherwise set to <c>en+1</c> at failure as in EISPACK
        /// </param>
        private static void qzit(int n, float[][] a, float[][] b, float eps1, bool matz, float[][] z, ref int ierr)
        {

            int i, j, k, l = 0;
            float r, s, t, a1, a2, a3 = 0;
            int k1, k2, l1, ll;
            float u1, u2, u3;
            float v1, v2, v3;
            float a11, a12, a21, a22, a33, a34, a43, a44;
            float b11, b12, b22, b33, b34, b44;
            int na, en, ld;
            float ep;
            float sh = 0;
            int km1, lm1 = 0;
            float ani, bni;
            int ish, itn, its, enm2, lor1;
            float epsa, epsb, anorm = 0, bnorm = 0;
            int enorn;
            bool notlas;


            ierr = 0;

            #region Compute epsa and epsb
            for (i = 0; i < n; ++i)
            {
                ani = 0.0f;
                bni = 0.0f;

                if (i != 0)
                    ani = (Math.Abs(a[i][(i - 1)]));

                for (j = i; j < n; ++j)
                {
                    ani += Math.Abs(a[i][j]);
                    bni += Math.Abs(b[i][j]);
                }

                if (ani > anorm) anorm = ani;
                if (bni > bnorm) bnorm = bni;
            }

            if (anorm == 0.0) anorm = 1.0f;
            if (bnorm == 0.0) bnorm = 1.0f;

            ep = eps1;
            if (ep == 0.0)
            {
                // Use round-off level if eps1 is zero
                ep = Epsilon(1.0f);
            }

            epsa = ep * anorm;
            epsb = ep * bnorm;
            #endregion


            // Reduce a to quasi-triangular form, while keeping b triangular
            lor1 = 0;
            enorn = n;
            en = n - 1;
            itn = n * 30;

        // Begin QZ step
        L60:
            if (en <= 1) goto L1001;
            if (!matz) enorn = en + 1;

            its = 0;
            na = en - 1;
            enm2 = na;

        L70:
            ish = 2;
            // Check for convergence or reducibility.
            for (ll = 0; ll <= en; ++ll)
            {
                lm1 = en - ll - 1;
                l = lm1 + 1;

                if (l + 1 == 1)
                    goto L95;

                if ((Math.Abs(a[l][lm1])) <= epsa)
                    break;
            }

        L90:
            a[l][lm1] = 0.0f;
            if (l < na) goto L95;

            // 1-by-1 or 2-by-2 block isolated
            en = lm1;
            goto L60;

        // Check for small top of b 
        L95:
            ld = l;

        L100:
            l1 = l + 1;
            b11 = b[l][l];

            if (Math.Abs(b11) > epsb) goto L120;

            b[l][l] = 0.0f;
            s = (Math.Abs(a[l][l]) + Math.Abs(a[l1][l]));
            u1 = a[l][l] / s;
            u2 = a[l1][l] / s;
            r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
            v1 = -(u1 + r) / r;
            v2 = -u2 / r;
            u2 = v2 / v1;

            for (j = l; j < enorn; ++j)
            {
                t = a[l][j] + u2 * a[l1][j];
                a[l][j] += t * v1;
                a[l1][j] += t * v2;

                t = b[l][j] + u2 * b[l1][j];
                b[l][j] += t * v1;
                b[l1][j] += t * v2;
            }

            if (l != 0)
                a[l][lm1] = -a[l][lm1];

            lm1 = l;
            l = l1;
            goto L90;

        L120:
            a11 = a[l][l] / b11;
            a21 = a[l1][l] / b11;
            if (ish == 1) goto L140;

            // Iteration strategy
            if (itn == 0) goto L1000;
            if (its == 10) goto L155;

            // Determine type of shift
            b22 = b[l1][l1];
            if (Math.Abs(b22) < epsb) b22 = epsb;
            b33 = b[na][na];
            if (Math.Abs(b33) < epsb) b33 = epsb;
            b44 = b[en][en];
            if (Math.Abs(b44) < epsb) b44 = epsb;
            a33 = a[na][na] / b33;
            a34 = a[na][en] / b44;
            a43 = a[en][na] / b33;
            a44 = a[en][en] / b44;
            b34 = b[na][en] / b44;
            t = (a43 * b34 - a33 - a44) * .5f;
            r = t * t + a34 * a43 - a33 * a44;
            if (r < 0.0) goto L150;

            // Determine single shift zero-th column of a
            ish = 1;
            r = Maths.Sqrt(r);
            sh = -t + r;
            s = -t - r;
            if (Math.Abs(s - a44) < Math.Abs(sh - a44))
                sh = s;

            // Look for two consecutive small sub-diagonal elements of a.
            for (ll = ld; ll + 1 <= enm2; ++ll)
            {
                l = enm2 + ld - ll - 1;

                if (l == ld)
                    goto L140;

                lm1 = l - 1;
                l1 = l + 1;
                t = a[l + 1][l + 1];

                if (Math.Abs(b[l][l]) > epsb)
                    t -= sh * b[l][l];

                if (Math.Abs(a[l][lm1]) <= (Math.Abs(t / a[l1][l])) * epsa)
                    goto L100;
            }

        L140:
            a1 = a11 - sh;
            a2 = a21;
            if (l != ld)
                a[l][lm1] = -a[l][lm1];
            goto L160;

        // Determine float shift zero-th column of a
        L150:
            a12 = a[l][l1] / b22;
            a22 = a[l1][l1] / b22;
            b12 = b[l][l1] / b22;
            a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + a12 - a11 * b12;
            a2 = a22 - a11 - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
            a3 = a[l1 + 1][l1] / b22;
            goto L160;

        // Ad hoc shift
        L155:
            a1 = 0.0f;
            a2 = 1.0f;
            a3 = 1.1605f;

        L160:
            ++its;
            --itn;

            if (!matz) lor1 = ld;

            // Main loop
            for (k = l; k <= na; ++k)
            {
                notlas = k != na && ish == 2;
                k1 = k + 1;
                k2 = k + 2;

                km1 = Math.Max(k, l + 1) - 1; // Computing MAX
                ll = Math.Min(en, k1 + ish);  // Computing MIN

                if (notlas) goto L190;

                // Zero a(k+1,k-1)
                if (k == l) goto L170;
                a1 = a[k][km1];
                a2 = a[k1][km1];

            L170:
                s = Math.Abs(a1) + Math.Abs(a2);
                if (s == 0.0) goto L70;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;

                    t = b[k][j] + u2 * b[k1][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                }

                if (k != l)
                    a[k1][km1] = 0.0f;
                goto L240;

            // Zero a(k+1,k-1) and a(k+2,k-1)
            L190:
                if (k == l) goto L200;
                a1 = a[k][km1];
                a2 = a[k1][km1];
                a3 = a[k2][km1];

            L200:
                s = Math.Abs(a1) + Math.Abs(a2) + Math.Abs(a3);
                if (s == 0.0) goto L260;
                u1 = a1 / s;
                u2 = a2 / s;
                u3 = a3 / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j] + u3 * a[k2][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;
                    a[k2][j] += t * v3;

                    t = b[k][j] + u2 * b[k1][j] + u3 * b[k2][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                    b[k2][j] += t * v3;
                }

                if (k == l) goto L220;
                a[k1][km1] = 0.0f;
                a[k2][km1] = 0.0f;

            // Zero b(k+2,k+1) and b(k+2,k)
            L220:
                s = (Math.Abs(b[k2][k2])) + (Math.Abs(b[k2][k1])) + (Math.Abs(b[k2][k]));
                if (s == 0.0) goto L240;
                u1 = b[k2][k2] / s;
                u2 = b[k2][k1] / s;
                u3 = b[k2][k] / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k2] + u2 * a[i][k1] + u3 * a[i][k];
                    a[i][k2] += t * v1;
                    a[i][k1] += t * v2;
                    a[i][k] += t * v3;

                    t = b[i][k2] + u2 * b[i][k1] + u3 * b[i][k];
                    b[i][k2] += t * v1;
                    b[i][k1] += t * v2;
                    b[i][k] += t * v3;
                }

                b[k2][k] = 0.0f;
                b[k2][k1] = 0.0f;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k2] + u2 * z[i][k1] + u3 * z[i][k];
                        z[i][k2] += t * v1;
                        z[i][k1] += t * v2;
                        z[i][k] += t * v3;
                    }
                }

            // Zero b(k+1,k)
            L240:
                s = (Math.Abs(b[k1][k1])) + (Math.Abs(b[k1][k]));
                if (s == 0.0) goto L260;
                u1 = b[k1][k1] / s;
                u2 = b[k1][k] / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k1] + u2 * a[i][k];
                    a[i][k1] += t * v1;
                    a[i][k] += t * v2;

                    t = b[i][k1] + u2 * b[i][k];
                    b[i][k1] += t * v1;
                    b[i][k] += t * v2;
                }

                b[k1][k] = 0.0f;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k1] + u2 * z[i][k];
                        z[i][k1] += t * v1;
                        z[i][k] += t * v2;
                    }
                }

            L260:
                ;
            }

            goto L70; // End QZ step

        // Set error -- all eigenvalues have not converged after 30*n iterations
        L1000:
            ierr = en + 1;

        // Save epsb for use by qzval and qzvec
        L1001:
            if (n > 1)
                b[n - 1][0] = epsb;
            return;
        }
        /// <summary>
        /// Third QZ stage: extract generalized eigenvalue numerators/denominators from (S,T) after convergence.
        /// </summary>
        /// <remarks>
        /// Port of EISPACK <c>QZVAL</c>. Given S (quasi-triangular) and T (upper triangular), this routine
        /// reduces any remaining 2×2 blocks in S consistently with T, and fills arrays such that
        /// λ_j = (alfr[j] + i·alfi[j]) / beta[j].
        /// <para>
        /// If <paramref name="matz"/> is true, right transformations are accumulated into Z for later use by <see cref="qzvec"/>.
        /// </para>
        /// </remarks>
        /// <param name="n">Matrix order</param>
        /// <param name="a">On entry: S from <see cref="qzit"/>; may be locally modified. On exit: still quasi-triangular</param>
        /// <param name="b">On entry: T from <see cref="qzit"/>; may be locally modified. On exit: still upper triangular</param>
        /// <param name="alfr">Real parts of α_j (numerators) for generalized λ_j = α_j / β_j. Length n. Written by the routine</param>
        /// <param name="alfi">Imag parts of α_j. Length n. Written by the routine (zero for real eigenvalues)</param>
        /// <param name="beta">Denominators β_j (nonnegative). Length n. Written by the routine</param>
        /// <param name="matz">If true, accumulate right transformations into <paramref name="z"/></param>
        /// <param name="z">Right transformation accumulator Z (updated if <paramref name="matz"/> is true)</param>
        private static void qzval(int n, float[][] a, float[][] b, float[] alfr, float[] alfi, float[] beta, bool matz, float[][] z)
        {
            int i, j;
            int na, en, nn;
            float c, d, e = 0;
            float r, s, t;
            float a1, a2, u1, u2, v1, v2;
            float a11, a12, a21, a22;
            float b11, b12, b22;
            float di, ei;
            float an = 0, bn;
            float cq, dr;
            float cz, ti, tr;
            float a1i, a2i, a11i, a12i, a22i, a11r, a12r, a22r;
            float sqi, ssi, sqr, szi, ssr, szr;

            float epsb = b[n - 1][0];
            int isw = 1;


            // Find eigenvalues of quasi-triangular matrices.
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;

                if (isw == 2) goto L505;
                if (en == 0) goto L410;
                if (a[en][na] != 0.0) goto L420;

                // 1-by-1 block, one real root
                L410:
                alfr[en] = a[en][en];
                if (b[en][en] < 0.0)
                {
                    alfr[en] = -alfr[en];
                }
                beta[en] = (Math.Abs(b[en][en]));
                alfi[en] = 0.0f;
                goto L510;

            // 2-by-2 block
            L420:
                if (Math.Abs(b[na][na]) <= epsb) goto L455;
                if (Math.Abs(b[en][en]) > epsb) goto L430;
                a1 = a[en][en];
                a2 = a[en][na];
                bn = 0.0f;
                goto L435;

            L430:
                an = Math.Abs(a[na][na]) + Math.Abs(a[na][en]) + Math.Abs(a[en][na]) + Math.Abs(a[en][en]);
                bn = Math.Abs(b[na][na]) + Math.Abs(b[na][en]) + Math.Abs(b[en][en]);
                a11 = a[na][na] / an;
                a12 = a[na][en] / an;
                a21 = a[en][na] / an;
                a22 = a[en][en] / an;
                b11 = b[na][na] / bn;
                b12 = b[na][en] / bn;
                b22 = b[en][en] / bn;
                e = a11 / b11;
                ei = a22 / b22;
                s = a21 / (b11 * b22);
                t = (a22 - e * b22) / b22;

                if (Math.Abs(e) <= Math.Abs(ei))
                    goto L431;

                e = ei;
                t = (a11 - e * b11) / b11;

            L431:
                c = (t - s * b12) * .5f;
                d = c * c + s * (a12 - e * b12);
                if (d < 0.0) goto L480;

                // Two real roots. Zero both a(en,na) and b(en,na)
                e += c + Sign(Maths.Sqrt(d), c);
                a11 -= e * b11;
                a12 -= e * b12;
                a22 -= e * b22;

                if (Math.Abs(a11) + Math.Abs(a12) < Math.Abs(a21) + Math.Abs(a22))
                    goto L432;

                a1 = a12;
                a2 = a11;
                goto L435;

            L432:
                a1 = a22;
                a2 = a21;

            // Choose and apply real z
            L435:
                s = Math.Abs(a1) + Math.Abs(a2);
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = 0; i <= en; ++i)
                {
                    t = a[i][en] + u2 * a[i][na];
                    a[i][en] += t * v1;
                    a[i][na] += t * v2;

                    t = b[i][en] + u2 * b[i][na];
                    b[i][en] += t * v1;
                    b[i][na] += t * v2;
                }

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][en] + u2 * z[i][na];
                        z[i][en] += t * v1;
                        z[i][na] += t * v2;
                    }
                }

                if (bn == 0.0) goto L475;
                if (an < System.Math.Abs(e) * bn) goto L455;
                a1 = b[na][na];
                a2 = b[en][na];
                goto L460;

            L455:
                a1 = a[na][na];
                a2 = a[en][na];

            // Choose and apply real q
            L460:
                s = System.Math.Abs(a1) + System.Math.Abs(a2);
                if (s == 0.0) goto L475;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Sign(Maths.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = na; j < n; ++j)
                {
                    t = a[na][j] + u2 * a[en][j];
                    a[na][j] += t * v1;
                    a[en][j] += t * v2;

                    t = b[na][j] + u2 * b[en][j];
                    b[na][j] += t * v1;
                    b[en][j] += t * v2;
                }

            L475:
                a[en][na] = 0.0f;
                b[en][na] = 0.0f;
                alfr[na] = a[na][na];
                alfr[en] = a[en][en];

                if (b[na][na] < 0.0)
                    alfr[na] = -alfr[na];

                if (b[en][en] < 0.0)
                    alfr[en] = -alfr[en];

                beta[na] = (System.Math.Abs(b[na][na]));
                beta[en] = (System.Math.Abs(b[en][en]));
                alfi[en] = 0.0f;
                alfi[na] = 0.0f;
                goto L505;

            // Two complex roots
            L480:
                e += c;
                ei = Maths.Sqrt(-d);
                a11r = a11 - e * b11;
                a11i = ei * b11;
                a12r = a12 - e * b12;
                a12i = ei * b12;
                a22r = a22 - e * b22;
                a22i = ei * b22;

                if (System.Math.Abs(a11r) + System.Math.Abs(a11i) +
                    System.Math.Abs(a12r) + System.Math.Abs(a12i) <
                    System.Math.Abs(a21) + System.Math.Abs(a22r)
                    + System.Math.Abs(a22i))
                    goto L482;

                a1 = a12r;
                a1i = a12i;
                a2 = -a11r;
                a2i = -a11i;
                goto L485;

            L482:
                a1 = a22r;
                a1i = a22i;
                a2 = -a21;
                a2i = 0.0f;

            // Choose complex z
            L485:
                cz = Maths.Sqrt(a1 * a1 + a1i * a1i);
                if (cz == 0.0) goto L487;
                szr = (a1 * a2 + a1i * a2i) / cz;
                szi = (a1 * a2i - a1i * a2) / cz;
                r = Maths.Sqrt(cz * cz + szr * szr + szi * szi);
                cz /= r;
                szr /= r;
                szi /= r;
                goto L490;

            L487:
                szr = 1.0f;
                szi = 0.0f;

            L490:
                if (an < (System.Math.Abs(e) + ei) * bn) goto L492;
                a1 = cz * b11 + szr * b12;
                a1i = szi * b12;
                a2 = szr * b22;
                a2i = szi * b22;
                goto L495;

            L492:
                a1 = cz * a11 + szr * a12;
                a1i = szi * a12;
                a2 = cz * a21 + szr * a22;
                a2i = szi * a22;

            // Choose complex q
            L495:
                cq = Maths.Sqrt(a1 * a1 + a1i * a1i);
                if (cq == 0.0) goto L497;
                sqr = (a1 * a2 + a1i * a2i) / cq;
                sqi = (a1 * a2i - a1i * a2) / cq;
                r = Maths.Sqrt(cq * cq + sqr * sqr + sqi * sqi);
                cq /= r;
                sqr /= r;
                sqi /= r;
                goto L500;

            L497:
                sqr = 1.0f;
                sqi = 0.0f;

            // Compute diagonal elements that would result if transformations were applied
            L500:
                ssr = sqr * szr + sqi * szi;
                ssi = sqr * szi - sqi * szr;
                i = 0;
                tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
                ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
                dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
                di = cq * szi * b12 + ssi * b22;
                goto L503;

            L502:
                i = 1;
                tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
                ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
                dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
                di = -ssi * b11 - sqi * cz * b12;

            L503:
                t = ti * dr - tr * di;
                j = na;

                if (t < 0.0)
                    j = en;

                r = Maths.Sqrt(dr * dr + di * di);
                beta[j] = bn * r;
                alfr[j] = an * (tr * dr + ti * di) / r;
                alfi[j] = an * t / r;
                if (i == 0) goto L502;

                L505:
                isw = 3 - isw;

            L510:
                ;
            }

            b[n - 1][0] = epsb;

            return;
        }
        /// <summary>
        /// Optional fourth QZ stage: compute (right) generalized eigenvectors and map them back via Z.
        /// </summary>
        /// <remarks>
        /// Port of EISPACK <c>QZVEC</c>. Solves the triangular generalized eigenvalue problems for (S,T),
        /// performing back-substitution to obtain eigenvectors in the deflated coordinate system and then
        /// transforms them to the original coordinates via Z.
        /// <para>
        /// On entry, (alfr, alfi, beta) must come from <see cref="qzval"/>. On exit, Z contains eigenvectors
        /// (real columns for real λ; 2 columns for each complex-conjugate pair).
        /// </para>
        /// </remarks>
        /// <param name="n">Matrix order</param>
        /// <param name="a">Quasi-triangular S (read-only here). Notation follows EISPACK; may be referenced</param>
        /// <param name="b">Upper triangular T (read-only here). Notation follows EISPACK; may be referenced</param>
        /// <param name="alfr">Real parts of α_j from <see cref="qzval"/></param>
        /// <param name="alfi">Imag parts of α_j from <see cref="qzval"/></param>
        /// <param name="beta">β_j from <see cref="qzval"/></param>
        /// <param name="z">
        /// On entry: right transformation accumulator from prior stages. On exit: columns replaced by the (right)
        /// generalized eigenvectors in the original coordinate system. Modified in place
        /// </param>
        private static void qzvec(int n, float[][] a, float[][] b, float[] alfr, float[] alfi, float[] beta, float[][] z)
        {
            int i, j, k, m;
            int na, ii, en, jj, nn, enm2;
            float d, q;
            float r = 0, s = 0, t, w, x = 0, y, t1, t2, w1, x1 = 0, z1 = 0, di;
            float ra, dr, sa;
            float ti, rr, tr, zz = 0;
            float alfm, almi, betm, almr;

            float epsb = b[n - 1][0];
            int isw = 1;


            // for en=n step -1 until 1 do --
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;
                if (isw == 2) goto L795;
                if (alfi[en] != 0.0) goto L710;

                // Real vector
                m = en;
                b[en][en] = 1.0f;
                if (na == -1) goto L800;
                alfm = alfr[m];
                betm = beta[m];

                // for i=en-1 step -1 until 1 do --
                for (ii = 0; ii <= na; ++ii)
                {
                    i = en - ii - 1;
                    w = betm * a[i][i] - alfm * b[i][i];
                    r = 0.0f;

                    for (j = m; j <= en; ++j)
                        r += (betm * a[i][j] - alfm * b[i][j]) * b[j][en];

                    if (i == 0 || isw == 2)
                        goto L630;

                    if (betm * a[i][i - 1] == 0.0)
                        goto L630;

                    zz = w;
                    s = r;
                    goto L690;

                L630:
                    m = i;
                    if (isw == 2) goto L640;

                    // Real 1-by-1 block
                    t = w;
                    if (w == 0.0)
                        t = epsb;
                    b[i][en] = -r / t;
                    goto L700;

                // Real 2-by-2 block
                L640:
                    x = betm * a[i][i + 1] - alfm * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    q = w * zz - x * y;
                    t = (x * s - zz * r) / q;
                    b[i][en] = t;
                    if (Math.Abs(x) <= Math.Abs(zz)) goto L650;
                    b[i + 1][en] = (-r - w * t) / x;
                    goto L690;

                L650:
                    b[i + 1][en] = (-s - y * t) / zz;

                L690:
                    isw = 3 - isw;

                L700:
                    ;
                }
                // End real vector
                goto L800;

            // Complex vector
            L710:
                m = na;
                almr = alfr[m];
                almi = alfi[m];
                betm = beta[m];

                // last vector component chosen imaginary so that eigenvector matrix is triangular
                y = betm * a[en][na];
                b[na][na] = -almi * b[en][en] / y;
                b[na][en] = (almr * b[en][en] - betm * a[en][en]) / y;
                b[en][na] = 0.0f;
                b[en][en] = 1.0f;
                enm2 = na;
                if (enm2 == 0) goto L795;

                // for i=en-2 step -1 until 1 do --
                for (ii = 0; ii < enm2; ++ii)
                {
                    i = na - ii - 1;
                    w = betm * a[i][i] - almr * b[i][i];
                    w1 = -almi * b[i][i];
                    ra = 0.0f;
                    sa = 0.0f;

                    for (j = m; j <= en; ++j)
                    {
                        x = betm * a[i][j] - almr * b[i][j];
                        x1 = -almi * b[i][j];
                        ra = ra + x * b[j][na] - x1 * b[j][en];
                        sa = sa + x * b[j][en] + x1 * b[j][na];
                    }

                    if (i == 0 || isw == 2) goto L770;
                    if (betm * a[i][i - 1] == 0.0) goto L770;

                    zz = w;
                    z1 = w1;
                    r = ra;
                    s = sa;
                    isw = 2;
                    goto L790;

                L770:
                    m = i;
                    if (isw == 2) goto L780;

                    // Complex 1-by-1 block 
                    tr = -ra;
                    ti = -sa;

                L773:
                    dr = w;
                    di = w1;

                // Complex divide (t1,t2) = (tr,ti) / (dr,di)
                L775:
                    if (Math.Abs(di) > Math.Abs(dr)) goto L777;
                    rr = di / dr;
                    d = dr + di * rr;
                    t1 = (tr + ti * rr) / d;
                    t2 = (ti - tr * rr) / d;

                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                L777:
                    rr = dr / di;
                    d = dr * rr + di;
                    t1 = (tr * rr + ti) / d;
                    t2 = (ti * rr - tr) / d;
                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                // Complex 2-by-2 block 
                L780:
                    x = betm * a[i][i + 1] - almr * b[i][i + 1];
                    x1 = -almi * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    tr = y * ra - w * r + w1 * s;
                    ti = y * sa - w * s - w1 * r;
                    dr = w * zz - w1 * z1 - x * y;
                    di = w * z1 + w1 * zz - x1 * y;
                    if (dr == 0.0 && di == 0.0)
                        dr = epsb;
                    goto L775;

                L782:
                    b[i + 1][na] = t1;
                    b[i + 1][en] = t2;
                    isw = 1;
                    if (Math.Abs(y) > Math.Abs(w) + Math.Abs(w1))
                        goto L785;
                    tr = -ra - x * b[(i + 1)][na] + x1 * b[(i + 1)][en];
                    ti = -sa - x * b[(i + 1)][en] - x1 * b[(i + 1)][na];
                    goto L773;

                L785:
                    t1 = (-r - zz * b[(i + 1)][na] + z1 * b[(i + 1)][en]) / y;
                    t2 = (-s - zz * b[(i + 1)][en] - z1 * b[(i + 1)][na]) / y;

                L787:
                    b[i][na] = t1;
                    b[i][en] = t2;

                L790:
                    ;
                }

            // End complex vector
            L795:
                isw = 3 - isw;

            L800:
                ;
            }

            // End back substitution. Transform to original coordinate system.
            for (jj = 0; jj < n; ++jj)
            {
                j = n - jj - 1;

                for (i = 0; i < n; ++i)
                {
                    zz = 0.0f;
                    for (k = 0; k <= j; ++k)
                        zz += z[i][k] * b[k][j];
                    z[i][j] = zz;
                }
            }

            // Normalize so that modulus of largest component of each vector is 1.
            // (isw is 1 initially from before)
            for (j = 0; j < n; ++j)
            {
                d = 0.0f;
                if (isw == 2) goto L920;
                if (alfi[j] != 0.0) goto L945;

                for (i = 0; i < n; ++i)
                {
                    if ((Math.Abs(z[i][j])) > d)
                        d = (Math.Abs(z[i][j]));
                }

                for (i = 0; i < n; ++i)
                    z[i][j] /= d;

                goto L950;

            L920:
                for (i = 0; i < n; ++i)
                {
                    r = System.Math.Abs(z[i][j - 1]) + System.Math.Abs(z[i][j]);
                    if (r != 0.0)
                    {
                        // Computing 2nd power
                        float u1 = z[i][j - 1] / r;
                        float u2 = z[i][j] / r;
                        r *= Maths.Sqrt(u1 * u1 + u2 * u2);
                    }
                    if (r > d)
                        d = r;
                }

                for (i = 0; i < n; ++i)
                {
                    z[i][j - 1] /= d;
                    z[i][j] /= d;
                }

            L945:
                isw = 3 - isw;

            L950:
                ;
            }

            return;
        }
        /// <summary>
        /// Estimates unit round-off (machine epsilon scaled by |x|) using the classic EISPACK strategy.
        /// </summary>
        /// <remarks>
        /// Port of EISPACK <c>EPSLON</c>. It returns an estimate of the smallest ε such that 1+ε ≠ 1,
        /// multiplied by <paramref name="x"/> in magnitude. Used when a caller passes 0 to request
        /// “use machine precision” as the iteration tolerance.
        /// </remarks>
        /// <param name="x">Reference magnitude. The result is approximately eps_machine * |x|</param>
        /// <returns>Estimated round-off scaled to |x|</returns>
        private static float Epsilon(float x)
        {
            float a, b, c, eps;

            a = 1.3333333333333333f;

        L10:
            b = a - 1.0f;
            c = b + b + b;
            eps = Math.Abs(c - 1.0f);

            if (eps == 0.0)
                goto L10;

            return eps * System.Math.Abs(x);
        }
        /// <summary>
        /// Returns <paramref name="a"/> with the sign of <paramref name="b"/> (sign transfer / copysign).
        /// </summary>
        /// <remarks>
        /// Equivalent to C/C++ <c>copysignf</c>: result is |a| if b ≥ 0, otherwise −|a|.
        /// Used throughout to stabilize reflector/rotation constructions.
        /// </remarks>
        /// <param name="a">Magnitude donor</param>
        /// <param name="b">Sign donor</param>
        /// <returns>|a| with the sign of b</returns>
        private static float Sign(float a, float b)
        {
            float x = a >= 0 ? a : -a;
            return b >= 0 ? x : -x;
        }
        #endregion
    }
}
