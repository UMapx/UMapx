using System;
using UMapx.Core;
using UMapx.Decomposition;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class for solving equations using the spectral decomposition of a matrix.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.mathworks.com/help/matlab/ref/roots.html
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Roots
    {
        #region Private data
        private EVD eig;
        private double eps;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class of equations using the spectral decomposition of a matrix.
        /// </summary>
        /// <param name="eps">Epsilon [0, 1]</param>
        public Roots(double eps = 1e-16)
        {
            this.Eps = eps;
        }
        /// <summary>
        /// Gets or sets an error [0, 1].
        /// </summary>
        public double Eps
        {
            get
            {
                return this.eps;
            }
            set
            {
                this.eps = Maths.Double(value);
            }
        }
        /// <summary>
        /// Returns a column vector corresponding to the numerical solution of the polynomial: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="polynomial">Polynomial</param>
        /// <returns>Array</returns>
        public Complex[] Compute(double[] polynomial)
        {
            // MATLAB roots method
            // represented by Valery Asiryan, 2018.
            // properties of polynomial:
            int length = polynomial.Length;
            int i, index = -1;

            // finding non-zero element:
            for (i = 0; i < length; i++)
            {
                if (polynomial[i] != 0)
                {
                    index = i;
                    break;
                }
            }

            // return null array:
            if (index == -1)
            {
                return new Complex[0];
            }

            // get scaling factor:
            int m = length - index - 1;
            double scale = polynomial[index];
            double[] c = new double[m];

            // create new polynomial:
            for (i = 0; i < m; i++)
            {
                c[i] = polynomial[i + index + 1] / scale;
            }

            // Eigen-value decomposition for
            // companion matrix:
            eig = new EVD(Matrice.Companion(c), this.eps);

            // Complex result:
            return eig.D;
        }
        /// <summary>
        /// Returns a column vector of polynomial coefficients: p(1)*x^n + ... + p(n)*x + p(n+1) = 0.
        /// </summary>
        /// <param name="roots">Roots</param>
        /// <returns>Array</returns>
        public double[] Compute(Complex[] roots)
        {
            // MATLAB roots method
            // represented by Valery Asiryan, 2018.
            // properties of polynomial:
            int length = roots.Length, m = length + 1, j, i;

            // arrays:
            Complex[] v = new Complex[length];
            Complex[] p = new Complex[m];

            // point:
            p[0] = 1.0;

            // create new polynomial:
            for (j = 0; j < length; j++)
            {
                // right part:
                for (i = 0; i <= j; i++)
                {
                    v[i] = roots[j] * p[i];
                }
                // left part:
                for (i = 0; i <= j; i++)
                {
                    p[i + 1] -= v[i];
                }
            }

            // Real result:
            return p.Real();
        }
        #endregion
    }
}
