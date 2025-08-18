using System;
using System.Text;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a Pade approximant.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Pad%C3%A9_approximant.
    /// </remarks>
    /// Example: exp(x) = 1 + x + x^2/2! + x^3/3! + ...
    /// <code>
    /// float[] taylorExp = new float[] { 1.0f, 1.0f, 0.5f, 1.0f/6.0f, 1.0f/24.0f, 1.0f/120.0f, 1.0f/720.0f };
    /// </code>
    /// </summary>
    [Serializable]
    public class Pade
    {
        #region Private data
        private int m;
        private int n;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the Pade approximant.
        /// <param name="m">The degree of the numerator of a rational function</param>
        /// <param name="n">The degree of the denominator of a rational function</param>
        /// </summary>
        public Pade(int m = 2, int n = 2)
        {
            this.M = m;
            this.N = n;
        }
        /// <summary>
        /// Gets or sets the degree of the numerator of a rational function.
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the degree of the denominator of a rational function
        /// </summary>
        public int N
        {
            get
            {
                return this.n;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Invalid argument value");

                this.n = value;
            }
        }
        /// <summary>
        /// Returns the Pade approximant.
        /// </summary>
        /// <param name="taylorCoeffs">Taylor series coefficients</param>
        /// <exception cref="ArgumentException">Exception</exception>
        /// <returns>Coeffs</returns>
        public (float[] NumeratorCoeffs, float[] DenominatorCoeffs) Compute(float[] taylorCoeffs)
        {
            if (taylorCoeffs.Length < m + n + 1)
                throw new ArgumentException("Not enough Taylor series coefficients for the specified orders");

            float[,] A = new float[n, n];
            float[] b = new float[n];

            for (int i = 0; i < n; i++)
            {
                b[i] = -taylorCoeffs[m + i + 1];

                for (int j = 0; j < n; j++)
                {
                    A[i, j] = taylorCoeffs[m + i - j];
                }
            }

            float[] q = Matrice.Solve(A, b);

            float[] denominatorCoeffs = new float[n + 1];
            denominatorCoeffs[0] = 1.0f;

            for (int i = 0; i < n; i++)
            {
                denominatorCoeffs[i + 1] = q[i];
            }

            float[] numeratorCoeffs = new float[m + 1];

            for (int k = 0; k <= m; k++)
            {
                float sum = 0.0f;

                for (int j = 0; j <= Math.Min(k, n); j++)
                {
                    sum += denominatorCoeffs[j] * taylorCoeffs[k - j];
                }
                numeratorCoeffs[k] = sum;
            }

            return (numeratorCoeffs, denominatorCoeffs);
        }
        /// <summary>
        /// Returns the Pade approximant.
        /// </summary>
        /// <param name="taylorCoeffs">Taylor series coefficients</param>
        /// <exception cref="ArgumentException">Exception</exception>
        /// <returns>Coeffs</returns>
        public (Complex32[] NumeratorCoeffs, Complex32[] DenominatorCoeffs) Compute(Complex32[] taylorCoeffs)
        {
            if (taylorCoeffs.Length < m + n + 1)
                throw new ArgumentException("Not enough Taylor series coefficients for the specified orders");

            Complex32[,] A = new Complex32[n, n];
            Complex32[] b = new Complex32[n];

            for (int i = 0; i < n; i++)
            {
                b[i] = -taylorCoeffs[m + i + 1];

                for (int j = 0; j < n; j++)
                {
                    A[i, j] = taylorCoeffs[m + i - j];
                }
            }

            Complex32[] q = Matrice.Solve(A, b);

            Complex32[] denominatorCoeffs = new Complex32[n + 1];
            denominatorCoeffs[0] = 1.0f;

            for (int i = 0; i < n; i++)
            {
                denominatorCoeffs[i + 1] = q[i];
            }

            Complex32[] numeratorCoeffs = new Complex32[m + 1];

            for (int k = 0; k <= m; k++)
            {
                Complex32 sum = 0.0f;

                for (int j = 0; j <= Math.Min(k, n); j++)
                {
                    sum += denominatorCoeffs[j] * taylorCoeffs[k - j];
                }
                numeratorCoeffs[k] = sum;
            }

            return (numeratorCoeffs, denominatorCoeffs);
        }
        /// <summary>
        /// Evaluates a function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns>Value</returns>
        public float Compute(float x, float[] numeratorCoeffs, float[] denominatorCoeffs)
        {
            float num = 0;
            float den = 0;

            for (int i = numeratorCoeffs.Length - 1; i >= 0; i--)
            {
                num = num * x + numeratorCoeffs[i];
            }

            for (int i = denominatorCoeffs.Length - 1; i >= 0; i--)
            {
                den = den * x + denominatorCoeffs[i];
            }

            return num / den;
        }
        /// <summary>
        /// Evaluates a function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns>Value</returns>
        public Complex32 Compute(Complex32 x, float[] numeratorCoeffs, float[] denominatorCoeffs)
        {
            Complex32 num = 0;
            Complex32 den = 0;

            for (int i = numeratorCoeffs.Length - 1; i >= 0; i--)
            {
                num = num * x + numeratorCoeffs[i];
            }

            for (int i = denominatorCoeffs.Length - 1; i >= 0; i--)
            {
                den = den * x + denominatorCoeffs[i];
            }

            return num / den;
        }
        /// <summary>
        /// Evaluates a function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns>Value</returns>
        public Complex32 Compute(float x, Complex32[] numeratorCoeffs, Complex32[] denominatorCoeffs)
        {
            Complex32 num = 0;
            Complex32 den = 0;

            for (int i = numeratorCoeffs.Length - 1; i >= 0; i--)
            {
                num = num * x + numeratorCoeffs[i];
            }

            for (int i = denominatorCoeffs.Length - 1; i >= 0; i--)
            {
                den = den * x + denominatorCoeffs[i];
            }

            return num / den;
        }
        /// <summary>
        /// Evaluates a function.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns>Value</returns>
        public Complex32 Compute(Complex32 x, Complex32[] numeratorCoeffs, Complex32[] denominatorCoeffs)
        {
            Complex32 num = 0;
            Complex32 den = 0;

            for (int i = numeratorCoeffs.Length - 1; i >= 0; i--)
            {
                num = num * x + numeratorCoeffs[i];
            }

            for (int i = denominatorCoeffs.Length - 1; i >= 0; i--)
            {
                den = den * x + denominatorCoeffs[i];
            }

            return num / den;
        }
        /// <summary>
        /// Returns the equation of a Pade approximant represented as a string.
        /// </summary>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns></returns>
        /// <exception cref="ArgumentNullException">Exception</exception>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public string Equation(float[] numeratorCoeffs, float[] denominatorCoeffs)
        {
            if (numeratorCoeffs == null || denominatorCoeffs == null)
                throw new ArgumentNullException();

            const float eps = 1e-12f;

            string num = FormatPolynomial(numeratorCoeffs, "X", eps);
            string den = FormatPolynomial(denominatorCoeffs, "X", eps);

            if (string.IsNullOrEmpty(num)) num = "0";
            if (string.IsNullOrEmpty(den)) den = "0";

            return $"({num}) / ({den})";
        }
        /// <summary>
        /// Formats polynomial c0 + c1*x + ... + cK*x^K as a readable string.
        /// </summary>
        /// <param name="c">Polynomial</param>
        /// <param name="var">Variable</param>
        /// <param name="eps">Epsilon</param>
        private static string FormatPolynomial(float[] c, string var, float eps)
        {
            var sb = new StringBuilder();

            // find highest non-negligible power
            int kmax = -1;
            for (int k = c.Length - 1; k >= 0; k--)
                if (Math.Abs(c[k]) > eps) { kmax = k; break; }

            if (kmax < 0) return string.Empty; // all ~zero

            for (int k = kmax; k >= 0; k--)
            {
                float a = c[k];
                if (Math.Abs(a) <= eps) continue;

                bool firstTerm = sb.Length == 0;

                // sign
                if (!firstTerm)
                {
                    sb.Append(a >= 0 ? " + " : " - ");
                }
                else
                {
                    if (a < 0) sb.Append("-");
                }

                float absA = Math.Abs(a);

                // build term
                if (k == 0)
                {
                    // constant term: always print coefficient
                    sb.Append(absA.ToString());
                }
                else
                {
                    // coefficient (hide 1)
                    if (absA > 1 + 1e-12 || absA < 1 - 1e-12)
                    {
                        sb.Append(absA.ToString());
                    }

                    sb.Append(var);

                    if (k >= 2)
                    {
                        sb.Append("^");
                        sb.Append(k.ToString());
                    }
                }
            }

            return sb.ToString();
        }
        /// <summary>
        /// Returns the equation of a Pade approximant represented as a string.
        /// </summary>
        /// <param name="numeratorCoeffs">Numerator coeffs</param>
        /// <param name="denominatorCoeffs">Denominator coeffs</param>
        /// <returns></returns>
        /// <exception cref="ArgumentNullException">Exception</exception>
        /// <returns>Text as a sequence of Unicode characters</returns>
        public string Equation(Complex32[] numeratorCoeffs, Complex32[] denominatorCoeffs)
        {
            if (numeratorCoeffs == null || denominatorCoeffs == null)
                throw new ArgumentNullException();

            const float eps = 1e-12f;

            string num = FormatPolynomial(numeratorCoeffs, "X", eps);
            string den = FormatPolynomial(denominatorCoeffs, "X", eps);

            if (string.IsNullOrEmpty(num)) num = "0";
            if (string.IsNullOrEmpty(den)) den = "0";

            return $"({num}) / ({den})";
        }
        /// <summary>
        /// Formats polynomial c0 + c1*x + ... + cK*x^K as a readable string.
        /// </summary>
        /// <param name="c">Polynomial</param>
        /// <param name="var">Variable</param>
        /// <param name="eps">Epsilon</param>
        private static string FormatPolynomial(Complex32[] c, string var, float eps)
        {
            var sb = new StringBuilder();

            // highest non-negligible power
            int kmax = -1;
            for (int k = c.Length - 1; k >= 0; k--)
                if (!IsZero(c[k], eps)) { kmax = k; break; }

            if (kmax < 0) return string.Empty; // all ~zero

            for (int k = kmax; k >= 0; k--)
            {
                var a = c[k];
                if (IsZero(a, eps)) continue;

                bool first = sb.Length == 0;

                float re = a.Real;
                float im = a.Imag;

                bool isRealish = Math.Abs(im) <= eps;
                bool isImagOnly = Math.Abs(re) <= eps && Math.Abs(im) > eps;

                // ---------- SIGN handling between terms ----------
                if (!first)
                {
                    if (isRealish)
                    {
                        sb.Append(re >= 0 ? " + " : " - ");
                    }
                    else if (isImagOnly)
                    {
                        sb.Append(im >= 0 ? " + " : " - ");
                    }
                    else
                    {
                        // general complex: keep '+' and let the coefficient carry inner sign
                        sb.Append(" + ");
                    }
                }
                else
                {
                    if (isRealish && re < 0) sb.Append("-");
                    if (isImagOnly && im < 0) sb.Append("-");
                    // for general complex first term, no leading sign; it will be inside parentheses if needed
                }

                // ---------- TERM body ----------
                // absolute value (for printing when we already emitted a sign)
                float absRe = Math.Abs(re);
                float absIm = Math.Abs(im);

                if (k == 0)
                {
                    // constant term: print full coefficient (with magnitude already signed above if pure real/imag)
                    if (isRealish)
                    {
                        sb.Append(absRe.ToString());
                    }
                    else if (isImagOnly)
                    {
                        if (NearlyOne(absIm, eps)) sb.Append("i");
                        else { sb.Append(absIm.ToString()).Append("i"); }
                    }
                    else
                    {
                        sb.Append(a);
                    }
                }
                else
                {
                    // non-constant: may hide ±1 (only when imag≈0)
                    bool isPlusMinusOneReal = isRealish && NearlyOne(absRe, eps);

                    if (isPlusMinusOneReal)
                    {
                        // coefficient is ±1 (imag≈0):
                        // we already printed a sign; so just the variable/power
                        sb.Append(var);
                        if (k >= 2) sb.Append("^").Append(k.ToString());
                    }
                    else if (isImagOnly && NearlyOne(absIm, eps))
                    {
                        // pure ±i: we printed sign already; just "ix^k"
                        sb.Append("i").Append(var);
                        if (k >= 2) sb.Append("^").Append(k.ToString());
                    }
                    else
                    {
                        // general coefficient (real, pure imag with |im|!=1, or complex)
                        string coeffText;

                        if (isRealish)
                            coeffText = absRe.ToString();
                        else if (isImagOnly)
                            coeffText = absIm.ToString() + "i";
                        else
                            coeffText = a.ToString();

                        sb.Append(coeffText).Append(var);
                        if (k >= 2) sb.Append("^").Append(k.ToString());
                    }
                }
            }

            return sb.ToString();
        }
        /// <summary>
        /// True if |z|≈0 (both parts small).
        /// </summary>
        /// <param name="z">Value</param>
        /// <param name="eps">Epsilon</param>
        /// <returns>Boolean</returns>
        private static bool IsZero(Complex32 z, float eps) => Math.Abs(z.Real) <= eps && Math.Abs(z.Imag) <= eps;
        /// <summary>
        /// True if |x - 1| ≤ eps.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="eps">Epsilon</param>
        /// <returns>Boolean</returns>
        private static bool NearlyOne(float x, float eps) => Math.Abs(x - 1f) <= eps;
        #endregion
    }

}
