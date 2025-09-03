using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines the least squares approximation class.
    /// <remarks>
    /// This class is a solution to the problem of finding the function A (x) ≈ F (x), where F (x) is the original function.
    /// More information can be found on the website:
    /// http://simenergy.ru/math-analysis/digital-processing/85-ordinary_least_squares
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Approximation
    {
        #region Private data
        private ApproximationMethod method;
        private int power;
        #endregion

        #region Approximation components
        /// <summary>
        /// Initializes the least squares approximation class.
        /// </summary>
        /// <param name="power">Polynomial degree</param>
        /// <param name="method">Approximation method</param>
        public Approximation(int power = 1, ApproximationMethod method = ApproximationMethod.Polynomial)
        {
            this.Power = power;
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the degree of the polynomial.
        /// </summary>
        public int Power
        {
            get
            {
                return this.power;
            }
            set
            {
                if (value < 1)
                    throw new ArgumentException("Invalid argument value");

                this.power = value;
            }
        }
        /// <summary>
        /// Gets or sets the approximation method.
        /// </summary>
        public ApproximationMethod MethodType
        {
            get
            {
                return this.method;
            }
            set
            {
                this.method = value;
            }
        }
        #endregion

        #region Public voids
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <returns>Array</returns>
        public float[] Compute(float[] x, float[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out _, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out _, out _, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out _, out _, out _);

                default:
                    return Approximation.Powr(x, y, power, out _, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public float[] Compute(float[] x, float[] y, out float[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out _, out _);

                default:
                    return Approximation.Powr(x, y, power, out cf, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <returns>Array</returns>
        public float[] Compute(float[] x, float[] y, out float[] cf, out float similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out similarity, out _);

                default:
                    return Approximation.Powr(x, y, power, out cf, out similarity, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public float[] Compute(float[] x, float[] y, out float[] cf, out float similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out similarity, out equation);

                default:
                    return Approximation.Powr(x, y, power, out cf, out similarity, out equation);
            }
        }

        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <returns>Array</returns>
        public ComplexF[] Compute(ComplexF[] x, ComplexF[] y)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out _, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out _, out _, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out _, out _, out _);

                default:
                    return Approximation.Powr(x, y, power, out _, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <returns>Array</returns>
        public ComplexF[] Compute(ComplexF[] x, ComplexF[] y, out ComplexF[] cf)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out _, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out _, out _);

                default:
                    return Approximation.Powr(x, y, power, out cf, out _, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <returns>Array</returns>
        public ComplexF[] Compute(ComplexF[] x, ComplexF[] y, out ComplexF[] cf, out ComplexF similarity)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out similarity, out _);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out similarity, out _);

                default:
                    return Approximation.Powr(x, y, power, out cf, out similarity, out _);
            }
        }
        /// <summary>
        /// Returns the approximation value.
        /// </summary>
        /// <param name="x">Array of argument values</param>
        /// <param name="y">Array of function values</param>
        /// <param name="cf">Approximation coefficients</param>
        /// <param name="similarity">Similarity</param>
        /// <param name="equation">Equation</param>
        /// <returns>Array</returns>
        public ComplexF[] Compute(ComplexF[] x, ComplexF[] y, out ComplexF[] cf, out ComplexF similarity, out string equation)
        {
            // chose method of approximation
            switch (method)
            {
                case ApproximationMethod.Polynomial:
                    return Approximation.Poly(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Logarithmic:
                    return Approximation.Logc(x, y, power, out cf, out similarity, out equation);

                case ApproximationMethod.Exponential:
                    return Approximation.Expn(x, y, power, out cf, out similarity, out equation);

                default:
                    return Approximation.Powr(x, y, power, out cf, out similarity, out equation);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Polynomial least-squares fit: y(x) ≈ Σ_{k=0}^{power} c_k x^k.
        /// </summary>
        /// <remarks>
        /// Builds a degree-<paramref name="power"/> polynomial model using least squares.
        /// The number of basis terms is m = power + 1. Returns the fitted values ŷ at the
        /// sample points <paramref name="x"/> and also outputs the coefficients, an error metric,
        /// and a human-readable equation string.
        /// </remarks>
        /// <param name="x">Sample abscissas</param>
        /// <param name="y">Sample ordinates</param>
        /// <param name="power">Polynomial degree (≥ 1)</param>
        /// <param name="cf">Output: polynomial coefficients c[0..power]</param>
        /// <param name="error">Output: fit error (as computed by LeastSquaresOptions.Error)</param>
        /// <param name="equation">Output: formatted equation string</param>
        /// <returns>Fitted values ŷ at points x (same length as <paramref name="y"/>)</returns>
        private static float[] Poly(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            float[] ya = LeastSquaresOptions.Polynomial(x, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf);
            return ya;
        }
        /// <summary>
        /// Polynomial least-squares fit: y(x) ≈ Σ_{k=0}^{power} c_k x^k.
        /// </summary>
        /// <remarks>
        /// Builds a degree-<paramref name="power"/> polynomial model using least squares.
        /// The number of basis terms is m = power + 1. Returns the fitted values ŷ at the
        /// sample points <paramref name="x"/> and also outputs the coefficients, an error metric,
        /// and a human-readable equation string.
        /// </remarks>
        /// <param name="x">Sample abscissas</param>
        /// <param name="y">Sample ordinates</param>
        /// <param name="power">Polynomial degree (≥ 1)</param>
        /// <param name="cf">Output: polynomial coefficients c[0..power]</param>
        /// <param name="error">Output: fit error (as computed by LeastSquaresOptions.Error)</param>
        /// <param name="equation">Output: formatted equation string</param>
        /// <returns>Fitted values ŷ at points x (same length as <paramref name="y"/>)</returns>
        private static ComplexF[] Poly(ComplexF[] x, ComplexF[] y, int power, out ComplexF[] cf, out ComplexF error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            cf = LeastSquaresOptions.Coefficients(x, y, m);
            ComplexF[] ya = LeastSquaresOptions.Polynomial(x, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf);
            return ya;
        }
        /// <summary>
        /// Logarithmic least-squares fit (real): y(x) ≈ Σ_{k=0}^{power} c_k [log(x)]^k.
        /// </summary>
        /// <remarks>
        /// Fits a polynomial in the transformed variable u = log(x).
        /// <para><b>Domain:</b> requires x[i] &gt; 0 for all i.</para>
        /// </remarks>
        /// <param name="x">Sample abscissas (must be &gt; 0)</param>
        /// <param name="y">Sample ordinates</param>
        /// <param name="power">Degree in log-domain (≥ 1)</param>
        /// <param name="cf">Output: coefficients in the log-domain</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string using “* Log(x)^k”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static float[] Logc(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            float[] xa = new float[n];
            float[] ya;

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = MathF.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * Log(x)^");
            return ya;
        }
        /// <summary>
        /// Logarithmic least-squares fit (real): y(x) ≈ Σ_{k=0}^{power} c_k [log(x)]^k.
        /// </summary>
        /// <remarks>
        /// Fits a polynomial in the transformed variable u = log(x).
        /// <para><b>Domain:</b> requires x[i] &gt; 0 for all i.</para>
        /// </remarks>
        /// <param name="x">Sample abscissas (must be &gt; 0)</param>
        /// <param name="y">Sample ordinates</param>
        /// <param name="power">Degree in log-domain (≥ 1)</param>
        /// <param name="cf">Output: coefficients in the log-domain</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string using “* Log(x)^k”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static ComplexF[] Logc(ComplexF[] x, ComplexF[] y, int power, out ComplexF[] cf, out ComplexF error, out string equation)
        {
            // Options:
            int n = x.Length, i;
            int m = (power < 1) ? 2 : power + 1;
            ComplexF[] xa = new ComplexF[n];
            ComplexF[] ya;

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = MathF.Log(x[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, y, m);
            ya = LeastSquaresOptions.Polynomial(xa, cf);
            error = LeastSquaresOptions.Error(ya, y);
            equation = LeastSquaresOptions.Equation(cf, " * Log(x)^");
            return ya;
        }
        /// <summary>
        /// Exponential least-squares fit (real): y(x) ≈ exp( Σ_{k=0}^{power} c_k x^k ).
        /// </summary>
        /// <remarks>
        /// Applies a log transform to the response: v = log(y), fits a polynomial v ≈ Σ c_k x^k,
        /// then maps back by exponentiation ŷ = exp( Σ c_k x^k ).
        /// <para><b>Domain:</b> requires y[i] &gt; 0 for all i (log defined).</para>
        /// </remarks>
        /// <param name="x">Sample abscissas</param>
        /// <param name="y">Sample ordinates (must be &gt; 0)</param>
        /// <param name="power">Polynomial degree inside the exponent (≥ 1)</param>
        /// <param name="cf">Output: polynomial coefficients for log(y)</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string “Exp( … )”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static float[] Expn(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            float[] ya = new float[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Log(y[i], MathF.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            float[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Pow(MathF.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf)})";
            return ya;
        }
        /// <summary>
        /// Exponential least-squares fit (real): y(x) ≈ exp( Σ_{k=0}^{power} c_k x^k ).
        /// </summary>
        /// <remarks>
        /// Applies a log transform to the response: v = log(y), fits a polynomial v ≈ Σ c_k x^k,
        /// then maps back by exponentiation ŷ = exp( Σ c_k x^k ).
        /// <para><b>Domain:</b> requires y[i] &gt; 0 for all i (log defined).</para>
        /// </remarks>
        /// <param name="x">Sample abscissas</param>
        /// <param name="y">Sample ordinates (must be &gt; 0)</param>
        /// <param name="power">Polynomial degree inside the exponent (≥ 1)</param>
        /// <param name="cf">Output: polynomial coefficients for log(y)</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string “Exp( … )”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static ComplexF[] Expn(ComplexF[] x, ComplexF[] y, int power, out ComplexF[] cf, out ComplexF error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            ComplexF[] ya = new ComplexF[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Log(y[i], MathF.E);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(x, ya, m);
            ComplexF[] p = LeastSquaresOptions.Polynomial(x, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Pow(MathF.E, p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf)})";
            return ya;
        }
        /// <summary>
        /// Power-law–type least-squares fit (real): y(x) ≈ Exp( Σ_{k=0}^{power} c_k [Log(x)]^k ).
        /// </summary>
        /// <remarks>
        /// Applies log transforms to both variables: u = log(x), v = log(y), fits v ≈ Σ c_k u^k,
        /// then maps back as ŷ = exp( Σ c_k [log(x)]^k ).
        /// <para><b>Domain:</b> requires x[i] &gt; 0 and y[i] &gt; 0 for all i.</para>
        /// </remarks>
        /// <param name="x">Sample abscissas (must be &gt; 0)</param>
        /// <param name="y">Sample ordinates (must be &gt; 0)</param>
        /// <param name="power">Polynomial degree in log-log domain (≥ 1)</param>
        /// <param name="cf">Output: coefficients in the log-log domain</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string using “Exp( … * Log(x)^k … )”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static float[] Powr(float[] x, float[] y, int power, out float[] cf, out float error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            float[] xa = new float[n];
            float[] ya = new float[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = MathF.Log(x[i]);
                ya[i] = MathF.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            float[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf, " * Log(x)^")})";
            return ya;
        }
        /// <summary>
        /// Power-law–type least-squares fit (real): y(x) ≈ Exp( Σ_{k=0}^{power} c_k [Log(x)]^k ).
        /// </summary>
        /// <remarks>
        /// Applies log transforms to both variables: u = log(x), v = log(y), fits v ≈ Σ c_k u^k,
        /// then maps back as ŷ = exp( Σ c_k [log(x)]^k ).
        /// <para><b>Domain:</b> requires x[i] &gt; 0 and y[i] &gt; 0 for all i.</para>
        /// </remarks>
        /// <param name="x">Sample abscissas (must be &gt; 0)</param>
        /// <param name="y">Sample ordinates (must be &gt; 0)</param>
        /// <param name="power">Polynomial degree in log-log domain (≥ 1)</param>
        /// <param name="cf">Output: coefficients in the log-log domain</param>
        /// <param name="error">Output: fit error on original y</param>
        /// <param name="equation">Output: equation string using “Exp( … * Log(x)^k … )”</param>
        /// <returns>Fitted values ŷ at x</returns>
        private static ComplexF[] Powr(ComplexF[] x, ComplexF[] y, int power, out ComplexF[] cf, out ComplexF error, out string equation)
        {
            // Options:
            int m = (power < 1) ? 2 : power + 1;
            int n = x.Length, i;
            ComplexF[] xa = new ComplexF[n];
            ComplexF[] ya = new ComplexF[n];

            // log-scale:
            for (i = 0; i < n; i++)
            {
                xa[i] = MathF.Log(x[i]);
                ya[i] = MathF.Log(y[i]);
            }

            // approximation:
            cf = LeastSquaresOptions.Coefficients(xa, ya, m);
            ComplexF[] p = LeastSquaresOptions.Polynomial(xa, cf);

            // exponential-scale:
            for (i = 0; i < n; i++)
            {
                ya[i] = MathF.Exp(p[i]);
            }

            error = LeastSquaresOptions.Error(ya, y);
            equation = $"Exp({LeastSquaresOptions.Equation(cf, " * Log(x)^")})";
            return ya;
        }
        #endregion
    }
}
