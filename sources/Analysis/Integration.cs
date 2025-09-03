using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements numerical integration.
    /// <remarks>
    /// This class is a solution to the problem of finding the value of the integral of the function F(x) within the values of a and b.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Integration
    {
        #region Private data
        private IntegrationMethod method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements numerical integration.
        /// </summary>
        /// <param name="method">Integration method</param>
        public Integration(IntegrationMethod method = IntegrationMethod.Rectangle)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the integration method.
        /// </summary>
        public IntegrationMethod MethodType
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
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Value</returns>
        public float Compute(IFloat function, float a, float b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case IntegrationMethod.Midpoint:
                    return Integration.Midp(function, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.Trap(function, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.Simp(function, a, b, n);

                case IntegrationMethod.Romberg:
                    return Integration.Romb(function, a, b, n);

                default:
                    return Integration.Rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Value</returns>
        public float Compute(float[] y, float a, float b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case IntegrationMethod.Midpoint:
                    return Integration.Midp(y, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.Trap(y, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.Simp(y, a, b, n);

                case IntegrationMethod.Romberg:
                    return Integration.Romb(y, a, b, n);

                default:
                    return Integration.Rect(y, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="function">Continuous function delegate</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex32 Compute(IComplex32 function, Complex32 a, Complex32 b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case IntegrationMethod.Midpoint:
                    return Integration.Midp(function, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.Trap(function, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.Simp(function, a, b, n);

                case IntegrationMethod.Romberg:
                    return Integration.Romb(function, a, b, n);

                default:
                    return Integration.Rect(function, a, b, n);
            }
        }
        /// <summary>
        /// Returns the value of the integral of a function.
        /// </summary>
        /// <param name="y">Function vector</param>
        /// <param name="a">Lower limit</param>
        /// <param name="b">Upper limit</param>
        /// <param name="n">Number of splits</param>
        /// <returns>Complex number</returns>
        public Complex32 Compute(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            // chose method of integration
            switch (method)
            {
                case IntegrationMethod.Midpoint:
                    return Integration.Midp(y, a, b, n);

                case IntegrationMethod.Trapezoidal:
                    return Integration.Trap(y, a, b, n);

                case IntegrationMethod.Simpson:
                    return Integration.Simp(y, a, b, n);

                case IntegrationMethod.Romberg:
                    return Integration.Romb(y, a, b, n);

                default:
                    return Integration.Rect(y, a, b, n);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Left Riemann sum (rectangle rule) over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses left endpoints: x_i = a + i*h, h = (b - a)/n, i = 0..n-1.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals (rectangles)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Rect(IFloat f, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// Left Riemann sum (rectangle rule) for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Expects samples y[i] ≈ f(a + i*h) at left endpoints with h = (b - a)/n for i = 0..n-1.
        /// Assumes y.Length ≥ n and n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Sampled values y[i] at x_i = a + i*h (left endpoints)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals/samples used (left endpoints)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Rect(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// Midpoint rule over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses midpoints: x_i = a + (i + 0.5) * h, h = (b - a)/n, i = 0..n-1.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Midp(IFloat f, float a, float b, int n)
        {
            // Midpoint
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5f) * h);
            }
            return sum;
        }
        /// <summary>
        /// Trapezoidal rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid x_i = a + i*h with h = (b - a)/(n - 1).
        /// Applies the trapezoidal rule over adjacent pairs. Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Midp(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5f * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// Trapezoidal rule over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses endpoints of each subinterval: x_i = a + i*h, h = (b - a)/n.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Trap(IFloat f, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5f * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// Trapezoidal rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid with h = (b - a)/(n - 1).
        /// Applies trapezoids between successive samples. Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static float Trap(float[] y, float a, float b, int n)
        {
            float sum = 0.0f;
            float h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5f * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// Simpson's rule over a function f on [a, b] with n subintervals (n even preferred).
        /// </summary>
        /// <remarks>
        /// If n is even, applies Simpson 1/3 on all subintervals. If n is odd, applies a 3/8 segment
        /// over the first 3 subintervals and Simpson 1/3 over the remaining (n - 3).
        /// Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals (even for pure Simpson 1/3)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx (NaN if n &lt; 2)</returns>
        private static float Simp(IFloat f, float a, float b, int n)
        {
            if (n < 2) return float.NaN;
            float h = (b - a) / n;
            float sum;

            if (n % 2 == 0)
            {
                float s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n; i += 2) s_odd += f(a + i * h);
                for (int i = 2; i < n; i += 2) s_even += f(a + i * h);
                sum = (h / 3f) * (f(a) + f(b) + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (f(a) + 3f * f(a + h) + 3f * f(a + 2 * h) + f(a + 3 * h));
                int m = n - 3;
                float s_odd = 0f, s_even = 0f;
                for (int j = 1; j < m; j += 2) s_odd += f(a + (3 + j) * h);
                for (int j = 2; j < m; j += 2) s_even += f(a + (3 + j) * h);
                sum += (h / 3f) * (f(a + 3 * h) + f(b) + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }
        /// <summary>
        /// Simpson's rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid with spacing h = (b - a)/(n - 1).
        /// If (n - 1) is even, applies Simpson 1/3 on the full range; otherwise uses a 3/8 segment
        /// on the first three subintervals and Simpson 1/3 on the remainder. Requires n ≥ 3.
        /// Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx (NaN if n &lt; 3)</returns>
        private static float Simp(float[] y, float a, float b, int n)
        {
            if (n < 3) return float.NaN;

            float h = (b - a) / (n - 1);
            float sum;

            if ((n - 1) % 2 == 0)
            {
                float s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n - 1; i += 2) s_odd += y[i];
                for (int i = 2; i < n - 1; i += 2) s_even += y[i];
                sum = (h / 3f) * (y[0] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (y[0] + 3f * y[1] + 3f * y[2] + y[3]);
                float s_odd = 0f, s_even = 0f;
                for (int i = 4; i < n - 1; i += 2) s_odd += y[i - 1];
                for (int i = 5; i < n - 1; i += 2) s_even += y[i - 1];
                sum += (h / 3f) * (y[3] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }

        /// <summary>
        /// Romberg integration over a function f on [a, b] using Richardson extrapolation.
        /// </summary>
        /// <remarks>
        /// Builds the Romberg tableau R(k, j), k = 0..maxK-1, j = 0..k, starting from the trapezoidal
        /// rule and refining by halving the step each level. Stops early if successive diagonal entries
        /// satisfy a relative tolerance. Complexity per level grows geometrically.
        /// </remarks>
        /// <param name="f">Integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="maxK">Maximum number of Romberg levels (table size)</param>
        /// <param name="eps">Relative tolerance for early stopping</param>
        /// <returns>Romberg estimate of ∫_a^b f(x) dx</returns>
        private static float Romb(IFloat f, float a, float b, int maxK, float eps = 1e-8f)
        {
            if (maxK < 1) throw new ArgumentException();
            float[,] R = new float[maxK, maxK];

            float h = b - a;
            R[0, 0] = 0.5f * h * (f(a) + f(b));
            int n = 1;

            for (int k = 1; k < maxK; k++)
            {
                n *= 2; h *= 0.5f;
                float sum = 0f;
                for (int i = 1; i < n; i += 2) sum += f(a + i * h);
                R[k, 0] = 0.5f * R[k - 1, 0] + h * sum;

                float factor = 4f;
                for (int j = 1; j <= k; j++, factor *= 4f)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1f);

                if (Math.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * Math.Abs(R[k, k]))
                    return R[k, k];
            }
            return R[maxK - 1, maxK - 1];
        }
        /// <summary>
        /// Romberg integration for tabulated samples y on a uniform grid over [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes y[i] ≈ f(a + i*h0) on a uniform grid with h0 = (b - a)/(N - 1), N = y.Length.
        /// Builds a Romberg tableau using composite trapezoidal rules with M = 2^k panels,
        /// which requires that (N - 1) is divisible by 2^k at each level. The number of usable
        /// levels is limited by the highest power of two dividing (N - 1) and by <paramref name="maxK"/>.
        /// Early exit occurs when successive diagonal entries satisfy the relative tolerance.
        /// </remarks>
        /// <param name="y">Samples y[0..N-1] on a uniform grid from a to b (inclusive)</param>
        /// <param name="a">Lower limit of integration (corresponds to y[0])</param>
        /// <param name="b">Upper limit of integration (corresponds to y[N-1])</param>
        /// <param name="maxK">Maximum number of Romberg levels to build (≥ 1)</param>
        /// <param name="eps">Relative tolerance for early stopping</param>
        /// <returns>Romberg estimate of ∫_a^b f(x) dx</returns>
        private static float Romb(float[] y, float a, float b, int maxK, float eps = 1e-8f)
        {
            if (y == null || y.Length < 2) throw new ArgumentException("Function must have at least 2 samples");
            if (maxK < 1) throw new ArgumentException(nameof(maxK));
            int N = y.Length;
            int S = N - 1; // number of subintervals on the finest grid

            // Determine how many dyadic levels are actually available from the tabulated grid:
            int usableLevels = 1; // we can always build level k = 0 (M = 1)
            while (usableLevels < maxK && (S % (1 << usableLevels) == 0)) usableLevels++;
            // We'll build k = 0..usableLevels-1

            float[,] R = new float[usableLevels, usableLevels];

            for (int k = 0; k < usableLevels; k++)
            {
                int M = 1 << k;      // number of panels at level k
                int stride = S / M;  // index step to pick samples for this M
                float h = (b - a) / M;

                // Composite trapezoid with nodes every 'stride' samples:
                float interiorSum = 0f;
                for (int m = 1; m < M; m++)
                {
                    int idx = m * stride;   // guaranteed integer by construction
                    interiorSum += y[idx];
                }
                R[k, 0] = 0.5f * h * (y[0] + y[N - 1]) + h * interiorSum;

                // Richardson extrapolation along the row:
                float factor = 4f;
                for (int j = 1; j <= k; j++, factor *= 4f)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1f);

                if (k > 0 && Math.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * Math.Abs(R[k, k]))
                    return R[k, k];
            }

            return R[usableLevels - 1, usableLevels - 1];
        }

        /// <summary>
        /// Left Riemann sum (rectangle rule) over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses left endpoints: x_i = a + i*h, h = (b - a)/n, i = 0..n-1.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals (rectangles)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Rect(IComplex32 f, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + i * h);
            }
            return sum;
        }
        /// <summary>
        /// Left Riemann sum (rectangle rule) for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Expects samples y[i] ≈ f(a + i*h) at left endpoints with h = (b - a)/n for i = 0..n-1.
        /// Assumes y.Length ≥ n and n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Sampled values y[i] at x_i = a + i*h (left endpoints)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals/samples used (left endpoints)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Rect(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * y[i];
            }
            return sum;
        }

        /// <summary>
        /// Midpoint rule over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses midpoints: x_i = a + (i + 0.5) * h, h = (b - a)/n, i = 0..n-1.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Midp(IComplex32 f, Complex32 a, Complex32 b, int n)
        {
            // Midpoint
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += h * f(a + (i + 0.5) * h);
            }
            return sum;
        }
        /// <summary>
        /// Trapezoidal rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid x_i = a + i*h with h = (b - a)/(n - 1).
        /// Applies the trapezoidal rule over adjacent pairs. Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Midp(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += h * 0.5 * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// Trapezoidal rule over a function f on [a, b] with n subintervals.
        /// </summary>
        /// <remarks>
        /// Uses endpoints of each subinterval: x_i = a + i*h, h = (b - a)/n.
        /// Assumes n ≥ 1. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Trap(IComplex32 f, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / n;
            for (int i = 0; i < n; i++)
            {
                sum += 0.5 * h * (f(a + i * h) + f(a + (i + 1) * h));
            }
            return sum;
        }
        /// <summary>
        /// Trapezoidal rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid with h = (b - a)/(n - 1).
        /// Applies trapezoids between successive samples. Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx</returns>
        private static Complex32 Trap(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            Complex32 sum = 0.0;
            Complex32 h = (b - a) / (n - 1);
            for (int i = 0; i < (n - 1); i++)
            {
                sum += 0.5 * h * (y[i] + y[i + 1]);
            }
            return sum;
        }

        /// <summary>
        /// Simpson's rule over a function f on [a, b] with n subintervals (n even preferred).
        /// </summary>
        /// <remarks>
        /// If n is even, applies Simpson 1/3 on all subintervals. If n is odd, applies a 3/8 segment
        /// over the first 3 subintervals and Simpson 1/3 over the remaining (n - 3).
        /// Requires n ≥ 2. Complexity: O(n).
        /// </remarks>
        /// <param name="f">Continuous integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of subintervals (even for pure Simpson 1/3)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx (NaN if n &lt; 2)</returns>
        private static Complex32 Simp(IComplex32 f, Complex32 a, Complex32 b, int n)
        {
            if (n < 2) return Complex32.NaN;
            Complex32 h = (b - a) / n;
            Complex32 sum;

            if (n % 2 == 0)
            {
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n; i += 2) s_odd += f(a + i * h);
                for (int i = 2; i < n; i += 2) s_even += f(a + i * h);
                sum = (h / 3f) * (f(a) + f(b) + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (f(a) + 3f * f(a + h) + 3f * f(a + 2 * h) + f(a + 3 * h));
                int m = n - 3;
                Complex32 s_odd = 0f, s_even = 0f;
                for (int j = 1; j < m; j += 2) s_odd += f(a + (3 + j) * h);
                for (int j = 2; j < m; j += 2) s_even += f(a + (3 + j) * h);
                sum += (h / 3f) * (f(a + 3 * h) + f(b) + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }
        /// <summary>
        /// Simpson's rule for tabulated samples on [a, b].
        /// </summary>
        /// <remarks>
        /// Assumes n samples y[0..n-1] on a uniform grid with spacing h = (b - a)/(n - 1).
        /// If (n - 1) is even, applies Simpson 1/3 on the full range; otherwise uses a 3/8 segment
        /// on the first three subintervals and Simpson 1/3 on the remainder. Requires n ≥ 3.
        /// Complexity: O(n).
        /// </remarks>
        /// <param name="y">Samples y[i] at uniform points on [a, b]</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="n">Number of samples (grid points)</param>
        /// <returns>Approximation of ∫_a^b f(x) dx (NaN if n &lt; 3)</returns>
        private static Complex32 Simp(Complex32[] y, Complex32 a, Complex32 b, int n)
        {
            if (n < 3) return Complex32.NaN;

            Complex32 h = (b - a) / (n - 1);
            Complex32 sum;

            if ((n - 1) % 2 == 0)
            {
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 1; i < n - 1; i += 2) s_odd += y[i];
                for (int i = 2; i < n - 1; i += 2) s_even += y[i];
                sum = (h / 3f) * (y[0] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            else
            {
                sum = (3f * h / 8f) * (y[0] + 3f * y[1] + 3f * y[2] + y[3]);
                Complex32 s_odd = 0f, s_even = 0f;
                for (int i = 4; i < n - 1; i += 2) s_odd += y[i - 1];
                for (int i = 5; i < n - 1; i += 2) s_even += y[i - 1];
                sum += (h / 3f) * (y[3] + y[n - 1] + 4f * s_odd + 2f * s_even);
            }
            return sum;
        }

        /// <summary>
        /// Romberg integration over a function f on [a, b] using Richardson extrapolation.
        /// </summary>
        /// <remarks>
        /// Builds the Romberg tableau R(k, j), k = 0..maxK-1, j = 0..k, starting from the trapezoidal
        /// rule and refining by halving the step each level. Stops early if successive diagonal entries
        /// satisfy a relative tolerance. Complexity per level grows geometrically.
        /// </remarks>
        /// <param name="f">Integrand f(x)</param>
        /// <param name="a">Lower limit of integration</param>
        /// <param name="b">Upper limit of integration</param>
        /// <param name="maxK">Maximum number of Romberg levels (table size)</param>
        /// <param name="eps">Relative tolerance for early stopping</param>
        /// <returns>Romberg estimate of ∫_a^b f(x) dx</returns>
        private static Complex32 Romb(IComplex32 f, Complex32 a, Complex32 b, int maxK, float eps = 1e-8f)
        {
            if (maxK < 1) throw new ArgumentException();
            Complex32[,] R = new Complex32[maxK, maxK];

            Complex32 h = b - a;
            R[0, 0] = 0.5f * h * (f(a) + f(b));
            int n = 1;

            for (int k = 1; k < maxK; k++)
            {
                n *= 2; h *= 0.5f;
                Complex32 sum = 0f;
                for (int i = 1; i < n; i += 2) sum += f(a + i * h);
                R[k, 0] = 0.5f * R[k - 1, 0] + h * sum;

                Complex32 factor = 4f;
                for (int j = 1; j <= k; j++, factor *= 4f)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1f);

                if (MathF.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * MathF.Abs(R[k, k]))
                    return R[k, k];
            }
            return R[maxK - 1, maxK - 1];
        }
        /// <summary>
        /// Romberg integration for complex tabulated samples y on a uniform grid over [a, b].
        /// </summary>
        /// <remarks>
        /// Same assumptions as the real-valued overload: y[i] ≈ f(a + i*h0), uniform grid,
        /// and (N - 1) must be divisible by 2^k for each level k to be usable.
        /// Early exit uses a relative tolerance based on complex magnitude.
        /// </remarks>
        /// <param name="y">Complex samples y[0..N-1] on a uniform grid from a to b</param>
        /// <param name="a">Lower limit (complex), corresponds to y[0]</param>
        /// <param name="b">Upper limit (complex), corresponds to y[N-1]</param>
        /// <param name="maxK">Maximum number of Romberg levels (≥ 1)</param>
        /// <param name="eps">Relative tolerance for early stopping</param>
        /// <returns>Romberg estimate of ∫_a^b f(x) dx in Complex32</returns>
        private static Complex32 Romb(Complex32[] y, Complex32 a, Complex32 b, int maxK, float eps = 1e-8f)
        {
            if (y == null || y.Length < 2) throw new ArgumentException("y must have at least 2 samples.");
            if (maxK < 1) throw new ArgumentException(nameof(maxK));
            int N = y.Length;
            int S = N - 1;

            int usableLevels = 1;
            while (usableLevels < maxK && (S % (1 << usableLevels) == 0)) usableLevels++;

            Complex32[,] R = new Complex32[usableLevels, usableLevels];

            for (int k = 0; k < usableLevels; k++)
            {
                int M = 1 << k;
                int stride = S / M;
                Complex32 h = (b - a) / M;

                Complex32 interiorSum = 0.0;
                for (int m = 1; m < M; m++)
                {
                    int idx = m * stride;
                    interiorSum += y[idx];
                }
                R[k, 0] = 0.5 * h * (y[0] + y[N - 1]) + h * interiorSum;

                Complex32 factor = 4.0;
                for (int j = 1; j <= k; j++, factor *= 4.0)
                    R[k, j] = (factor * R[k, j - 1] - R[k - 1, j - 1]) / (factor - 1.0);

                if (k > 0 && MathF.Abs(R[k, k] - R[k - 1, k - 1]) <= eps * MathF.Abs(R[k, k]))
                    return R[k, k];
            }

            return R[usableLevels - 1, usableLevels - 1];
        }
        #endregion
    }
}
