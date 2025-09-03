using System;
using UMapx.Core;

namespace UMapx.Analysis
{
    /// <summary>
    /// Defines a class that implements interpolation.
    /// <remarks>
    /// This class is a solution to the problem of finding an intermediate value of the function F(x).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Interpolation
    {
        #region Private data
        private InterpolationMethod method;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes a class that implements interpolation.
        /// </summary>
        /// <param name="method">Interpolation method</param>
        public Interpolation(InterpolationMethod method = InterpolationMethod.Lagrange)
        {
            this.method = method;
        }
        /// <summary>
        /// Gets or sets the interpolation method.
        /// </summary>
        public InterpolationMethod MethodType
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
        /// Returns the value of a function at a point.
        /// <remarks>
        /// In this case, only bilinear interpolation is used.
        /// </remarks>
        /// </summary>
        /// <param name="x">Array of values of the first argument</param>
        /// <param name="y">Array of values of the second argument</param>
        /// <param name="z">Function matrix</param>
        /// <param name="xl">The value of the first argument to calculate</param>
        /// <param name="yl">The value of the second argument to calculate</param>
        /// <returns>Value</returns>
        public float Compute(float[] x, float[] y, float[,] z, float xl, float yl)
        {
            return Bilinear(x, y, z, xl, yl);
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Value</returns>
        public float Compute(float[] x, float[] y, float xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case InterpolationMethod.Newton:
                    return Interpolation.Newto(x, y, xl);

                case InterpolationMethod.Barycentric:
                    return Interpolation.Baryc(x, y, xl);

                case InterpolationMethod.Lagrange:
                    return Interpolation.Lagra(x, y, xl);

                case InterpolationMethod.Linear:
                default:
                    return Interpolation.Linear(x, y, xl);
            }
        }
        /// <summary>
        /// Returns the value of a function at a point.
        /// </summary>
        /// <param name="x">Array of values of the argument</param>
        /// <param name="y">Array of values of the function</param>
        /// <param name="xl">The value of the argument to calculate</param>
        /// <returns>Complex number</returns>
        public ComplexF Compute(ComplexF[] x, ComplexF[] y, ComplexF xl)
        {
            // chose method of interpolation
            switch (method)
            {
                case InterpolationMethod.Newton:
                    return Interpolation.Newto(x, y, xl);

                case InterpolationMethod.Barycentric:
                    return Interpolation.Baryc(x, y, xl);

                case InterpolationMethod.Lagrange:
                default:
                    return Interpolation.Lagra(x, y, xl);

                case InterpolationMethod.Linear:
                    throw new NotSupportedException("Linear interpolation is not defined for complex-valued functions");
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Piecewise-linear interpolation on a sorted 1D grid using binary search.
        /// </summary>
        /// <remarks>
        /// - Expects strictly increasing <paramref name="x"/> with the same length as <paramref name="y"/>.<br/>
        /// - Returns endpoint values for out-of-range <paramref name="xl"/> (clamped extrapolation).<br/>
        /// - Time complexity: O(log n) due to binary search.
        /// </remarks>
        /// <param name="x">Sorted grid nodes x[0..n-1], strictly increasing</param>
        /// <param name="y">Function samples y[i] = f(x[i]) of the same length as x</param>
        /// <param name="xl">Query point</param>
        /// <returns>Interpolated value at xl (clamped to the nearest endpoint if outside [x0, x_{n-1}])</returns>
        /// <exception cref="ArgumentException">Thrown if arrays are null, lengths mismatch, or empty</exception>
        private static float Linear(float[] x, float[] y, float xl)
        {
            int n = x?.Length ?? 0;
            if (n == 0 || y == null || y.Length != n) throw new ArgumentException();
            if (xl <= x[0]) return y[0];
            if (xl >= x[n - 1]) return y[n - 1];

            int lo = 0, hi = n - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) >> 1;
                if (xl >= x[mid]) lo = mid; else hi = mid;
            }
            float t = (xl - x[lo]) / (x[lo + 1] - x[lo]);
            return y[lo] + t * (y[lo + 1] - y[lo]);
        }
        /// <summary>
        /// Bilinear interpolation on a 2D rectilinear grid with clamped behavior at the edges.
        /// </summary>
        /// <remarks>
        /// - Expects strictly increasing 1D grids <paramref name="x"/> (size nx) and <paramref name="y"/> (size ny).<br/>
        /// - <paramref name="z"/> must be an nx-by-ny matrix with z[i,j] = f(x[i], y[j]).<br/>
        /// - If (<paramref name="xval"/>, <paramref name="yval"/>) is outside the grid rectangle, the value from the nearest
        ///   boundary cell is returned (clamped).<br/>
        /// - Time complexity: O(log nx + log ny) due to binary searches.
        /// </remarks>
        /// <param name="x">X-grid nodes (length nx), strictly increasing</param>
        /// <param name="y">Y-grid nodes (length ny), strictly increasing</param>
        /// <param name="z">Function values, shape [nx, ny]</param>
        /// <param name="xval">Query x-coordinate</param>
        /// <param name="yval">Query y-coordinate</param>
        /// <returns>Interpolated value at (xval, yval)</returns>
        /// <exception cref="ArgumentException">Thrown if grid sizes are invalid or inconsistent</exception>

        private static float Bilinear(float[] x, float[] y, float[,] z, float xval, float yval)
        {
            int nx = x?.Length ?? 0, ny = y?.Length ?? 0;
            if (nx < 2 || ny < 2 || z == null || z.GetLength(0) != nx || z.GetLength(1) != ny)
                throw new ArgumentException();

            if (xval <= x[0]) return (yval <= y[0]) ? z[0, 0] : (yval >= y[ny - 1] ? z[0, ny - 1] : z[0, LowerIndex(y, yval)]);
            if (xval >= x[nx - 1]) return (yval <= y[0]) ? z[nx - 1, 0] : (yval >= y[ny - 1] ? z[nx - 1, ny - 1] : z[nx - 1, LowerIndex(y, yval)]);

            int ix = LowerIndex(x, xval); // x[ix] <= xval < x[ix+1]
            int iy = LowerIndex(y, yval); // y[iy] <= yval < y[iy+1]

            float x0 = x[ix], x1 = x[ix + 1];
            float y0 = y[iy], y1 = y[iy + 1];
            float tx = (xval - x0) / (x1 - x0);
            float ty = (yval - y0) / (y1 - y0);

            float z00 = z[ix, iy], z10 = z[ix + 1, iy], z01 = z[ix, iy + 1], z11 = z[ix + 1, iy + 1];
            float z0 = z00 + tx * (z10 - z00);
            float z1 = z01 + tx * (z11 - z01);
            return z0 + ty * (z1 - z0);

            static int LowerIndex(float[] a, float v)
            {
                int lo = 0, hi = a.Length - 1;
                if (v <= a[0]) return 0;
                if (v >= a[hi]) return hi - 1;
                while (hi - lo > 1)
                {
                    int mid = (lo + hi) >> 1;
                    if (v >= a[mid]) lo = mid; else hi = mid;
                }
                return lo;
            }
        }
        /// <summary>
        /// Lagrange polynomial interpolation (naive O(n²) evaluation).
        /// </summary>
        /// <remarks>
        /// - Expects pairwise distinct nodes <paramref name="x"/> (not necessarily uniform).<br/>
        /// - Numerically unstable for large n; prefer barycentric form for better stability.
        /// </remarks>
        /// <param name="x">Interpolation nodes x[0..n-1]</param>
        /// <param name="y">Function samples y[i] = f(x[i])</param>
        /// <param name="xval">Query point</param>
        /// <returns>Interpolated value at xval</returns>
        private static float Lagra(float[] x, float[] y, float xval)
        {
            float yval = 0.0f;
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                float products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += products;
            }
            return yval;
        }
        /// <summary>
        /// Newton interpolation in divided differences with Horner-like evaluation.
        /// </summary>
        /// <remarks>
        /// - Builds the divided-difference table in-place (O(n²)), then evaluates in O(n).<br/>
        /// - More numerically stable than the naïve Lagrange form; nodes need not be uniform.
        /// </remarks>
        /// <param name="x">Interpolation nodes x[0..n-1] (distinct)</param>
        /// <param name="y">Function samples y[i] = f(x[i])</param>
        /// <param name="xval">Query point</param>
        /// <returns>Interpolated value at xval</returns>
        private static float Newto(float[] x, float[] y, float xval)
        {
            int n = x.Length;
            var a = (float[])y.Clone();
            for (int i = 1; i < n; i++)
                for (int j = n - 1; j >= i; j--)
                    a[j] = (a[j] - a[j - 1]) / (x[j] - x[j - i]);

            float res = a[n - 1];
            for (int i = n - 2; i >= 0; i--)
                res = a[i] + (xval - x[i]) * res;
            return res;
        }
        /// <summary>
        /// First-form barycentric Lagrange interpolation (precompute weights, O(n²); evaluate O(n)).
        /// </summary>
        /// <remarks>
        /// - Requires pairwise distinct nodes; throws if duplicates are detected.<br/>
        /// - More numerically robust than naïve Lagrange; for repeated queries, cache weights <c>w[i]</c>.
        /// </remarks>
        /// <param name="x">Interpolation nodes (distinct)</param>
        /// <param name="y">Function samples at nodes</param>
        /// <param name="xval">Query point; if equal to a node, returns the corresponding sample exactly</param>
        /// <returns>Interpolated value at xval</returns>
        /// <exception cref="ArgumentException">Thrown when duplicate nodes are detected</exception>
        private static float Baryc(float[] x, float[] y, float xval)
        {
            int n = x.Length;
            for (int i = 0; i < n; i++) if (xval == x[i]) return y[i];

            var w = new float[n];
            for (int i = 0; i < n; i++)
            {
                float prod = 1f;
                for (int j = 0; j < n; j++) if (i != j) prod *= (x[i] - x[j]);
                if (prod == 0f) throw new ArgumentException("Duplicate nodes");
                w[i] = 1f / prod;
            }
            float num = 0f, den = 0f;
            for (int i = 0; i < n; i++)
            {
                float d = xval - x[i];
                float wi = w[i] / d;
                num += wi * y[i];
                den += wi;
            }
            return num / den;
        }
        /// <summary>
        /// Lagrange polynomial interpolation (naive O(n²) evaluation).
        /// </summary>
        /// <remarks>
        /// - Expects pairwise distinct nodes <paramref name="x"/> (not necessarily uniform).<br/>
        /// - Numerically unstable for large n; prefer barycentric form for better stability.
        /// </remarks>
        /// <param name="x">Interpolation nodes x[0..n-1]</param>
        /// <param name="y">Function samples y[i] = f(x[i])</param>
        /// <param name="xval">Query point</param>
        /// <returns>Interpolated value at xval</returns>
        private static ComplexF Lagra(ComplexF[] x, ComplexF[] y, ComplexF xval)
        {
            ComplexF yval = 0.0;
            int length = x.Length;
            int i, j;

            for (i = 0; i < length; i++)
            {
                ComplexF products = y[i];
                for (j = 0; j < length; j++)
                {
                    if (i != j)
                    {
                        products *= (xval - x[j]) / (x[i] - x[j]);
                    }
                }
                yval += products;
            }
            return yval;
        }
        /// <summary>
        /// Newton interpolation in divided differences with Horner-like evaluation.
        /// </summary>
        /// <remarks>
        /// - Builds the divided-difference table in-place (O(n²)), then evaluates in O(n).<br/>
        /// - More numerically stable than the naïve Lagrange form; nodes need not be uniform.
        /// </remarks>
        /// <param name="x">Interpolation nodes x[0..n-1] (distinct)</param>
        /// <param name="y">Function samples y[i] = f(x[i])</param>
        /// <param name="xval">Query point</param>
        /// <returns>Interpolated value at xval</returns>
        private static ComplexF Newto(ComplexF[] x, ComplexF[] y, ComplexF xval)
        {
            int n = x.Length;
            var a = (ComplexF[])y.Clone();
            for (int i = 1; i < n; i++)
                for (int j = n - 1; j >= i; j--)
                    a[j] = (a[j] - a[j - 1]) / (x[j] - x[j - i]);

            ComplexF res = a[n - 1];
            for (int i = n - 2; i >= 0; i--)
                res = a[i] + (xval - x[i]) * res;
            return res;
        }
        /// <summary>
        /// First-form barycentric Lagrange interpolation (precompute weights, O(n²); evaluate O(n)).
        /// </summary>
        /// <remarks>
        /// - Requires pairwise distinct nodes; throws if duplicates are detected.<br/>
        /// - More numerically robust than naïve Lagrange; for repeated queries, cache weights <c>w[i]</c>.
        /// </remarks>
        /// <param name="x">Interpolation nodes (distinct)</param>
        /// <param name="y">Function samples at nodes</param>
        /// <param name="xval">Query point; if equal to a node, returns the corresponding sample exactly</param>
        /// <returns>Interpolated value at xval</returns>
        /// <exception cref="ArgumentException">Thrown when duplicate nodes are detected</exception>
        private static ComplexF Baryc(ComplexF[] x, ComplexF[] y, ComplexF xval)
        {
            int n = x.Length;
            for (int i = 0; i < n; i++) if (xval == x[i]) return y[i];

            var w = new ComplexF[n];
            for (int i = 0; i < n; i++)
            {
                ComplexF prod = 1f;
                for (int j = 0; j < n; j++) if (i != j) prod *= (x[i] - x[j]);
                if (prod == 0f) throw new ArgumentException("Duplicate nodes");
                w[i] = 1f / prod;
            }
            ComplexF num = 0f, den = 0f;
            for (int i = 0; i < n; i++)
            {
                ComplexF d = xval - x[i];
                ComplexF wi = w[i] / d;
                num += wi * y[i];
                den += wi;
            }
            return num / den;
        }
        #endregion
    }
}
