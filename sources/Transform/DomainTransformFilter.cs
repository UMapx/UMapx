using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the domain transform filter.
    /// <remarks>
    /// This filter is a computationally effective analogue of a bilateral filter.
    /// More information can be found on the website:
    /// http://www.inf.ufrgs.br/~eslgastal/DomainTransform/Gastal_Oliveira_SIGGRAPH2011_Domain_Transform.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class DomainTransformFilter : IFilter
    {
        #region Private data
        float sigma_s;
        float sigma_r;
        int iterations;
        float factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the domain transform filter.
        /// </summary>
        /// <param name="sigma_s">σs</param>
        /// <param name="sigma_r">σr</param>
        /// <param name="iterations">Number of iterations</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public DomainTransformFilter(float sigma_s, float sigma_r, int iterations = 3, float factor = -1.0f)
        {
            SigmaS = sigma_s;
            SigmaR = sigma_r;
            Iterations = iterations;
            Factor = factor;
        }
        /// <summary>
        /// Gets or sets the value of σs.
        /// </summary>
        public float SigmaS
        {
            get
            {
                return this.sigma_s;
            }
            set
            {
                this.sigma_s = value;
            }
        }
        /// <summary>
        /// Gets or sets the value of σr.
        /// </summary>
        public float SigmaR
        {
            get
            {
                return this.sigma_r;
            }
            set
            {
                this.sigma_r = value;
            }
        }
        /// <summary>
        /// Gets or sets the number of iterations.
        /// </summary>
        public int Iterations
        {
            get
            {
                return this.iterations;
            }
            set
            {
                this.iterations = value;
            }
        }
        /// <summary>
        /// Gets or sets the factor [-1, 1].
        /// </summary>
        public float Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        #endregion

        #region Public apply voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                float[] copy = (float[])data.Clone();
                DomainTransformFilter.domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0f + this.factor) * (data[i] - copy[i]) + copy[i];
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                float[,] copy = (float[,])data.Clone();
                DomainTransformFilter.domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0f + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                Complex[] copy = (Complex[])data.Clone();
                DomainTransformFilter.domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                Complex[,] copy = (Complex[,])data.Clone();
                DomainTransformFilter.domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }

            return;
        }
        #endregion

        #region Domain private voids
        // **************************************************
        //              DOMAIN TRANSFORM FILTER
        // **************************************************
        // ORIGINALS: Eduardo S.L. Gastal, Manuel M. Oliveira.
        // Domain Transform for Edge-Aware Image and Video 
        // Processing. ACM Transactions on Graphics. 
        // Volume 30 (2011), Number 4.Proceedings of SIGGRAPH 
        // 2011, Article 69.
        // 
        // Designed by Valery Asiryan (c), 2015-2020
        // Moscow, Russia.
        // **************************************************

        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        internal static void domainfilter(float[,] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            float sigma_H_i;
            int i, j;

            // get differences
            float[,] dIcdx = DomainTransformFilter.Diff(I, 1, Direction.Horizontal);
            float[,] dIcdy = DomainTransformFilter.Diff(I, 1, Direction.Vertical);

            // shift patterns
            float[,] dIdx = new float[h, w];
            float[,] dIdy = new float[h, w];

            for (i = 0; i < h; i++)
                for (j = 1; j < w; j++)
                    dIdx[i, j] = Math.Abs(dIcdx[i, j - 1]);

            for (i = 1; i < h; i++)
                for (j = 0; j < w; j++)
                    dIdy[i, j] = Math.Abs(dIcdy[i - 1, j]);

            // sigma patterns and result image
            for (i = 0; i < h; i++)
            {
                for (j = 0; j < w; j++)
                {
                    dIdx[i, j] = 1 + sigma_s / sigma_r * dIdx[i, j];
                    dIdy[i, j] = 1 + sigma_s / sigma_r * dIdy[i, j];
                }
            }

            // iterations
            for (i = 0; i < iterations; i++)
            {
                sigma_H_i = sigma_s * Maths.Sqrt(3) * Maths.Pow(2, (iterations - (i + 1))) / Maths.Sqrt(Maths.Pow(4, iterations) - 1);

                // 2D filter
                tdrf_h(I, dIdx, sigma_H_i);
                tdrf_v(I, dIdy, sigma_H_i);
            }

            return;
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        internal static void domainfilter(Complex[,] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            float sigma_H_i;
            int i, j;

            // get differences
            Complex[,] dIcdx = DomainTransformFilter.Diff(I, 1, Direction.Horizontal);
            Complex[,] dIcdy = DomainTransformFilter.Diff(I, 1, Direction.Vertical);

            // shift patterns
            Complex[,] dIdx = new Complex[h, w];
            Complex[,] dIdy = new Complex[h, w];

            for (i = 0; i < h; i++)
                for (j = 1; j < w; j++)
                    dIdx[i, j] = Maths.Abs(dIcdx[i, j - 1]);

            for (i = 1; i < h; i++)
                for (j = 0; j < w; j++)
                    dIdy[i, j] = Maths.Abs(dIcdy[i - 1, j]);

            // sigma patterns and result image
            for (i = 0; i < h; i++)
            {
                for (j = 0; j < w; j++)
                {
                    dIdx[i, j] = 1 + sigma_s / sigma_r * dIdx[i, j];
                    dIdy[i, j] = 1 + sigma_s / sigma_r * dIdy[i, j];
                }
            }

            // iterations
            for (i = 0; i < iterations; i++)
            {
                sigma_H_i = sigma_s * Maths.Sqrt(3) * Maths.Pow(2, (iterations - (i + 1))) / Maths.Sqrt(Maths.Pow(4, iterations) - 1);

                // 2D filter
                tdrf_h(I, dIdx, sigma_H_i);
                tdrf_v(I, dIdy, sigma_H_i);
            }

            return;
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        internal static void domainfilter(float[] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            float sigma_H_i;
            int i;

            // get differences
            float[] dIcdy = DomainTransformFilter.Diff(I, 1);

            // shift patterns
            float[] dIdy = new float[h];

            for (i = 1; i < h; i++)
                dIdy[i] = Math.Abs(dIcdy[i - 1]);

            // sigma patterns and result image
            for (i = 0; i < h; i++)
            {
                dIdy[i] = 1 + sigma_s / sigma_r * dIdy[i];
            }

            // iterations
            for (i = 0; i < iterations; i++)
            {
                sigma_H_i = sigma_s * Maths.Sqrt(3) * Maths.Pow(2, (iterations - (i + 1))) / Maths.Sqrt(Maths.Pow(4, iterations) - 1);

                // 1D filter
                tdrf(I, dIdy, sigma_H_i);
            }

            return;
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        internal static void domainfilter(Complex[] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            float sigma_H_i;
            int i;

            // get differences
            Complex[] dIcdy = DomainTransformFilter.Diff(I, 1);

            // shift patterns
            Complex[] dIdy = new Complex[h];

            for (i = 1; i < h; i++)
                dIdy[i] = Maths.Abs(dIcdy[i - 1]);

            // sigma patterns and result image
            for (i = 0; i < h; i++)
            {
                dIdy[i] = 1 + sigma_s / sigma_r * dIdy[i];
            }

            // iterations
            for (i = 0; i < iterations; i++)
            {
                sigma_H_i = sigma_s * Maths.Sqrt(3) * Maths.Pow(2, (iterations - (i + 1))) / Maths.Sqrt(Maths.Pow(4, iterations) - 1);

                // 1D filter
                tdrf(I, dIdy, sigma_H_i);
            }

            return;
        }

        /// <summary>
        /// Transformed domain recursive filter (horizontal).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf_h(float[,] F, float[,] D, float sigma)
        {
            // params
            float a = (float)Math.Exp(-Math.Sqrt(2) / sigma);
            float[,] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int w = F.GetLength(1);
            int i, j;

            // Left -> Right filter.
            for (i = 0; i < h; i++)
                for (j = 1; j < w; j++)
                    F[i, j] = F[i, j] + V[i, j] * (F[i, j - 1] - F[i, j]);

            // Right -> Left filter.
            for (i = 0; i < h; i++)
                for (j = w - 2; j >= 0; j--)
                    F[i, j] = F[i, j] + V[i, j + 1] * (F[i, j + 1] - F[i, j]);

            return;
        }
        /// <summary>
        /// Transformed domain recursive filter (vertical).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf_v(float[,] F, float[,] D, float sigma)
        {
            // params
            float a = Maths.Exp(-Maths.Sqrt2 / sigma);
            float[,] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int w = F.GetLength(1);
            int i, j;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                for (j = 0; j < w; j++)
                    F[i, j] = F[i, j] + V[i, j] * (F[i - 1, j] - F[i, j]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                for (j = 0; j < w; j++)
                    F[i, j] = F[i, j] + V[i + 1, j] * (F[i + 1, j] - F[i, j]);

            return;
        }
        /// <summary>
        /// Transformed domain recursive filter (horizontal).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf_h(Complex[,] F, Complex[,] D, float sigma)
        {
            // params
            float a = Maths.Exp(-Maths.Sqrt2 / sigma);
            Complex[,] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int w = F.GetLength(1);
            int i, j;

            // Left -> Right filter.
            for (i = 0; i < h; i++)
                for (j = 1; j < w; j++)
                    F[i, j] = F[i, j] + V[i, j] * (F[i, j - 1] - F[i, j]);

            // Right -> Left filter.
            for (i = 0; i < h; i++)
                for (j = w - 2; j >= 0; j--)
                    F[i, j] = F[i, j] + V[i, j + 1] * (F[i, j + 1] - F[i, j]);

            return;
        }
        /// <summary>
        /// Transformed domain recursive filter (vertical).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf_v(Complex[,] F, Complex[,] D, float sigma)
        {
            // params
            float a = Maths.Exp(-Maths.Sqrt2 / sigma);
            Complex[,] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int w = F.GetLength(1);
            int i, j;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                for (j = 0; j < w; j++)
                    F[i, j] = F[i, j] + V[i, j] * (F[i - 1, j] - F[i, j]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                for (j = 0; j < w; j++)
                    F[i, j] = F[i, j] + V[i + 1, j] * (F[i + 1, j] - F[i, j]);

            return;
        }

        /// <summary>
        /// Transformed domain recursive filter.
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf(float[] F, float[] D, float sigma)
        {
            // params
            float a = Maths.Exp(-Maths.Sqrt2 / sigma);
            float[] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int i;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                F[i] = F[i] + V[i] * (F[i - 1] - F[i]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                F[i] = F[i] + V[i + 1] * (F[i + 1] - F[i]);

            return;
        }
        /// <summary>
        /// Transformed domain recursive filter.
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        internal static void tdrf(Complex[] F, Complex[] D, float sigma)
        {
            // params
            float a = Maths.Exp(-Maths.Sqrt2 / sigma);
            Complex[] V = Matrice.Pow(a, D);
            int h = F.GetLength(0);
            int i;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                F[i] = F[i] + V[i] * (F[i - 1] - F[i]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                F[i] = F[i] + V[i + 1] * (F[i + 1] - F[i]);

            return;
        }
        #endregion

        #region Diff private voids
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Order</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static float[,] Diff(float[,] a, int n, Direction direction, bool reverse = false)
        {
            // start
            int i, r, c;

            // direction of processing
            if (direction == Direction.Horizontal)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    c = a.GetLength(1);

                    if (c == 0)
                        break;

                    a = DiffHorizontal(a, reverse);
                }
            }
            else if (direction == Direction.Vertical)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);

                    if (r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                }
            }
            else
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);
                    c = a.GetLength(1);

                    if (c == 0 || r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                    a = DiffHorizontal(a, reverse);
                }
            }

            return a;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static float[,] DiffVertical(float[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new float[0, m];

            // new array
            float[,] y = new float[r, m];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i - 1, k] - a[i, k];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i, k] - a[i - 1, k];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static float[,] DiffHorizontal(float[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new float[m, 0];

            // new array
            float[,] y = new float[m, c];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i - 1] - a[k, i];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i] - a[k, i - 1];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="n">Order</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static Complex[,] Diff(Complex[,] a, int n, Direction direction, bool reverse = false)
        {
            // start
            int i, r, c;

            // direction of processing
            if (direction == Direction.Horizontal)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    c = a.GetLength(1);

                    if (c == 0)
                        break;

                    a = DiffHorizontal(a, reverse);
                }
            }
            else if (direction == Direction.Vertical)
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);

                    if (r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                }
            }
            else
            {
                // do job
                for (i = 0; i < n; i++)
                {
                    r = a.GetLength(0);
                    c = a.GetLength(1);

                    if (c == 0 || r == 0)
                        break;

                    a = DiffVertical(a, reverse);
                    a = DiffHorizontal(a, reverse);
                }
            }

            return a;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static Complex[,] DiffVertical(Complex[,] a, bool reverse = false)
        {
            // vertical direction 
            // of processing
            int r = a.GetLength(0) - 1;
            int m = a.GetLength(1);
            if (r == 0)
                return new Complex[0, m];

            // new array
            Complex[,] y = new Complex[r, m];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i - 1, k] - a[i, k];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = r; i > 0; i--)
                        y[i - 1, k] = a[i, k] - a[i - 1, k];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of matrix elements.
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Matrix</returns>
        private static Complex[,] DiffHorizontal(Complex[,] a, bool reverse = false)
        {
            // horizontal direction 
            // of processing
            int c = a.GetLength(1) - 1;
            int m = a.GetLength(0);
            if (c == 0)
                return new Complex[m, 0];

            // new array
            Complex[,] y = new Complex[m, c];
            int i, k;

            // do job
            if (reverse)
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i - 1] - a[k, i];
            }
            else
            {
                for (k = 0; k < m; k++)
                    for (i = c; i > 0; i--)
                        y[k, i - 1] = a[k, i] - a[k, i - 1];
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="n">Order</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Array</returns>
        private static float[] Diff(float[] v, int n, bool reverse = false)
        {
            // start
            float[] z;
            float[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new float[0];

                y = new float[length];

                if (reverse)
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i - 1] - z[i];
                }
                else
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i] - z[i - 1];
                }
            }

            return y;
        }
        /// <summary>
        /// Returns the difference of vector elements.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="n">Order</param>
        /// <param name="reverse">Reverse processing or not</param>
        /// <returns>Array</returns>
        private static Complex[] Diff(Complex[] v, int n, bool reverse = false)
        {
            // start
            Complex[] z;
            Complex[] y = v;
            int i, j, length;

            // do job
            for (j = 0; j < n; j++)
            {
                z = y;
                length = z.Length - 1;

                if (length == 0)
                    return new Complex[0];

                y = new Complex[length];

                if (reverse)
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i - 1] - z[i];
                }
                else
                {
                    for (i = length; i > 0; i--)
                        y[i - 1] = z[i] - z[i - 1];
                }
            }

            return y;
        }
        #endregion
    }
}
