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
        double sigma_s;
        double sigma_r;
        int iterations;
        double factor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the domain transform filter.
        /// </summary>
        /// <param name="sigma_s">σs</param>
        /// <param name="sigma_r">σr</param>
        /// <param name="iterations">Number of iterations</param>
        /// <param name="factor">Factor [-1, 1]</param>
        public DomainTransformFilter(double sigma_s, double sigma_r, int iterations = 3, double factor = -1.0)
        {
            SigmaS = sigma_s;
            SigmaR = sigma_r;
            Iterations = iterations;
            Factor = factor;
        }
        /// <summary>
        /// Gets or sets the value of σs.
        /// </summary>
        public double SigmaS
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
        public double SigmaR
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
        public double Factor
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
        public void Apply(double[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                double[] copy = (double[])data.Clone();
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
        public void Apply(double[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                double[,] copy = (double[,])data.Clone();
                DomainTransformFilter.domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
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

        #region Private voids
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
        internal static void domainfilter(double[,] I, double sigma_s, double sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            double sigma_H_i;
            int i, j;

            // get differences
            double[,] dIcdx = Matrice.Diff(I, 1, Direction.Horizontal);
            double[,] dIcdy = Matrice.Diff(I, 1, Direction.Vertical);

            // shift patterns
            double[,] dIdx = new double[h, w];
            double[,] dIdy = new double[h, w];

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
                sigma_H_i = sigma_s * Math.Sqrt(3) * Math.Pow(2, (iterations - (i + 1))) / Math.Sqrt(Math.Pow(4, iterations) - 1);

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
        internal static void domainfilter(Complex[,] I, double sigma_s, double sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            double sigma_H_i;
            int i, j;

            // get differences
            Complex[,] dIcdx = Matrice.Diff(I, 1, Direction.Horizontal);
            Complex[,] dIcdy = Matrice.Diff(I, 1, Direction.Vertical);

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
                sigma_H_i = sigma_s * Math.Sqrt(3) * Math.Pow(2, (iterations - (i + 1))) / Math.Sqrt(Math.Pow(4, iterations) - 1);

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
        internal static void domainfilter(double[] I, double sigma_s, double sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            double sigma_H_i;
            int i;

            // get differences
            double[] dIcdy = Matrice.Diff(I, 1);

            // shift patterns
            double[] dIdy = new double[h];

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
                sigma_H_i = sigma_s * Math.Sqrt(3) * Math.Pow(2, (iterations - (i + 1))) / Math.Sqrt(Math.Pow(4, iterations) - 1);

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
        internal static void domainfilter(Complex[] I, double sigma_s, double sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            double sigma_H_i;
            int i;

            // get differences
            Complex[] dIcdy = Matrice.Diff(I, 1);

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
                sigma_H_i = sigma_s * Math.Sqrt(3) * Math.Pow(2, (iterations - (i + 1))) / Math.Sqrt(Math.Pow(4, iterations) - 1);

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
        internal static void tdrf_h(double[,] F, double[,] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
            double[,] V = Matrice.Pow(a, D);
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
        internal static void tdrf_v(double[,] F, double[,] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
            double[,] V = Matrice.Pow(a, D);
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
        internal static void tdrf_h(Complex[,] F, Complex[,] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
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
        internal static void tdrf_v(Complex[,] F, Complex[,] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
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
        internal static void tdrf(double[] F, double[] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
            double[] V = Matrice.Pow(a, D);
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
        internal static void tdrf(Complex[] F, Complex[] D, double sigma)
        {
            // params
            double a = Math.Exp(-Math.Sqrt(2) / sigma);
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
    }
}
