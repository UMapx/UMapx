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
                DomainTransformFilter.Domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0f + this.factor) * (data[i] - copy[i]) + copy[i];
            }
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
                DomainTransformFilter.Domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0f + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex32[] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int i;

                // guided filter
                Complex32[] copy = (Complex32[])data.Clone();
                DomainTransformFilter.Domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    data[i] = (1.0 + this.factor) * (data[i] - copy[i]) + copy[i];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex32[,] data)
        {
            // enhancement or not?
            if (this.factor != 0)
            {
                // params
                int l0 = data.GetLength(0);
                int l1 = data.GetLength(1);
                int i, j;

                // guided filter
                Complex32[,] copy = (Complex32[,])data.Clone();
                DomainTransformFilter.Domainfilter(copy, this.sigma_s, this.sigma_r, this.iterations);

                // process
                for (i = 0; i < l0; i++)
                    for (j = 0; j < l1; j++)
                        data[i, j] = (1.0 + this.factor) * (data[i, j] - copy[i, j]) + copy[i, j];
            }
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
        // **************************************************

        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        private static void Domainfilter(float[,] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            float sigma_H_i;
            int i, j;

            // get differences
            float[,] dIcdx = MatrixF.Diff(I, 1, Direction.Horizontal);
            float[,] dIcdy = MatrixF.Diff(I, 1, Direction.Vertical);

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
                sigma_H_i = sigma_s * MathF.Sqrt(3) * MathF.Pow(2, (iterations - (i + 1))) / MathF.Sqrt(MathF.Pow(4, iterations) - 1);

                // 2D filter
                Tdrf_h(I, dIdx, sigma_H_i);
                Tdrf_v(I, dIdy, sigma_H_i);
            }
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        private static void Domainfilter(Complex32[,] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            int w = I.GetLength(1);
            float sigma_H_i;
            int i, j;

            // get differences
            Complex32[,] dIcdx = MatrixF.Diff(I, 1, Direction.Horizontal);
            Complex32[,] dIcdy = MatrixF.Diff(I, 1, Direction.Vertical);

            // shift patterns
            Complex32[,] dIdx = new Complex32[h, w];
            Complex32[,] dIdy = new Complex32[h, w];

            for (i = 0; i < h; i++)
                for (j = 1; j < w; j++)
                    dIdx[i, j] = MathF.Abs(dIcdx[i, j - 1]);

            for (i = 1; i < h; i++)
                for (j = 0; j < w; j++)
                    dIdy[i, j] = MathF.Abs(dIcdy[i - 1, j]);

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
                sigma_H_i = sigma_s * MathF.Sqrt(3) * MathF.Pow(2, (iterations - (i + 1))) / MathF.Sqrt(MathF.Pow(4, iterations) - 1);

                // 2D filter
                Tdrf_h(I, dIdx, sigma_H_i);
                Tdrf_v(I, dIdy, sigma_H_i);
            }
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        private static void Domainfilter(float[] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            float sigma_H_i;
            int i;

            // get differences
            float[] dIcdy = MatrixF.Diff(I, 1);

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
                sigma_H_i = sigma_s * MathF.Sqrt(3) * MathF.Pow(2, (iterations - (i + 1))) / MathF.Sqrt(MathF.Pow(4, iterations) - 1);

                // 1D filter
                Tdrf(I, dIdy, sigma_H_i);
            }
        }
        /// <summary>
        /// Domain transform filter.
        /// </summary>
        /// <param name="I">Input signal</param>
        /// <param name="sigma_s">High sigma</param>
        /// <param name="sigma_r">Low sigma</param>
        /// <param name="iterations">Number of iterations</param>
        private static void Domainfilter(Complex32[] I, float sigma_s, float sigma_r, int iterations = 3)
        {
            // params
            int h = I.GetLength(0);
            float sigma_H_i;
            int i;

            // get differences
            Complex32[] dIcdy = MatrixF.Diff(I, 1);

            // shift patterns
            Complex32[] dIdy = new Complex32[h];

            for (i = 1; i < h; i++)
                dIdy[i] = MathF.Abs(dIcdy[i - 1]);

            // sigma patterns and result image
            for (i = 0; i < h; i++)
            {
                dIdy[i] = 1 + sigma_s / sigma_r * dIdy[i];
            }

            // iterations
            for (i = 0; i < iterations; i++)
            {
                sigma_H_i = sigma_s * MathF.Sqrt(3) * MathF.Pow(2, (iterations - (i + 1))) / MathF.Sqrt(MathF.Pow(4, iterations) - 1);

                // 1D filter
                Tdrf(I, dIdy, sigma_H_i);
            }
        }

        /// <summary>
        /// Transformed domain recursive filter (horizontal).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf_h(float[,] F, float[,] D, float sigma)
        {
            // params
            float a = (float)Math.Exp(-Math.Sqrt(2) / sigma);
            float[,] V = MatrixF.Pow(a, D);
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
        }
        /// <summary>
        /// Transformed domain recursive filter (vertical).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf_v(float[,] F, float[,] D, float sigma)
        {
            // params
            float a = MathF.Exp(-MathF.Sqrt2 / sigma);
            float[,] V = MatrixF.Pow(a, D);
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
        }
        /// <summary>
        /// Transformed domain recursive filter (horizontal).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf_h(Complex32[,] F, Complex32[,] D, float sigma)
        {
            // params
            float a = MathF.Exp(-MathF.Sqrt2 / sigma);
            Complex32[,] V = MatrixF.Pow(a, D);
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
        }
        /// <summary>
        /// Transformed domain recursive filter (vertical).
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf_v(Complex32[,] F, Complex32[,] D, float sigma)
        {
            // params
            float a = MathF.Exp(-MathF.Sqrt2 / sigma);
            Complex32[,] V = MatrixF.Pow(a, D);
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
        }

        /// <summary>
        /// Transformed domain recursive filter.
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf(float[] F, float[] D, float sigma)
        {
            // params
            float a = MathF.Exp(-MathF.Sqrt2 / sigma);
            float[] V = MatrixF.Pow(a, D);
            int h = F.GetLength(0);
            int i;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                F[i] = F[i] + V[i] * (F[i - 1] - F[i]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                F[i] = F[i] + V[i + 1] * (F[i + 1] - F[i]);
        }
        /// <summary>
        /// Transformed domain recursive filter.
        /// </summary>
        /// <param name="F">Input signal</param>
        /// <param name="D">Difference</param>
        /// <param name="sigma">Sigma</param>
        private static void Tdrf(Complex32[] F, Complex32[] D, float sigma)
        {
            // params
            float a = MathF.Exp(-MathF.Sqrt2 / sigma);
            Complex32[] V = MatrixF.Pow(a, D);
            int h = F.GetLength(0);
            int i;

            // Left -> Right filter.
            for (i = 1; i < h; i++)
                F[i] = F[i] + V[i] * (F[i - 1] - F[i]);

            // Right -> Left filter.
            for (i = h - 2; i >= 0; i--)
                F[i] = F[i] + V[i + 1] * (F[i + 1] - F[i]);
        }
        #endregion
    }
}
