using System;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines the Zak transform.
    /// </summary>
    [Serializable]
    public class ZakTransform : ITransform
    {
        #region Private data
        private FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        private FastFourierTransform FFT = new FastFourierTransform(false, Direction.Vertical);
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Zak transform.
        /// </summary>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        public ZakTransform(int m)
        {
            M = m;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
        /// </remarks>
        /// </summary>
        public int M { get; set; }
        #endregion

        #region Zak transform
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = A.Length;
            double[] vort = new double[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = A[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = (sum / L2).Real;
            }

            return vort;
        }
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = A.Length;
            Complex[] vort = new Complex[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = A[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = sum / L2;
            }

            return vort;
        }
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Zak transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Zak transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
}
