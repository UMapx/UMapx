using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast Weyl-Heisenberg transform.
    /// </summary>
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    [Serializable]
    public partial class FastWeylHeisenbergTransform : WeylHeisenbergTransform, IWindowTransform, ITransform
    {
        // Go to "FastWeylHeisenbergTransform.Asiryan.cs" file to find out more about the fast O(N log N) implementation
        // designed and developed by Valery Asiryan, Yerevan, Armenia (2025)

        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N/2]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N = A.Length;
            var cache = PolyphaseCache.Build(N, this.m, this.window);
            return FWHT(A, cache);
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N = B.Length / 2;
            var cache = PolyphaseCache.Build(N, this.m, this.window);
            return IFWHT(B, cache);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Forward(Complex32[,] A)
        {
            int N = A.GetLength(0); // rows
            int M = A.GetLength(1); // cols

            PolyphaseCache cacheCols = null; // U for length N (vertical)
            PolyphaseCache cacheRows = null; // V for length M (horizontal)

            if (Direction == Direction.Both || Direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N, this.m, this.window);

            if (Direction == Direction.Both || Direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M, this.m, this.window);

            if (Direction == Direction.Both)
            {
                // 1) Left-multiply: U^H · A  → shape: (2N × M)
                var tmp = new Complex32[2 * N, M];

                Parallel.For(0, M, j =>
                {
                    var col = new Complex32[N];

                    for (int i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    var tr = FWHT(col, cacheCols); // length 2N

                    for (int i = 0; i < 2 * N; i++)
                    {
                        tmp[i, j] = tr[i];
                    }
                });

                // 2) Right-multiply by V via conjugation trick:
                //    row * V = conj( FWHT( conj(row) ) )  → shape: (2N × 2M)
                var B = new Complex32[2 * N, 2 * M];

                Parallel.For(0, 2 * N, i =>
                {
                    var row = new Complex32[M];

                    for (int j = 0; j < M; j++)
                    {
                        row[j] = tmp[i, j];
                    }

                    // conj → FWHT → conj
                    for (int j = 0; j < M; j++)
                    {
                        row[j] = new Complex32(row[j].Real, -row[j].Imag);
                    }

                    var tr = FWHT(row, cacheRows); // V^H * conj(row)

                    for (int j = 0; j < 2 * M; j++)
                    {
                        var c = tr[j];
                        B[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return B;
            }

            if (Direction == Direction.Vertical)
            {
                // B = U^H · A  → (2N × M)
                var B = new Complex32[2 * N, M];

                Parallel.For(0, M, j =>
                {
                    var col = new Complex32[N];

                    for (int i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    var tr = FWHT(col, cacheCols); // 2N

                    for (int i = 0; i < 2 * N; i++)
                    {
                        B[i, j] = tr[i];
                    }
                });

                return B;
            }
            else // Horizontal
            {
                // B = A · V  → (N × 2M)  (via the same conjugation trick per row)
                var B = new Complex32[N, 2 * M];

                Parallel.For(0, N, i =>
                {
                    var row = new Complex32[M];

                    for (int j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    for (int j = 0; j < M; j++)
                    {
                        row[j] = new Complex32(row[j].Real, -row[j].Imag);
                    }

                    var tr = FWHT(row, cacheRows); // V^H * conj(row)

                    for (int j = 0; j < 2 * M; j++)
                    {
                        var c = tr[j];
                        B[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return B;
            }
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Backward(Complex32[,] B)
        {
            int N2 = B.GetLength(0); // possibly 2N
            int M2 = B.GetLength(1); // possibly 2M

            PolyphaseCache cacheCols = null; // for target N = N2/2 (vertical inverse)
            PolyphaseCache cacheRows = null; // for target M = M2/2 (horizontal inverse)

            if (Direction == Direction.Both || Direction == Direction.Vertical)
                cacheCols = PolyphaseCache.Build(N2 / 2, this.m, this.window);

            if (Direction == Direction.Both || Direction == Direction.Horizontal)
                cacheRows = PolyphaseCache.Build(M2 / 2, this.m, this.window);

            if (Direction == Direction.Both)
            {
                // Inverse order of Forward(Both):
                // 1) Undo right multiply (· V) per row using the conjugation trick:
                //    row = conj( IFWHT( conj(B_row) ) )  → shape: (N2 × M2/2)
                var tmp = new Complex32[N2, M2 / 2];

                Parallel.For(0, N2, i =>
                {
                    var row = new Complex32[M2];

                    for (int j = 0; j < M2; j++)
                    {
                        row[j] = B[i, j];
                    }

                    for (int j = 0; j < M2; j++)
                    {
                        row[j] = new Complex32(row[j].Real, -row[j].Imag); // conj(B_row)
                    }

                    var inv = IFWHT(row, cacheRows); // length M2/2

                    for (int j = 0; j < M2 / 2; j++)
                    {
                        var c = inv[j];
                        tmp[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                // 2) Undo left multiply (U^H ·) per column: 2N → N
                var A = new Complex32[N2 / 2, M2 / 2];

                Parallel.For(0, M2 / 2, j =>
                {
                    var col = new Complex32[N2];

                    for (int i = 0; i < N2; i++)
                    {
                        col[i] = tmp[i, j];
                    }

                    var inv = IFWHT(col, cacheCols); // length N2/2

                    for (int i = 0; i < N2 / 2; i++)
                    {
                        A[i, j] = inv[i];
                    }
                });

                return A;
            }

            if (Direction == Direction.Vertical)
            {
                // A = (U^H)^{-1} · B : columns 2N → N
                var A = new Complex32[N2 / 2, M2];

                Parallel.For(0, M2, j =>
                {
                    var col = new Complex32[N2];

                    for (int i = 0; i < N2; i++)
                    {
                        col[i] = B[i, j];
                    }

                    var inv = IFWHT(col, cacheCols); // N2/2
                    
                    for (int i = 0; i < N2 / 2; i++)
                    {
                        A[i, j] = inv[i];
                    }
                });

                return A;
            }
            else // Horizontal
            {
                // A = B · V^{-1} : rows 2M → M (via conjugation trick)
                var A = new Complex32[N2, M2 / 2];

                Parallel.For(0, N2, i =>
                {
                    var row = new Complex32[M2];

                    for (int j = 0; j < M2; j++)
                    {
                        row[j] = B[i, j];
                    }

                    for (int j = 0; j < M2; j++)
                    {
                        row[j] = new Complex32(row[j].Real, -row[j].Imag); // conj(B_row)
                    }

                    var inv = IFWHT(row, cacheRows); // M2/2

                    for (int j = 0; j < M2 / 2; j++)
                    {
                        var c = inv[j];
                        A[i, j] = new Complex32(c.Real, -c.Imag); // conj(...)
                    }
                });

                return A;
            }
        }

        #endregion
    }
}
