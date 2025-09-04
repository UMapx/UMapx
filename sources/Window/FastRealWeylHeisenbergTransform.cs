using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast real Weyl-Heisenberg transform.
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete real orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastRealWeylHeisenbergTransform : RealWeylHeisenbergTransform, IWindowTransform, ITransform
    {
        #region Initialize
        /// <summary>
        /// Initializes fast real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N/4]</param>
        /// <param name="direction">Processing direction</param>
        public FastRealWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical) : base(window, m, direction) { }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N2 = A.Length;
            int N = N2 / 2;

            if (N == this.m)
                throw new ArgumentException("M could not be equal N/2");

            var a = new Complex32[N];

            for (int i = 0; i < N; i++)
            {
                a[i] = new Complex32(A[i], A[i + N]);
            }

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);
            var b = FastWeylHeisenbergTransform.FWHT(a, cache);
            var B = new float[N2];

            for (int i = 0; i < N2; i++)
            {
                B[i] = b[i].Real;
            }

            return B;
        }
        /// <summary>
        /// Backward real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override float[] Backward(float[] B)
        {
            int N2 = B.Length;
            int N = N2 / 2;

            if (N == this.m)
                throw new ArgumentException("M could not be equal N/2");

            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);
            var c = new Complex32[N2];

            for (int i = 0; i < N2; i++)
            {
                c[i] = new Complex32(B[i], 0);
            }

            var a = FastWeylHeisenbergTransform.IFWHT(c, cache);
            var A = new float[N2];
            
            for (int i = 0; i < N; i++)
            {
                A[i + 0] = a[i].Real;
                A[i + N] = a[i].Imag;
            }
            
            return A.Mul(2);
        }
        /// <summary>
        /// Forward real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override float[,] Forward(float[,] A)
        {
            int rows2 = A.GetLength(0);
            int cols2 = A.GetLength(1);
            int rows = rows2 / 2;
            int cols = cols2 / 2;

            if (rows == this.m || cols == this.m)
                throw new ArgumentException("M could not be equal N/2");

            FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
            FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

            if (direction == Direction.Both)
            {
                var tmp = new float[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var col = new Complex32[rows];

                    for (int i = 0; i < rows; i++)
                    {
                        col[i] = new Complex32(A[i, j], A[i + rows, j]);
                    }

                    var tr = FastWeylHeisenbergTransform.FWHT(col, cacheCols);

                    for (int i = 0; i < rows2; i++)
                    {
                        tmp[i, j] = tr[i].Real;
                    }
                });

                var B = new float[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var row = new Complex32[cols];

                    for (int j = 0; j < cols; j++)
                    {
                        row[j] = new Complex32(tmp[i, j], tmp[i, j + cols]);
                    }

                    var tr = FastWeylHeisenbergTransform.FWHT(row, cacheRows);

                    for (int j = 0; j < cols2; j++)
                    {
                        B[i, j] = tr[j].Real;
                    }
                });

                return B;
            }

            if (direction == Direction.Vertical)
            {
                var B = new float[rows2, cols2];
                Parallel.For(0, cols2, j =>
                {
                    var col = new Complex32[rows];

                    for (int i = 0; i < rows; i++)
                    {
                        col[i] = new Complex32(A[i, j], A[i + rows, j]);
                    }

                    var tr = FastWeylHeisenbergTransform.FWHT(col, cacheCols);

                    for (int i = 0; i < rows2; i++)
                    {
                        B[i, j] = tr[i].Real;
                    }
                });
                return B;
            }
            else // Horizontal
            {
                var B = new float[rows2, cols2];
                Parallel.For(0, rows2, i =>
                {
                    var row = new Complex32[cols];

                    for (int j = 0; j < cols; j++)
                    {
                        row[j] = new Complex32(A[i, j], A[i, j + cols]);
                    }
                    
                    var tr = FastWeylHeisenbergTransform.FWHT(row, cacheRows);

                    for (int j = 0; j < cols2; j++)
                    {
                        B[i, j] = tr[j].Real;
                    }
                });
                return B;
            }
        }
        /// <summary>
        /// Backward real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public override float[,] Backward(float[,] B)
        {
            int rows2 = B.GetLength(0);
            int cols2 = B.GetLength(1);
            int rows = rows2 / 2;
            int cols = cols2 / 2;

            if (rows == this.m || cols == this.m)
                throw new ArgumentException("M could not be equal N/2");

            FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
            FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

            if (direction == Direction.Both)
            {
                var tmp = new float[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var col = new Complex32[rows2];

                    for (int i = 0; i < rows2; i++)
                    {
                        col[i] = new Complex32(B[i, j], 0);
                    }

                    var inv = FastWeylHeisenbergTransform.IFWHT(col, cacheCols).Mul(2);

                    for (int i = 0; i < rows; i++)
                    {
                        tmp[i, j] = inv[i].Real;
                        tmp[i + rows, j] = inv[i].Imag;
                    }
                });

                var A = new float[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var row = new Complex32[cols2];

                    for (int j = 0; j < cols2; j++)
                    {
                        row[j] = new Complex32(tmp[i, j], 0);
                    }

                    var inv = FastWeylHeisenbergTransform.IFWHT(row, cacheRows).Mul(2);

                    for (int j = 0; j < cols; j++)
                    {
                        A[i, j] = inv[j].Real;
                        A[i, j + cols] = inv[j].Imag;
                    }
                });

                return A;
            }

            if (direction == Direction.Vertical)
            {
                var A = new float[rows2, cols2];
                Parallel.For(0, cols2, j =>
                {
                    var col = new Complex32[rows2];

                    for (int i = 0; i < rows2; i++)
                    {
                        col[i] = new Complex32(B[i, j], 0);
                    }

                    var inv = FastWeylHeisenbergTransform.IFWHT(col, cacheCols).Mul(2);
                    
                    for (int i = 0; i < rows; i++)
                    {
                        A[i, j] = inv[i].Real;
                        A[i + rows, j] = inv[i].Imag;
                    }
                });
                return A;
            }
            else // Horizontal
            {
                var A = new float[rows2, cols2];
                Parallel.For(0, rows2, i =>
                {
                    var row = new Complex32[cols2];

                    for (int j = 0; j < cols2; j++)
                    {
                        row[j] = new Complex32(B[i, j], 0);
                    }

                    var inv = FastWeylHeisenbergTransform.IFWHT(row, cacheRows).Mul(2);
                    
                    for (int j = 0; j < cols; j++)
                    {
                        A[i, j] = inv[j].Real;
                        A[i, j + cols] = inv[j].Imag;
                    }
                });
                return A;
            }
        }
        /// <summary>
        /// Forward real Weyl-Heisenberg transform for complex input.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Forward(Complex32[] A)
        {
            int N2 = A.Length;
            int N = N2 / 2;

            if (N == this.m)
                throw new ArgumentException("M could not be equal N/2");

            // Build cache for length N (since FWHT(N) -> length 2N)
            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            // Pack real and imaginary parts separately (top half -> Re, bottom half -> Im)
            var packRe = new Complex32[N];
            var packIm = new Complex32[N];

            for (int i = 0; i < N; i++)
            {
                packRe[i] = new Complex32(A[i].Real, A[i + N].Real);
                packIm[i] = new Complex32(A[i].Imag, A[i + N].Imag);
            }

            // Transform both channels
            var trRe = FastWeylHeisenbergTransform.FWHT(packRe, cache); // length = 2N
            var trIm = FastWeylHeisenbergTransform.FWHT(packIm, cache); // length = 2N

            // As in float[] Forward: take only the real part of FWHT outputs;
            // combine (Re-channel, Im-channel) back into Complex32
            var B = new Complex32[N2];

            for (int i = 0; i < N2; i++)
            {
                B[i] = new Complex32(trRe[i].Real, trIm[i].Real);
            }

            return B;
        }
        /// <summary>
        /// Backward real Weyl-Heisenberg transform for complex input.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public override Complex32[] Backward(Complex32[] B)
        {
            int N2 = B.Length;
            int N = N2 / 2;

            if (N == this.m)
                throw new ArgumentException("M could not be equal N/2");

            // Build cache for length N (since IFWHT(2N) -> length N)
            var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

            // Prepare two "real-only" inputs (Re and Im channels)
            var inRe = new Complex32[N2];
            var inIm = new Complex32[N2];

            for (int i = 0; i < N2; i++)
            {
                inRe[i] = new Complex32(B[i].Real, 0.0f);
                inIm[i] = new Complex32(B[i].Imag, 0.0f);
            }

            // Invert both channels
            var invRe = FastWeylHeisenbergTransform.IFWHT(inRe, cache); // length = N
            var invIm = FastWeylHeisenbergTransform.IFWHT(inIm, cache); // length = N

            // Unpack: (top half <- Real parts, bottom half <- Imag parts)
            var A = new Complex32[N2];

            for (int i = 0; i < N; i++)
            {
                A[i] = new Complex32(invRe[i].Real, invIm[i].Real);
                A[i + N] = new Complex32(invRe[i].Imag, invIm[i].Imag);
            }

            // Match the float[] Backward scaling (returns A * 2)
            for (int i = 0; i < N2; i++)
            {
                A[i] = new Complex32(A[i].Real, A[i].Imag);
            }

            return A.Mul(2);
        }
        /// <summary>
        /// Forward real Weyl-Heisenberg transform for complex matrices.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Forward(Complex32[,] A)
        {
            int rows2 = A.GetLength(0);
            int cols2 = A.GetLength(1);
            int rows = rows2 / 2;
            int cols = cols2 / 2;

            if (rows == this.m || cols == this.m)
                throw new ArgumentException("M could not be equal N/2");

            FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
            FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

            if (direction == Direction.Both)
            {
                var tmp = new Complex32[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var colRe = new Complex32[rows];
                    var colIm = new Complex32[rows];

                    for (int i = 0; i < rows; i++)
                    {
                        colRe[i] = new Complex32(A[i, j].Real, A[i + rows, j].Real);
                        colIm[i] = new Complex32(A[i, j].Imag, A[i + rows, j].Imag);
                    }

                    var trRe = FastWeylHeisenbergTransform.FWHT(colRe, cacheCols);
                    var trIm = FastWeylHeisenbergTransform.FWHT(colIm, cacheCols);

                    for (int i = 0; i < rows2; i++)
                    {
                        tmp[i, j] = new Complex32(trRe[i].Real, trIm[i].Real);
                    }
                });

                var B = new Complex32[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var rowRe = new Complex32[cols];
                    var rowIm = new Complex32[cols];

                    for (int j = 0; j < cols; j++)
                    {
                        rowRe[j] = new Complex32(tmp[i, j].Real, tmp[i, j + cols].Real);
                        rowIm[j] = new Complex32(tmp[i, j].Imag, tmp[i, j + cols].Imag);
                    }

                    var trRe = FastWeylHeisenbergTransform.FWHT(rowRe, cacheRows);
                    var trIm = FastWeylHeisenbergTransform.FWHT(rowIm, cacheRows);

                    for (int j = 0; j < cols2; j++)
                    {
                        B[i, j] = new Complex32(trRe[j].Real, trIm[j].Real);
                    }
                });

                return B;
            }

            if (direction == Direction.Vertical)
            {
                var B = new Complex32[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var colRe = new Complex32[rows];
                    var colIm = new Complex32[rows];

                    for (int i = 0; i < rows; i++)
                    {
                        colRe[i] = new Complex32(A[i, j].Real, A[i + rows, j].Real);
                        colIm[i] = new Complex32(A[i, j].Imag, A[i + rows, j].Imag);
                    }

                    var trRe = FastWeylHeisenbergTransform.FWHT(colRe, cacheCols);
                    var trIm = FastWeylHeisenbergTransform.FWHT(colIm, cacheCols);

                    for (int i = 0; i < rows2; i++)
                    {
                        B[i, j] = new Complex32(trRe[i].Real, trIm[i].Real);
                    }
                });

                return B;
            }
            else // Horizontal
            {
                var B = new Complex32[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var rowRe = new Complex32[cols];
                    var rowIm = new Complex32[cols];

                    for (int j = 0; j < cols; j++)
                    {
                        rowRe[j] = new Complex32(A[i, j].Real, A[i, j + cols].Real);
                        rowIm[j] = new Complex32(A[i, j].Imag, A[i, j + cols].Imag);
                    }

                    var trRe = FastWeylHeisenbergTransform.FWHT(rowRe, cacheRows);
                    var trIm = FastWeylHeisenbergTransform.FWHT(rowIm, cacheRows);

                    for (int j = 0; j < cols2; j++)
                    {
                        B[i, j] = new Complex32(trRe[j].Real, trIm[j].Real);
                    }
                });

                return B;
            }
        }
        /// <summary>
        /// Backward real Weyl-Heisenberg transform for complex matrices.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public override Complex32[,] Backward(Complex32[,] B)
        {
            int rows2 = B.GetLength(0);
            int cols2 = B.GetLength(1);
            int rows = rows2 / 2;
            int cols = cols2 / 2;

            if (rows == this.m || cols == this.m)
                throw new ArgumentException("M could not be equal N/2");

            FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
            FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

            if (direction == Direction.Both || direction == Direction.Vertical)
                cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
            if (direction == Direction.Both || direction == Direction.Horizontal)
                cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

            if (direction == Direction.Both)
            {
                var tmp = new Complex32[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var colRe = new Complex32[rows2];
                    var colIm = new Complex32[rows2];

                    for (int i = 0; i < rows2; i++)
                    {
                        colRe[i] = new Complex32(B[i, j].Real, 0);
                        colIm[i] = new Complex32(B[i, j].Imag, 0);
                    }

                    var invRe = FastWeylHeisenbergTransform.IFWHT(colRe, cacheCols).Mul(2);
                    var invIm = FastWeylHeisenbergTransform.IFWHT(colIm, cacheCols).Mul(2);

                    for (int i = 0; i < rows; i++)
                    {
                        tmp[i, j] = new Complex32(invRe[i].Real, invIm[i].Real);
                        tmp[i + rows, j] = new Complex32(invRe[i].Imag, invIm[i].Imag);
                    }
                });

                var A = new Complex32[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var rowRe = new Complex32[cols2];
                    var rowIm = new Complex32[cols2];

                    for (int j = 0; j < cols2; j++)
                    {
                        rowRe[j] = new Complex32(tmp[i, j].Real, 0);
                        rowIm[j] = new Complex32(tmp[i, j].Imag, 0);
                    }

                    var invRe = FastWeylHeisenbergTransform.IFWHT(rowRe, cacheRows).Mul(2);
                    var invIm = FastWeylHeisenbergTransform.IFWHT(rowIm, cacheRows).Mul(2);

                    for (int j = 0; j < cols; j++)
                    {
                        A[i, j] = new Complex32(invRe[j].Real, invIm[j].Real);
                        A[i, j + cols] = new Complex32(invRe[j].Imag, invIm[j].Imag);
                    }
                });

                return A;
            }

            if (direction == Direction.Vertical)
            {
                var A = new Complex32[rows2, cols2];

                Parallel.For(0, cols2, j =>
                {
                    var colRe = new Complex32[rows2];
                    var colIm = new Complex32[rows2];

                    for (int i = 0; i < rows2; i++)
                    {
                        colRe[i] = new Complex32(B[i, j].Real, 0);
                        colIm[i] = new Complex32(B[i, j].Imag, 0);
                    }

                    var invRe = FastWeylHeisenbergTransform.IFWHT(colRe, cacheCols).Mul(2);
                    var invIm = FastWeylHeisenbergTransform.IFWHT(colIm, cacheCols).Mul(2);

                    for (int i = 0; i < rows; i++)
                    {
                        A[i, j] = new Complex32(invRe[i].Real, invIm[i].Real);
                        A[i + rows, j] = new Complex32(invRe[i].Imag, invIm[i].Imag);
                    }
                });

                return A;
            }
            else // Horizontal
            {
                var A = new Complex32[rows2, cols2];

                Parallel.For(0, rows2, i =>
                {
                    var rowRe = new Complex32[cols2];
                    var rowIm = new Complex32[cols2];

                    for (int j = 0; j < cols2; j++)
                    {
                        rowRe[j] = new Complex32(B[i, j].Real, 0);
                        rowIm[j] = new Complex32(B[i, j].Imag, 0);
                    }

                    var invRe = FastWeylHeisenbergTransform.IFWHT(rowRe, cacheRows).Mul(2);
                    var invIm = FastWeylHeisenbergTransform.IFWHT(rowIm, cacheRows).Mul(2);

                    for (int j = 0; j < cols; j++)
                    {
                        A[i, j] = new Complex32(invRe[j].Real, invIm[j].Real);
                        A[i, j + cols] = new Complex32(invRe[j].Imag, invIm[j].Imag);
                    }
                });

                return A;
            }
        }
        #endregion
    }
}