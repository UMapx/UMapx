using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    /// <summary>
    /// Defines fast real Weyl-Heisenberg transform.
    /// </summary>
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete real orthogonal
    /// Weyl-Heisenberg transforms. This implementation was designed and developed by Valery Asiryan, Yerevan, Armenia (2025).
    /// </remarks>
    [Serializable]
    public partial class FastRealWeylHeisenbergTransform : RealWeylHeisenbergTransform, IWindowTransform, ITransform
    {
        // Go to "FastWeylHeisenbergTransform.Fourier.cs" and "FastRealWeylHeisenbergTransform.Hartley.cs" files
        // to find out more about the fast O(N log N) implementations.

        #region Initialize
        /// <summary>
        /// Initializes fast real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [2, N/2]</param>
        /// <param name="spectrumType">Spectrum type</param>
        /// <param name="direction">Processing direction</param>
        public FastRealWeylHeisenbergTransform(IWindow window, int m = 8, SpectrumType spectrumType = SpectrumType.Fourier, Direction direction = Direction.Vertical) 
            : base(window, m, spectrumType, direction) { }
        #endregion

        #region Real Weyl-Heisenberg Transform
        /// <summary>
        /// Forward real Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public override float[] Forward(float[] A)
        {
            int N2 = A.Length;
            int N = N2 / 2;

            if (SpectrumType == SpectrumType.Fourier)
            {
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
            else
            {
                var cache = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(N, this.m, this.window);
                var B = FastRealWeylHeisenbergTransform.FRWHT(A, cache);
                return B;
            }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
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
            else
            {
                var cache = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(N, this.m, this.window);
                var A = FastRealWeylHeisenbergTransform.IFRWHT(B, cache);
                return A;
            }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
                FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new float[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var col = new Complex32[rows];
                        for (int i = 0; i < rows; i++)
                            col[i] = new Complex32(A[i, j], A[i + rows, j]);

                        var tr = FastWeylHeisenbergTransform.FWHT(col, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            tmp[i, j] = tr[i].Real;
                    });

                    var B = new float[rows2, cols2];

                    Parallel.For(0, rows2, i =>
                    {
                        var row = new Complex32[cols];
                        for (int j = 0; j < cols; j++)
                            row[j] = new Complex32(tmp[i, j], tmp[i, j + cols]);

                        var tr = FastWeylHeisenbergTransform.FWHT(row, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = tr[j].Real;
                    });

                    return B;
                }

                if (Direction == Direction.Vertical)
                {
                    var B = new float[rows2, cols2];
                    Parallel.For(0, cols2, j =>
                    {
                        var col = new Complex32[rows];
                        for (int i = 0; i < rows; i++)
                            col[i] = new Complex32(A[i, j], A[i + rows, j]);

                        var tr = FastWeylHeisenbergTransform.FWHT(col, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            B[i, j] = tr[i].Real;
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
                            row[j] = new Complex32(A[i, j], A[i, j + cols]);

                        var tr = FastWeylHeisenbergTransform.FWHT(row, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = tr[j].Real;
                    });
                    return B;
                }
            }
            else // SpectrumType.Hartley
            {
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheCols = null;
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new float[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var y = new float[rows2]; // [top; bottom]
                        for (int i = 0; i < rows; i++)
                        {
                            y[i] = A[i, j];
                            y[i + rows] = A[i + rows, j];
                        }

                        var tr = FastRealWeylHeisenbergTransform.FRWHT(y, cacheCols); // len = rows2
                        for (int i = 0; i < rows2; i++)
                            tmp[i, j] = tr[i];
                    });

                    var B = new float[rows2, cols2];
                    Parallel.For(0, rows2, i =>
                    {
                        var y = new float[cols2];
                        for (int j = 0; j < cols; j++)
                        {
                            y[j] = tmp[i, j];
                            y[j + cols] = tmp[i, j + cols];
                        }

                        var tr = FastRealWeylHeisenbergTransform.FRWHT(y, cacheRows); // len = cols2
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = tr[j];
                    });

                    return B;
                }

                if (Direction == Direction.Vertical)
                {
                    var B = new float[rows2, cols2];
                    Parallel.For(0, cols2, j =>
                    {
                        var y = new float[rows2];
                        for (int i = 0; i < rows; i++)
                        {
                            y[i] = A[i, j];
                            y[i + rows] = A[i + rows, j];
                        }

                        var tr = FastRealWeylHeisenbergTransform.FRWHT(y, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            B[i, j] = tr[i];
                    });
                    return B;
                }
                else // Horizontal
                {
                    var B = new float[rows2, cols2];
                    Parallel.For(0, rows2, i =>
                    {
                        var y = new float[cols2];
                        for (int j = 0; j < cols; j++)
                        {
                            y[j] = A[i, j];
                            y[j + cols] = A[i, j + cols];
                        }

                        var tr = FastRealWeylHeisenbergTransform.FRWHT(y, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = tr[j];
                    });
                    return B;
                }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
                FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new float[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var col = new Complex32[rows2];
                        for (int i = 0; i < rows2; i++)
                            col[i] = new Complex32(B[i, j], 0);

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
                            row[j] = new Complex32(tmp[i, j], 0);

                        var inv = FastWeylHeisenbergTransform.IFWHT(row, cacheRows).Mul(2);
                        for (int j = 0; j < cols; j++)
                        {
                            A[i, j] = inv[j].Real;
                            A[i, j + cols] = inv[j].Imag;
                        }
                    });

                    return A;
                }

                if (Direction == Direction.Vertical)
                {
                    var A = new float[rows2, cols2];
                    Parallel.For(0, cols2, j =>
                    {
                        var col = new Complex32[rows2];
                        for (int i = 0; i < rows2; i++)
                            col[i] = new Complex32(B[i, j], 0);

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
                            row[j] = new Complex32(B[i, j], 0);

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
            else // SpectrumType.Hartley
            {
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheCols = null;
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new float[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var c = new float[rows2];
                        for (int i = 0; i < rows2; i++)
                            c[i] = B[i, j];

                        var y = FastRealWeylHeisenbergTransform.IFRWHT(c, cacheCols); // len = rows2
                        for (int i = 0; i < rows2; i++)
                            tmp[i, j] = y[i];
                    });

                    var A = new float[rows2, cols2];

                    Parallel.For(0, rows2, i =>
                    {
                        var c = new float[cols2];
                        for (int j = 0; j < cols2; j++)
                            c[j] = tmp[i, j];

                        var y = FastRealWeylHeisenbergTransform.IFRWHT(c, cacheRows); // len = cols2
                        for (int j = 0; j < cols2; j++)
                            A[i, j] = y[j];
                    });

                    return A;
                }

                if (Direction == Direction.Vertical)
                {
                    var A = new float[rows2, cols2];
                    Parallel.For(0, cols2, j =>
                    {
                        var c = new float[rows2];
                        for (int i = 0; i < rows2; i++)
                            c[i] = B[i, j];

                        var y = FastRealWeylHeisenbergTransform.IFRWHT(c, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            A[i, j] = y[i];
                    });
                    return A;
                }
                else // Horizontal
                {
                    var A = new float[rows2, cols2];
                    Parallel.For(0, rows2, i =>
                    {
                        var c = new float[cols2];
                        for (int j = 0; j < cols2; j++)
                            c[j] = B[i, j];

                        var y = FastRealWeylHeisenbergTransform.IFRWHT(c, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            A[i, j] = y[j];
                    });
                    return A;
                }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

                var packRe = new Complex32[N];
                var packIm = new Complex32[N];
                for (int i = 0; i < N; i++)
                {
                    packRe[i] = new Complex32(A[i].Real, A[i + N].Real);
                    packIm[i] = new Complex32(A[i].Imag, A[i + N].Imag);
                }

                var trRe = FastWeylHeisenbergTransform.FWHT(packRe, cache);
                var trIm = FastWeylHeisenbergTransform.FWHT(packIm, cache);

                var B = new Complex32[N2];
                for (int i = 0; i < N2; i++)
                    B[i] = new Complex32(trRe[i].Real, trIm[i].Real);

                return B;
            }
            else // SpectrumType.Hartley
            {
                var cache = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(N, this.m, this.window);

                var yRe = new float[N2];
                for (int i = 0; i < N; i++)
                {
                    yRe[i] = A[i].Real;
                    yRe[i + N] = A[i + N].Real;
                }
                var cRe = FastRealWeylHeisenbergTransform.FRWHT(yRe, cache);

                var yIm = new float[N2];
                for (int i = 0; i < N; i++)
                {
                    yIm[i] = A[i].Imag;
                    yIm[i + N] = A[i + N].Imag;
                }
                var cIm = FastRealWeylHeisenbergTransform.FRWHT(yIm, cache);

                var B = new Complex32[N2];
                for (int i = 0; i < N2; i++)
                    B[i] = new Complex32(cRe[i], cIm[i]);

                return B;
            }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                var cache = FastWeylHeisenbergTransform.PolyphaseCache.Build(N, this.m, this.window);

                var inRe = new Complex32[N2];
                var inIm = new Complex32[N2];
                for (int i = 0; i < N2; i++)
                {
                    inRe[i] = new Complex32(B[i].Real, 0.0f);
                    inIm[i] = new Complex32(B[i].Imag, 0.0f);
                }

                var invRe = FastWeylHeisenbergTransform.IFWHT(inRe, cache); // len = N
                var invIm = FastWeylHeisenbergTransform.IFWHT(inIm, cache); // len = N

                var A = new Complex32[N2];
                for (int i = 0; i < N; i++)
                {
                    A[i] = new Complex32(invRe[i].Real, invIm[i].Real);
                    A[i + N] = new Complex32(invRe[i].Imag, invIm[i].Imag);
                }
                return A.Mul(2);
            }
            else // SpectrumType.Hartley
            {
                var cache = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(N, this.m, this.window);

                var cRe = new float[N2];
                for (int i = 0; i < N2; i++) cRe[i] = B[i].Real;
                var yRe = FastRealWeylHeisenbergTransform.IFRWHT(cRe, cache); // len = N2

                var cIm = new float[N2];
                for (int i = 0; i < N2; i++) cIm[i] = B[i].Imag;
                var yIm = FastRealWeylHeisenbergTransform.IFRWHT(cIm, cache); // len = N2

                var A = new Complex32[N2];
                for (int i = 0; i < N; i++)
                {
                    A[i] = new Complex32(yRe[i], yIm[i]);
                    A[i + N] = new Complex32(yRe[i + N], yIm[i + N]);
                }
                return A;
            }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
                FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
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

                if (Direction == Direction.Vertical)
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
            else // SpectrumType.Hartley
            {
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheCols = null;
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new Complex32[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var yRe = new float[rows2];
                        var yIm = new float[rows2];
                        for (int i = 0; i < rows; i++)
                        {
                            yRe[i] = A[i, j].Real; yRe[i + rows] = A[i + rows, j].Real;
                            yIm[i] = A[i, j].Imag; yIm[i + rows] = A[i + rows, j].Imag;
                        }
                        var cRe = FastRealWeylHeisenbergTransform.FRWHT(yRe, cacheCols);
                        var cIm = FastRealWeylHeisenbergTransform.FRWHT(yIm, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            tmp[i, j] = new Complex32(cRe[i], cIm[i]);
                    });

                    var B = new Complex32[rows2, cols2];

                    Parallel.For(0, rows2, i =>
                    {
                        var yRe = new float[cols2];
                        var yIm = new float[cols2];
                        for (int j = 0; j < cols; j++)
                        {
                            yRe[j] = tmp[i, j].Real; yRe[j + cols] = tmp[i, j + cols].Real;
                            yIm[j] = tmp[i, j].Imag; yIm[j + cols] = tmp[i, j + cols].Imag;
                        }
                        var cRe = FastRealWeylHeisenbergTransform.FRWHT(yRe, cacheRows);
                        var cIm = FastRealWeylHeisenbergTransform.FRWHT(yIm, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = new Complex32(cRe[j], cIm[j]);
                    });

                    return B;
                }

                if (Direction == Direction.Vertical)
                {
                    var B = new Complex32[rows2, cols2];
                    Parallel.For(0, cols2, j =>
                    {
                        var yRe = new float[rows2];
                        var yIm = new float[rows2];
                        for (int i = 0; i < rows; i++)
                        {
                            yRe[i] = A[i, j].Real; yRe[i + rows] = A[i + rows, j].Real;
                            yIm[i] = A[i, j].Imag; yIm[i + rows] = A[i + rows, j].Imag;
                        }
                        var cRe = FastRealWeylHeisenbergTransform.FRWHT(yRe, cacheCols);
                        var cIm = FastRealWeylHeisenbergTransform.FRWHT(yIm, cacheCols);
                        for (int i = 0; i < rows2; i++)
                            B[i, j] = new Complex32(cRe[i], cIm[i]);
                    });
                    return B;
                }
                else // Horizontal
                {
                    var B = new Complex32[rows2, cols2];
                    Parallel.For(0, rows2, i =>
                    {
                        var yRe = new float[cols2];
                        var yIm = new float[cols2];
                        for (int j = 0; j < cols; j++)
                        {
                            yRe[j] = A[i, j].Real; yRe[j + cols] = A[i, j + cols].Real;
                            yIm[j] = A[i, j].Imag; yIm[j + cols] = A[i, j + cols].Imag;
                        }
                        var cRe = FastRealWeylHeisenbergTransform.FRWHT(yRe, cacheRows);
                        var cIm = FastRealWeylHeisenbergTransform.FRWHT(yIm, cacheRows);
                        for (int j = 0; j < cols2; j++)
                            B[i, j] = new Complex32(cRe[j], cIm[j]);
                    });
                    return B;
                }
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

            if (SpectrumType == SpectrumType.Fourier)
            {
                FastWeylHeisenbergTransform.PolyphaseCache cacheCols = null;
                FastWeylHeisenbergTransform.PolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastWeylHeisenbergTransform.PolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastWeylHeisenbergTransform.PolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
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

                if (Direction == Direction.Vertical)
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
            else // SpectrumType.Hartley
            {
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheCols = null;
                FastRealWeylHeisenbergTransform.RealPolyphaseCache cacheRows = null;

                if (Direction == Direction.Both || Direction == Direction.Vertical)
                    cacheCols = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(rows, this.m, this.window);
                if (Direction == Direction.Both || Direction == Direction.Horizontal)
                    cacheRows = FastRealWeylHeisenbergTransform.RealPolyphaseCache.Build(cols, this.m, this.window);

                if (Direction == Direction.Both)
                {
                    var tmp = new Complex32[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var cRe = new float[rows2];
                        var cIm = new float[rows2];
                        for (int i = 0; i < rows2; i++)
                        {
                            cRe[i] = B[i, j].Real;
                            cIm[i] = B[i, j].Imag;
                        }

                        var yRe = FastRealWeylHeisenbergTransform.IFRWHT(cRe, cacheCols); // [top; bottom]
                        var yIm = FastRealWeylHeisenbergTransform.IFRWHT(cIm, cacheCols);

                        for (int i = 0; i < rows; i++)
                        {
                            tmp[i, j] = new Complex32(yRe[i], yIm[i]);
                            tmp[i + rows, j] = new Complex32(yRe[i + rows], yIm[i + rows]);
                        }
                    });

                    var A = new Complex32[rows2, cols2];

                    Parallel.For(0, rows2, i =>
                    {
                        var cRe = new float[cols2];
                        var cIm = new float[cols2];
                        for (int j = 0; j < cols2; j++)
                        {
                            cRe[j] = tmp[i, j].Real;
                            cIm[j] = tmp[i, j].Imag;
                        }

                        var yRe = FastRealWeylHeisenbergTransform.IFRWHT(cRe, cacheRows);
                        var yIm = FastRealWeylHeisenbergTransform.IFRWHT(cIm, cacheRows);

                        for (int j = 0; j < cols; j++)
                        {
                            A[i, j] = new Complex32(yRe[j], yIm[j]);
                            A[i, j + cols] = new Complex32(yRe[j + cols], yIm[j + cols]);
                        }
                    });

                    return A;
                }

                if (Direction == Direction.Vertical)
                {
                    var A = new Complex32[rows2, cols2];

                    Parallel.For(0, cols2, j =>
                    {
                        var cRe = new float[rows2];
                        var cIm = new float[rows2];
                        for (int i = 0; i < rows2; i++)
                        {
                            cRe[i] = B[i, j].Real;
                            cIm[i] = B[i, j].Imag;
                        }

                        var yRe = FastRealWeylHeisenbergTransform.IFRWHT(cRe, cacheCols);
                        var yIm = FastRealWeylHeisenbergTransform.IFRWHT(cIm, cacheCols);

                        for (int i = 0; i < rows; i++)
                        {
                            A[i, j] = new Complex32(yRe[i], yIm[i]);
                            A[i + rows, j] = new Complex32(yRe[i + rows], yIm[i + rows]);
                        }
                    });

                    return A;
                }
                else // Horizontal
                {
                    var A = new Complex32[rows2, cols2];

                    Parallel.For(0, rows2, i =>
                    {
                        var cRe = new float[cols2];
                        var cIm = new float[cols2];
                        for (int j = 0; j < cols2; j++)
                        {
                            cRe[j] = B[i, j].Real;
                            cIm[j] = B[i, j].Imag;
                        }

                        var yRe = FastRealWeylHeisenbergTransform.IFRWHT(cRe, cacheRows);
                        var yIm = FastRealWeylHeisenbergTransform.IFRWHT(cIm, cacheRows);

                        for (int j = 0; j < cols; j++)
                        {
                            A[i, j] = new Complex32(yRe[j], yIm[j]);
                            A[i, j + cols] = new Complex32(yRe[j + cols], yIm[j + cols]);
                        }
                    });

                    return A;
                }
            }
        }
        #endregion
    }
}