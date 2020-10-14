using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines a local averaging filter.
    /// </summary>
    internal static class BoxFilterOptions
    {
        #region Public static voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="r">Radius</param>
        public static void boxf(double[,] data, Direction direction, int r)
        {
            int N = data.GetLength(0);
            int M = data.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    BoxFilterOptions.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    BoxFilterOptions.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    double[] col = new double[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    BoxFilterOptions.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    double[] row = new double[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    BoxFilterOptions.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        /// <param name="direction">Processing direction</param>
        /// <param name="r">Radius</param>
        public static void boxf(Complex[,] data, Direction direction, int r)
        {
            int N = data.GetLength(0);
            int M = data.GetLength(1);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    BoxFilterOptions.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                }
                );

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    BoxFilterOptions.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = data[i, j];
                    }

                    BoxFilterOptions.boxf(col, N, r);

                    for (i = 0; i < N; i++)
                    {
                        data[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = data[i, j];
                    }

                    BoxFilterOptions.boxf(row, M, r);

                    for (j = 0; j < M; j++)
                    {
                        data[i, j] = row[j];
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="l">Length</param>
        /// <param name="r">Radius</param>
        public static void boxf(double[] input, int l, int r)
        {
            if (l == 1) return;
            int h = r >= l ? l - 1 : r;

            int v = h >> 1;
            int dl = l - v;
            double s = 0;
            int i;

            for (i = 0; i < h; i++)
            {
                s += input[i];
            }
            for (i = 0; i < v; i++)
            {
                input[i] = s / h;
            }
            for (i = v; i < dl; i++)
            {
                s = s - input[i - v] + input[i + v];
                input[i] = s / h;
            }
            for (i = dl; i < l; i++)
            {
                s = s - input[i - v] + input[i];
                input[i] = s / h;
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="l">Length</param>
        /// <param name="r">Radius</param>
        public static void boxf(Complex[] input, int l, int r)
        {
            if (l == 1) return;
            int h = r >= l ? l - 1 : r;

            int v = h >> 1;
            int dl = l - v;
            Complex s = 0;
            int i;

            for (i = 0; i < h; i++)
            {
                s += input[i];
            }
            for (i = 0; i < v; i++)
            {
                input[i] = s / h;
            }
            for (i = v; i < dl; i++)
            {
                s = s - input[i - v] + input[i + v];
                input[i] = s / h;
            }
            for (i = dl; i < l; i++)
            {
                s = s - input[i - v] + input[i];
                input[i] = s / h;
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix array</param>
        /// <returns>Matrix</returns>
        public static double[,] boxf(double[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            double[,] sum = new double[r, c], cur;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                cur = data[i];

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize all signals:
                        sum[j, k] += cur[j, k] / length;
                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix array</param>
        /// <returns>Matrix</returns>
        public static Complex[,] boxf(Complex[][,] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0), c = data[0].GetLength(1);
            Complex[,] sum = new Complex[r, c], cur;
            int i, j, k;

            // process
            for (i = 0; i < length; i++)
            {
                cur = data[i];

                for (j = 0; j < r; j++)
                {
                    for (k = 0; k < c; k++)
                    {
                        // summarize all signals:
                        sum[j, k] += cur[j, k] / length;
                    }
                }
            }

            return sum;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Jagged array</param>
        /// <returns>Array</returns>
        public static double[] boxf(double[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            double[] sum = new double[r], cur;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                cur = data[i];

                for (j = 0; j < r; j++)
                {
                    // summarize all signals:
                    sum[j] += cur[j] / length;
                }
            }

            return sum;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Jagged array</param>
        /// <returns>Array</returns>
        public static Complex[] boxf(Complex[][] data)
        {
            // exception
            int length = data.Length;
            if (length == 0) return null;

            // data
            int r = data[0].GetLength(0);
            Complex[] sum = new Complex[r], cur;
            int i, j;

            // process
            for (i = 0; i < length; i++)
            {
                cur = data[i];

                for (j = 0; j < r; j++)
                {
                    // summarize all signals:
                    sum[j] += cur[j] / length;
                }
            }

            return sum;
        }
        #endregion
    }
}
