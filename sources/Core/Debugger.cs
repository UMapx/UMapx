using System;
using System.Text;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to debug data such as matrices and vectors.
    /// </summary>
    public static class Debugger
    {
        #region Private data
        /// <summary>
        /// Space between the columns.
        /// </summary>
        private const int Gap = 3;
        /// <summary>
        /// Two leading spaces.
        /// </summary>
        private const string RowIndent = "  ";
        /// <summary>
        /// Private tic.
        /// </summary>
        private static int tic;
        #endregion

        #region Print methods
        /// <summary>
        /// Print value to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Value</param>
        public static void Print<T>(T A)
        {
            Console.WriteLine(A);
            Console.WriteLine();
        }
        /// <summary>
        /// Print vector to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Array</param>
        /// <param name="vertical">Vertical or not</param>
        public static void Print<T>(T[] A, bool vertical = false)
        {
            int n = A?.Length ?? 0;

            // Empty case
            if (n == 0)
            {
                Console.WriteLine("Empty vector 0");
                return;
            }
            else
            {
                Console.WriteLine($"Vector {n}");
            }

            // Choose the direction
            if (vertical)
            {
                // Pre-format values and compute single-column width
                var s = new string[n];
                int colWidth = 0;
                for (int i = 0; i < n; i++)
                {
                    string str = A[i]?.ToString() ?? string.Empty;
                    s[i] = str;
                    if (str.Length > colWidth) colWidth = str.Length;
                }

                // Header for a single column (keeps style consistent with matrix print)
                Console.WriteLine("Columns 1 through 1");

                // Print each element on its own line, right-aligned
                for (int i = 0; i < n; i++)
                {
                    var line = new StringBuilder();
                    line.Append(RowIndent);              // uses same indent as your matrix/vector print
                    line.Append(s[i].PadLeft(colWidth)); // align to the right like MATLAB
                    Console.WriteLine(line.ToString());
                }
            }
            else
            {
                // Pre-format values and compute per-column widths
                var s = new string[n];
                var colWidth = new int[n];
                for (int j = 0; j < n; j++)
                {
                    string str = A[j]?.ToString() ?? string.Empty;
                    s[j] = str;
                    colWidth[j] = str.Length;
                }

                int consoleWidth;
                try { consoleWidth = Math.Max(20, Console.WindowWidth); }
                catch { consoleWidth = 80; } // fallback for hosts without a real console
                int startCol = 0;

                while (startCol < n)
                {
                    // Choose how many columns fit into the current console width
                    int widthUsed = RowIndent.Length;
                    int endCol = startCol;

                    while (endCol < n)
                    {
                        int add = (endCol > startCol ? Gap : 0) + colWidth[endCol];
                        if (widthUsed + add > consoleWidth)
                        {
                            if (endCol == startCol) endCol++; // ensure at least one column
                            break;
                        }
                        widthUsed += add;
                        endCol++;
                    }
                    if (endCol <= startCol) endCol = startCol + 1;

                    if (startCol > 0) Console.WriteLine();
                    Console.WriteLine($"Columns {startCol + 1} through {endCol}");

                    // Print rows for the chosen block
                    var line = new StringBuilder();
                    line.Append(RowIndent);

                    for (int j = startCol; j < endCol; j++)
                    {
                        if (j > startCol) line.Append(' ', Gap);
                        line.Append(s[j].PadLeft(colWidth[j]));
                    }
                    Console.WriteLine(line.ToString());

                    startCol = endCol;
                }
            }
            Console.WriteLine();
        }
        /// <summary>
        /// Print matrix to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Matrix</param>
        public static void Print<T>(T[,] A)
        {
            int m = A?.GetLength(0) ?? 0, n = A?.GetLength(1) ?? 0;

            // Empty cases
            if (m == 0 || n == 0)
            {
                Console.WriteLine($"Empty matrix {m} x {n}");
                return;
            }
            else
            {
                Console.WriteLine($"Matrix {m} x {n}");
            }

            var s = new string[m, n];
            var colWidth = new int[n];

            // Pre-format numbers and compute per-column widths
            for (int j = 0; j < n; j++)
            {
                int w = 0;
                for (int i = 0; i < m; i++)
                {
                    string str = A[i, j].ToString();
                    s[i, j] = str;
                    if (str.Length > w) w = str.Length;
                }
                colWidth[j] = w;
            }

            int consoleWidth;
            try { consoleWidth = Math.Max(20, Console.WindowWidth); }
            catch { consoleWidth = 80; } // fallback for hosts without a real console
            int startCol = 0;

            while (startCol < n)
            {
                // Choose how many columns fit into the current console width
                int widthUsed = RowIndent.Length;
                int endCol = startCol;

                while (endCol < n)
                {
                    int add = (endCol > startCol ? Gap : 0) + colWidth[endCol];
                    if (widthUsed + add > consoleWidth)
                    {
                        if (endCol == startCol) endCol++; // ensure at least one column
                        break;
                    }
                    widthUsed += add;
                    endCol++;
                }

                if (endCol <= startCol)
                    endCol = startCol + 1;

                if (startCol > 0)
                {
                    Console.WriteLine();
                }
                Console.WriteLine($"Columns {startCol + 1} through {endCol}");

                // Print rows for the chosen block
                for (int i = 0; i < m; i++)
                {
                    var line = new StringBuilder();
                    line.Append(RowIndent);
                    for (int j = startCol; j < endCol; j++)
                    {
                        if (j > startCol) line.Append(' ', Gap);
                        line.Append(s[i, j].PadLeft(colWidth[j]));
                    }
                    Console.WriteLine(line.ToString());
                }

                // Next block
                startCol = endCol;
            }
            Console.WriteLine();
        }
        #endregion

        #region Tic toc methods
        /// <summary>
        /// Starts elapsing time.
        /// </summary>
        public static void Tic()
        {
            tic = Environment.TickCount;
        }
        /// <summary>
        /// Returns elapsed time in miliseconds.
        /// </summary>
        /// <returns>Int</returns>
        public static int Toc()
        {
            return Environment.TickCount - tic;
        }
        #endregion
    }
}
