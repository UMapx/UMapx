using System;
using System.Linq;
using System.Reflection;
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
        /// Prints a value to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Value</param>
        public static void Print<T>(this T A)
        {
            Console.WriteLine(A);
            Console.WriteLine();
        }
        /// <summary>
        /// Prints a vector to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Array</param>
        /// <param name="vertical">Vertical or not</param>
        public static void Print<T>(this T[] A, bool vertical = false)
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
        /// Prints a matrix to console.
        /// </summary>
        /// <typeparam name="T">Type</typeparam>
        /// <param name="A">Matrix</param>
        public static void Print<T>(this T[,] A)
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

        #region Info methods
        /// <summary>
        /// Prints a simple reflection-based summary of an object's public properties and methods.
        /// </summary>
        /// <param name="T">
        /// The target instance. (Note: parameter name is uppercase by design here; typically it's named <c>obj</c>.)
        /// </param>
        /// <param name="includeInherited">
        /// If <c>true</c>, include members inherited from base types; otherwise, only members declared on the object's exact type are included.
        /// </param>
        /// <param name="includeStatic">
        /// If <c>true</c>, include static members in addition to instance members; otherwise, include instance members only.
        /// </param>
        /// <param name="withSignatures">
        /// If <c>true</c>, print method signatures (method name + parameter types); otherwise, print method names only and collapse overloads by name.
        /// </param>
        /// <remarks>
        /// This method writes to the console:
        /// 1) The runtime type of the object,
        /// 2) A comma-separated list of public property names,
        /// 3) A comma-separated list of public method names (or signatures).
        /// Special-name methods (property accessors, event add/remove, operators) and constructors are excluded.
        /// </remarks>
        public static void Info(this object T, bool includeInherited = false, bool includeStatic = false, bool withSignatures = false)
        {
            // No-op if the target is null
            if (T is null) return;

            var t = T.GetType();

            // Build BindingFlags based on options:
            // - Always include Public members
            // - Include instance members; optionally include static ones
            // - Optionally restrict to members declared only on this type (exclude inherited)
            var flags = BindingFlags.Public;
            flags |= includeStatic ? BindingFlags.Instance | BindingFlags.Static
                                   : BindingFlags.Instance;
            if (!includeInherited) flags |= BindingFlags.DeclaredOnly;

            // Collect property names
            var props = t.GetProperties(flags).Select(p => p.Name);

            // Collect methods excluding special-name helpers (get_/set_, add_/remove_) and constructors
            var methodsQuery = t.GetMethods(flags)
                .Where(m => !m.IsSpecialName && !m.IsConstructor);

            // Either print full signatures or just method names (distinct to collapse overloads)
            var methods = withSignatures
                ? methodsQuery.Select(m =>
                    $"{m.Name}({string.Join(", ", m.GetParameters().Select(p => p.ParameterType.Name))})")
                : methodsQuery.Select(m => m.Name).Distinct();

            // Output
            Console.WriteLine($"{T.GetType()}");
            Console.WriteLine($"Properties: {string.Join(", ", props)}");
            Console.WriteLine($"Methods: {string.Join(", ", methods)}");
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
        /// Prints elapsed time in miliseconds.
        /// </summary>
        /// <returns>Int</returns>
        public static void Toc()
        {
            Console.WriteLine($"Elapsed time is {Environment.TickCount - tic} milliseconds");
            Console.WriteLine();
        }
        #endregion
    }
}
