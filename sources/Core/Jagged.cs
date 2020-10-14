using System;

namespace UMapx.Core
{
    /// <summary>
    /// Uses to work with gear arrays.
    /// </summary>
    public static class Jagged
    {
        #region Conversions
        /// <summary>
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Jagged array</returns>
        public static double[][] ToJagged(this double[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            double[][] jagged = new double[ml][];
            double[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new double[mr];
                for (j = 0; j < mr; j++)
                {
                    data[j] = m[i, j];
                }
                jagged[i] = data;
            }
            return jagged;
        }
        /// <summary>
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
        public static double[,] FromJagged(this double[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            double[,] m = new double[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    m[i, j] = jagged[i][j];
                }
            }
            return m;
        }
        /// <summary>
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Jagged array</returns>
        public static Complex[][] ToJagged(this Complex[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex[][] jagged = new Complex[ml][];
            Complex[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new Complex[mr];
                for (j = 0; j < mr; j++)
                {
                    data[j] = m[i, j];
                }
                jagged[i] = data;
            }
            return jagged;
        }
        /// <summary>
        /// Returns matrix.
        /// </summary>
        /// <param name="jagged">Jagged array</param>
        /// <returns>Matrix</returns>
        public static Complex[,] FromJagged(this Complex[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            Complex[,] m = new Complex[ml, mr];
            int i, j;

            for (i = 0; i < ml; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    m[i, j] = jagged[i][j];
                }
            }
            return m;
        }
        #endregion

        #region Jagged array
        /// <summary>
        /// 
        /// </summary>
        private static Random rnd = new Random();
        /// <summary>
        /// Implements the construction of a vector of random numbers, the values of which are distributed according to a uniform distribution.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[][] Rand(int m, int l)
        {
            double[][] H = new double[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new double[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = rnd.NextDouble();
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a vector of random numbers, the values of which are distributed according to a uniform distribution.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Randc(int m, int l)
        {
            Complex[][] H = new Complex[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new Complex[l];
                for (j = 0; j < l; j++)
                {
                    H[i][j] = new Complex(rnd.NextDouble(), rnd.NextDouble());
                }
            }

            return H;
        }

        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[][] Randi(int m, int l)
        {
            return Randi(m, l, 1, l + 1);
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Matrix</returns>
        public static double[][] Randi(int m, int l, int a, int b)
        {
            double[][] H = new double[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new double[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = rnd.Next(a, b);
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Randic(int m, int l)
        {
            return Randic(m, l, 1, l + 1);
        }
        /// <summary>
        /// Implements the construction of a vector of integer random numbers.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Randic(int m, int l, int a, int b)
        {
            Complex[][] H = new Complex[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new Complex[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = new Complex(rnd.Next(a, b), rnd.Next(a, b));
                }
            }

            return H;
        }

        /// <summary>
        /// Implements the construction of a zero matrix.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[][] Zero(int m, int l)
        {
            double[][] H = new double[m][];
            int i;

            for (i = 0; i < m; i++)
            {
                H[i] = new double[l];
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of ones.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[][] One(int m, int l)
        {
            double[][] H = new double[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new double[l];
                for (j = 0; j < l; j++)
                {
                    H[i][j] = 1.0;
                }
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a eye matrix.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static double[][] Eye(int m, int l)
        {
            double[][] H = new double[m][];
            int i;

            for (i = 0; i < m; i++)
            {
                H[i] = new double[l];
                H[i][i] = 1.0;
            }

            return H;
        }
        #endregion

        #region Parse methods
        /// <summary>
        /// Parses the original string into a matrix of double numbers.
        /// <remarks>
        /// Example: "[1, 2, 3; 4, 5, 6; 7, 8, 9]";
        /// </remarks>
        /// </summary>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static double[][] Parse(this double[][] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums;
            int r = rows.Length, n;
            double[][] H = new double[r][];
            int i, j;

            // collecting rows:
            for (i = 0; i < r; i++)
            {
                nums = rows[i].Split('|');
                n = nums.Length;
                H[i] = new double[n];

                for (j = 0; j < n; j++)
                {
                    H[i][j] = double.Parse(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of double numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, ref double[][] result)
        {
            double[][] zero = null;
            try
            {
                result = Jagged.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }
        /// <summary>
        /// Parses the original string into a matrix of complex numbers.
        /// </summary>
        /// <remarks>
        /// Example: "[1 + 2i, 2 + 4i; 3 + 6i, 4 + 8i]";
        /// </remarks>
        /// <param name="a">Matrix</param>
        /// <param name="s">Input string</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Parse(this Complex[][] a, string s)
        {
            string[] rows = StringOptions.Matpar(s);
            string[] nums;
            int r = rows.Length, n;
            Complex[][] H = new Complex[r][];
            int i, j;

            // collecting rows:
            for (i = 0; i < r; i++)
            {
                nums = rows[i].Split('|');
                n = nums.Length;
                H[i] = new Complex[n];

                for (j = 0; j < n; j++)
                {
                    H[i][j] = StringOptions.Compar(nums[j]);
                }
            }
            return H;
        }
        /// <summary>
        /// Tries to parse the original row into a matrix of complex numbers.
        /// </summary>
        /// <param name="s">Input string</param>
        /// <param name="result">Matrix</param>
        /// <returns>Boolean</returns>
        public static bool TryParse(string s, ref Complex[][] result)
        {
            Complex[][] zero = null;
            try
            {
                result = Jagged.Parse(zero, s);
                return true;
            }
            catch (FormatException)
            {
                result = zero;
                return false;
            }
        }
        #endregion

        #region Matrix conversions
        /// <summary>
        /// Negates all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Negate(this double[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            double[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = -v[j];
                }
            }

            return H;
        }
        /// <summary>
        /// Negates all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[][] Negate(this Complex[][] m)
        {
            int r0 = m.GetLength(0), r1;
            Complex[][] H = new Complex[r0][];
            Complex[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = -v[j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a complex matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static Complex[][] ToComplex(this double[][] m)
        {
            int r0 = m.GetLength(0), r1;
            Complex[][] H = new Complex[r0][];
            double[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = v[j];
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a matrix whose values belong to the interval [0, 255].
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] ToByte(this double[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            double[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = Maths.Byte(v[j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Returns a matrix whose values belong to the interval [0, 1].
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] ToDouble(this double[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            double[] v; double c;
            double min = double.MaxValue, max = double.MinValue;
            int i, j;

            // find min/max
            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    c = v[j];
                    if (c > max) max = c;
                    if (c < min) min = c;
                }
            }

            // scaling to [0, 1] range:
            double range = max - min;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = (v[j] - min) / range;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Abs(this double[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            double[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = Math.Abs(v[j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Takes a module for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Abs(this Complex[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            Complex[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = v[j].Abs;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes an angle for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Angle(this Complex[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            Complex[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = v[j].Angle;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes the real part for all elements of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Real(this Complex[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            Complex[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = v[j].Real;
                }
            }

            return H;
        }
        /// <summary>
        /// Takes the imaginary part for all elements of the matrix.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static double[][] Imag(this Complex[][] m)
        {
            int r0 = m.GetLength(0), r1;
            double[][] H = new double[r0][];
            Complex[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = v[j].Imag;
                }
            }

            return H;
        }
        #endregion
    }
}
