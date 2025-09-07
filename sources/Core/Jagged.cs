using System;

namespace UMapx.Core
{
    /// <summary>
    /// Used to work with jagged arrays.
    /// </summary>
    public static partial class Jagged
    {
        #region Conversions
        /// <summary>
        /// Returns jagged array.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Jagged array</returns>
        public static float[][] ToJagged(this float[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            float[][] jagged = new float[ml][];
            float[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new float[mr];
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
        public static float[,] FromJagged(this float[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            float[,] m = new float[ml, mr];
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
        public static Complex32[][] ToJagged(this Complex32[,] m)
        {
            int ml = m.GetLength(0), mr = m.GetLength(1);
            Complex32[][] jagged = new Complex32[ml][];
            Complex32[] data;
            int i, j;

            for (i = 0; i < ml; i++)
            {
                data = new Complex32[mr];
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
        public static Complex32[,] FromJagged(this Complex32[][] jagged)
        {
            int ml = jagged.GetLength(0), mr = jagged[0].GetLength(0);
            Complex32[,] m = new Complex32[ml, mr];
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
        /// Random generator.
        /// </summary>
        private static readonly Random rnd = new Random();
        /// <summary>
        /// Constructs a matrix of random numbers with values uniformly distributed in the interval [0, 1).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[][] Rand(int m, int l)
        {
            float[][] H = new float[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new float[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = (float)rnd.NextDouble();
                }
            }

            return H;
        }
        /// <summary>
        /// Constructs a complex matrix of random numbers whose real and imaginary parts are uniformly distributed in the interval [0, 1).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex32[][] Randc(int m, int l)
        {
            Complex32[][] H = new Complex32[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new Complex32[l];
                for (j = 0; j < l; j++)
                {
                    H[i][j] = new Complex32((float)rnd.NextDouble(), (float)rnd.NextDouble());
                }
            }

            return H;
        }

        /// <summary>
        /// Constructs a matrix of integer random numbers uniformly distributed between 1 and l (inclusive).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[][] Randi(int m, int l)
        {
            return Randi(m, l, 1, l + 1);
        }
        /// <summary>
        /// Constructs a matrix of integer random numbers uniformly distributed in the interval [<paramref name="a"/>, <paramref name="b"/>).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound (exclusive)</param>
        /// <returns>Matrix</returns>
        public static float[][] Randi(int m, int l, int a, int b)
        {
            float[][] H = new float[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new float[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = rnd.Next(a, b);
                }
            }

            return H;
        }
        /// <summary>
        /// Constructs a complex matrix of integer random numbers whose real and imaginary parts are uniformly distributed between 1 and l (inclusive).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Complex32[][] Randic(int m, int l)
        {
            return Randic(m, l, 1, l + 1);
        }
        /// <summary>
        /// Constructs a complex matrix of integer random numbers whose real and imaginary parts are uniformly distributed in the interval [<paramref name="a"/>, <paramref name="b"/>).
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="a">Lower bound</param>
        /// <param name="b">Upper bound (exclusive)</param>
        /// <returns>Matrix</returns>
        public static Complex32[][] Randic(int m, int l, int a, int b)
        {
            Complex32[][] H = new Complex32[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new Complex32[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = new Complex32(rnd.Next(a, b), rnd.Next(a, b));
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
        public static float[][] Zero(int m, int l)
        {
            float[][] H = new float[m][];
            int i;

            for (i = 0; i < m; i++)
            {
                H[i] = new float[l];
            }

            return H;
        }
        /// <summary>
        /// Implements the construction of a matrix of ones.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static float[][] One(int m, int l)
        {
            float[][] H = new float[m][];
            int i, j;

            for (i = 0; i < m; i++)
            {
                H[i] = new float[l];

                for (j = 0; j < l; j++)
                {
                    H[i][j] = 1.0f;
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
        public static float[][] Eye(int m, int l)
        {
            float[][] H = new float[m][];

            for (int i = 0; i < m; i++)
            {
                H[i] = new float[l];

                for (int j = 0; j < l; j++)
                {
                    H[i][j] = i == j ? 1 : 0;
                }
            }

            return H;
        }
        #endregion

        #region Matrix conversions
        /// <summary>
        /// Negates all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[][] Negate(this float[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            float[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
        public static Complex32[][] Negate(this Complex32[][] m)
        {
            int r0 = m.GetLength(0), r1;
            Complex32[][] H = new Complex32[r0][];
            Complex32[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new Complex32[r1];

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
        public static Complex32[][] ToComplex(this float[][] m)
        {
            int r0 = m.GetLength(0), r1;
            Complex32[][] H = new Complex32[r0][];
            float[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new Complex32[r1];

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
        public static float[][] ToByte(this float[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            float[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
        public static float[][] ToFloat(this float[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            float[] v; float c;
            float min = float.MaxValue, max = float.MinValue;
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
            float range = max - min;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = (v[j] - min) / range;
                }
            }

            return H;
        }
        /// <summary>
        /// Calculates the modulus for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[][] Abs(this float[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            float[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

                for (j = 0; j < r1; j++)
                {
                    H[i][j] = Math.Abs(v[j]);
                }
            }

            return H;
        }
        /// <summary>
        /// Calculates the modulus for all matrix elements.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Matrix</returns>
        public static float[][] Abs(this Complex32[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            Complex32[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
        public static float[][] Angle(this Complex32[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            Complex32[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
        public static float[][] Real(this Complex32[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            Complex32[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
        public static float[][] Imag(this Complex32[][] m)
        {
            int r0 = m.GetLength(0), r1;
            float[][] H = new float[r0][];
            Complex32[] v;
            int i, j;

            for (i = 0; i < r0; i++)
            {
                v = m[i];
                r1 = v.GetLength(0);
                H[i] = new float[r1];

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
