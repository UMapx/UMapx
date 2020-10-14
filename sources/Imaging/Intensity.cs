using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with brightness represented as a value belonging to the interval [0, 1].
    /// </summary>
    public static class Intensity
    {
        #region Nonlinear methods components
        /// <summary>
        /// Implements the Single Scale Retinex algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double SingleScaleRetinex(double x, double xlow, double nbase, double a, double b)
        {
            // Singe scale retinex modified algorithm
            // by Valery Asiryan
            // 
            return Math.Exp(a * Math.Log(x / xlow, nbase) + b) - 0.5f;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] SingleScaleRetinex(double nbase, double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.SingleScaleRetinex(w, v, nbase, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the local contrast inversion algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double LocalContrastInversion(double x, double xlow, double a, double b)
        {
            return a * x / (xlow + b);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LocalContrastInversion(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LocalContrastInversion(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the local contrast enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="a">Factor [-1, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double LocalContrastEnhancement(double x, double xlow, double a)
        {
            return x + a * (x - xlow);

            //return (1.0 + a) * (x - xlow) + xlow;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LocalContrastEnhancement(double a, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LocalContrastEnhancement(w, v, a);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the homomorphic enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double HomomorphicEnhancement(double x, double mu, double a, double b)
        {
            return Math.Exp(Math.Log(x) - a * Math.Log(mu + b));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] HomomorphicEnhancement(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.HomomorphicEnhancement(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the ξ-contrast enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double KsiContrastEnchancement(double x, double mu, double a, double b)
        {
            // x ∈ [0, 1], μ ∈ [0, 1] - mean of x.
            // σ - variance, ξ - coefficient.
            //
            // σ = x - μ
            // ξ = σ / μ
            // result value:
            // x' = x + α * ξ + β, where α, β ∈ [-1, 1].

            double sigma = x - mu;
            double ksi = sigma / mu;
            return x + a * ksi + b;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] KsiContrastEnchancement(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.KsiContrastEnchancement(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the SAUCE algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="d">Degree of difference [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double SAUCE(double x, double mu, double d)
        {
            // Ravimal Bandara algorithm
            // implementation:

            double a = (mu - d / 2.0);
            double b = (mu + d / 2.0);

            if (x < a)
            {
                return 0;
            }
            else if (x >= b)
            {
                return 1;
            }
            return (x - a) / d;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] SAUCE(double a, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.SAUCE(w, v, a);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the Bradley threshold correction.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="difference">Difference [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bradley(double x, double xlow, double difference = 0.15)
        {
            // Bradley local threshold void.
            // Derek Bradley, Gerhard Roth (2005). Adaptive Thresholding Using the Integral Image.
            // Retrieved from http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf

            double z = 1.0 - difference;
            return (x < xlow * z) ? 0 : 1;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="difference">Difference [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] Bradley(double difference, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.Bradley(w, v, difference);
                }
            }

            return table;
        }
        #endregion

        #region LogStretrch Histogram method components
        /// <summary>
        /// Logarithm of 0.5.
        /// </summary>
        public const double log05 = -0.693147180559945;
        /// <summary>
        /// Logarithmic epsilon.
        /// </summary>
        public const double logEpsilon = 1e-9;
        /// <summary>
        /// Implements the logarithmic stretch algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogStretch(double x, double mu, double s, double l)
        {
            return Intensity.LogPow(x, Maths.Range(Intensity.log05 / Math.Log(mu), s, l));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LogStretch(double s, double l, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LogStretch(w, v, s, l);
                }
            }

            return table;
        }
        /// <summary>
        /// Returns the number raised to the logarithmic power.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="power">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogPow(double a, double power)
        {
            return Math.Exp(Math.Log(a) * power);
        }
        #endregion

        #region Linear methods components
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="g">Gamma</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Gamma(double g, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Gamma(x / (double)length, g);
            }
            return table;
        }
        /// <summary>
        /// Implements the gamma correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="g">Gamma</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gamma(double x, double g)
        {
            return Math.Pow(x, g);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Shift(double b, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Shift(x / (double)length, b);
            }
            return table;
        }
        /// <summary>
        /// Implements the shift correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <returns>Double precision floating point number</returns>
        public static double Shift(double x, double b)
        {
            double v = log05 / Math.Log(0.5 - b);
            return LogPow(x, v);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Bin(double threshold, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Bin(x / (double)length, threshold);
            }
            return table;
        }
        /// <summary>
        /// Implements the threshold correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bin(double x, double threshold)
        {
            return (x > threshold) ? 1 : 0;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="average">Average</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Exposure(double average, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Exposure(x / (double)length, average);
            }
            return table;
        }
        /// <summary>
        /// Implements the exposure correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="average">Average</param>
        /// <returns>Double precision floating point number</returns>
        public static double Exposure(double x, double average)
        {
            double T = 255.0 / average;
            return 1 - Math.Exp(-T * x);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Sin(double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Sin(x / (double)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the sine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sin(double x, double delta)
        {
            return 0.5 * Math.Sin((3.14 * x) - (3.14 / 2)) + 0.5 + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Cos(double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Cos(x / (double)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the cosine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cos(double x, double delta)
        {
            return 0.5 * Math.Cos((3.14 * x) - 3.14) + 0.5 + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Log(double a, double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Log(x / (double)length, a, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the logarithmic correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log(double x, double a, double delta)
        {
            return Math.Log((1.0 + (x + delta) / 0.5), a);
        }
        /// <summary>
        /// Returns the correction mask for formula: Y = (X + V).
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Add(double value, int length)
        {
            double[] table = new double[length];
            for (int x = 0; x < length; x++)
            {
                table[x] = x / (double)length + value;
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Contrast(double value, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Contrast(x / (double)length, value);
            }
            return table;
        }
        /// <summary>
        /// Implements the contrast correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="value">Contrast</param>
        /// <returns>Double precision floating point number</returns>
        public static double Contrast(double x, double value)
        {
            value = (1 + value);
            double xc = x;
            xc -= 0.5f;
            xc *= value;
            xc += 0.5f;
            return xc;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="power">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] LogContrast(double power, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LogContrast(x / (double)length, power);
            }
            return table;
        }
        /// <summary>
        /// Implements the log-contrast correction.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="power">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogContrast(double x, double power)
        {
            if (x <= 0.5)
            {
                return Intensity.LogPow(x * 2, power) * 0.5;
            }
            return 1.0 - Intensity.LogPow((1 - x) * 2, power) * 0.5;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Invert(int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Invert(x / (double)length);
            }
            return table;
        }
        /// <summary>
        /// Negates the value.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Invert(double x)
        {
            return 1.0 - x;
        }
        /// <summary>
        /// Equalizes a value relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <returns>Double precision floating point number</returns>
        public static double Equalize(double x, double min, double max)
        {
            double a = max - min;
            double b = x - min;
            double c = (a != 0) ? b / a : x;
            return Maths.Double(c);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Equalize(double min, double max, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Equalize(x / (double)length, min, max);
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="range">Pair of numbers Max и Min</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Linear(RangeDouble range, double delta, int length)
        {
            return Intensity.Linear(range.Max, range.Min, delta, length);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="xmax">Maximum value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Linear(double xmax, double xmin, double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Linear(x / (double)length, xmax, xmin, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the linear correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="xmax">Maximum value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Linear(double x, double xmax, double xmin, double delta)
        {
            return (x - xmin) / (xmax - xmin) + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="xmin">Minimum value of the input range</param>
        /// <param name="xmax">Maximum value of the input range</param>
        /// <param name="ymin">Minimum value of the output range</param>
        /// <param name="ymax">Maximum value of the output range</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Levels(double xmin, double xmax, double ymin, double ymax, int length)
        {
            double[] table = new double[length];
            double k = 0, b = 0;
            double v;

            if (xmax != xmin)
            {
                k = (ymax - ymin) / (xmax - xmin);
                b = (ymin) - k * xmin;
            }

            for (int i = 0; i < length; i++)
            {
                v = i / (double)length;

                if (v >= xmax)
                { v = ymax; }
                else if (v <= xmin)
                { v = ymin; }
                else
                { v = k * v + b; }

                table[i] = v;
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="input">Input values</param>
        /// <param name="output">Output values</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Levels(RangeDouble input, RangeDouble output, int length)
        {
            return Intensity.Levels(input.Min, input.Max, output.Min, output.Max, length);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Quantize(int levels, int length)
        {
            if (levels > length)
                throw new Exception("Number of levels cannot be greater than length");

            int interval = length / levels + 1;
            double[] table = new double[length];
            double min = double.MaxValue;
            double max = double.MinValue;
            double v;
            int q, i;

            // calculating
            for (i = 0; i < length; i++)
            {
                q = (i / interval) * interval;
                v = q / (length - 1.0);

                if (v < min)
                    min = v;

                if (v > max)
                    max = v;

                table[i] = v;
            }

            // normalizing
            v = (max - min);

            for (i = 0; i < length; i++)
            {
                table[i] = (table[i] - min) / v;
            }

            return table;
        }
        #endregion
    }
}
