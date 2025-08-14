using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Used to work with brightness represented as a value belonging to the interval [0, 1].
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
        /// <returns>Value</returns>
        public static float SingleScaleRetinex(float x, float xlow, float nbase, float a, float b)
        {
            // Singe scale retinex modified algorithm
            // by Valery Asiryan
            // 
            return a * (float)Math.Log(x / xlow, nbase) + b;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static float[,] SingleScaleRetinex(float nbase, float a, float b, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float LocalContrastInversion(float x, float xlow, float a, float b)
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
        public static float[,] LocalContrastInversion(float a, float b, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float LocalContrastEnhancement(float x, float xlow, float a)
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
        public static float[,] LocalContrastEnhancement(float a, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float HomomorphicEnhancement(float x, float mu, float a, float b)
        {
            return (float)Math.Exp(Math.Log(x) - a * Math.Log(mu + b));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static float[,] HomomorphicEnhancement(float a, float b, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float KsiContrastEnchancement(float x, float mu, float a, float b)
        {
            // x ∈ [0, 1], μ ∈ [0, 1] - mean of x.
            // σ - variance, ξ - coefficient.
            //
            // σ = x - μ
            // ξ = σ / μ
            // result value:
            // x' = x + α * ξ + β, where α, β ∈ [-1, 1].

            float sigma = x - mu;
            float ksi = sigma / mu;
            return x + a * ksi + b;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static float[,] KsiContrastEnchancement(float a, float b, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float SAUCE(float x, float mu, float d)
        {
            // Ravimal Bandara algorithm
            // implementation:

            float a = (mu - d / 2.0f);
            float b = (mu + d / 2.0f);

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
        public static float[,] SAUCE(float a, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float Bradley(float x, float xlow, float difference = 0.15f)
        {
            // Bradley local threshold void.
            // Derek Bradley, Gerhard Roth (2005). Adaptive Thresholding Using the Integral Image.
            // Retrieved from http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf

            float z = 1.0f - difference;
            return (x < xlow * z) ? 0 : 1;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="difference">Difference [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static float[,] Bradley(float difference, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        public const float log05 = -0.693147180559945f;
        /// <summary>
        /// Logarithmic epsilon.
        /// </summary>
        public const float logEpsilon = 1e-9f;
        /// <summary>
        /// Implements the logarithmic stretch algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <returns>Value</returns>
        public static float LogStretch(float x, float mu, float s, float l)
        {
            return Intensity.LogPow(x, Maths.Range(Intensity.log05 / (float)Math.Log(mu), s, l));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static float[,] LogStretch(float s, float l, int length)
        {
            float[,] table = new float[length, length];
            float w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (float)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (float)length;

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
        /// <returns>Value</returns>
        public static float LogPow(float a, float power)
        {
            return (float)Math.Exp(Math.Log(a) * power);
        }
        #endregion

        #region Linear methods components
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="g">Gamma</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Gamma(float g, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Gamma(x / (float)length, g);
            }
            return table;
        }
        /// <summary>
        /// Implements the gamma correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="g">Gamma</param>
        /// <returns>Value</returns>
        public static float Gamma(float x, float g)
        {
            return (float)Math.Pow(x, g);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Shift(float b, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Shift(x / (float)length, b);
            }
            return table;
        }
        /// <summary>
        /// Implements the shift correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <returns>Value</returns>
        public static float Shift(float x, float b)
        {
            float v = log05 / (float)Math.Log(0.5 - b);
            return LogPow(x, v);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Bin(float threshold, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Bin(x / (float)length, threshold);
            }
            return table;
        }
        /// <summary>
        /// Implements the threshold correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <returns>Value</returns>
        public static float Bin(float x, float threshold)
        {
            return (x > threshold) ? 1 : 0;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="average">Average</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Exposure(float average, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Exposure(x / (float)length, average);
            }
            return table;
        }
        /// <summary>
        /// Implements the exposure correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="average">Average</param>
        /// <returns>Value</returns>
        public static float Exposure(float x, float average)
        {
            float T = 255.0f / average;
            return 1 - (float)Math.Exp(-T * x);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Sin(float delta, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Sin(x / (float)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the sine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Value</returns>
        public static float Sin(float x, float delta)
        {
            return 0.5f * (float)Math.Sin((3.14 * x) - (3.14 / 2)) + 0.5f + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Cos(float delta, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Cos(x / (float)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the cosine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Value</returns>
        public static float Cos(float x, float delta)
        {
            return 0.5f * (float)Math.Cos((3.14 * x) - 3.14) + 0.5f + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Log(float a, float delta, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Log(x / (float)length, a, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the logarithmic correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <returns>Value</returns>
        public static float Log(float x, float a, float delta)
        {
            return (float)Math.Log((1.0 + (x + delta) / 0.5), a);
        }
        /// <summary>
        /// Returns the correction mask for formula: Y = (X + V).
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Add(float value, int length)
        {
            float[] table = new float[length];
            for (int x = 0; x < length; x++)
            {
                table[x] = x / (float)length + value;
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Contrast(float value, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Contrast(x / (float)length, value);
            }
            return table;
        }
        /// <summary>
        /// Implements the contrast correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="value">Contrast</param>
        /// <returns>Value</returns>
        public static float Contrast(float x, float value)
        {
            value = (1 + value);
            float xc = x;
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
        public static float[] LogContrast(float power, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LogContrast(x / (float)length, power);
            }
            return table;
        }
        /// <summary>
        /// Implements the log-contrast correction.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="power">Power</param>
        /// <returns>Value</returns>
        public static float LogContrast(float x, float power)
        {
            if (x <= 0.5)
            {
                return (float)Intensity.LogPow(x * 2, power) * 0.5f;
            }
            return 1.0f - (float)Intensity.LogPow((1 - x) * 2, power) * 0.5f;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Invert(int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Invert(x / (float)length);
            }
            return table;
        }
        /// <summary>
        /// Negates the value.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Value</returns>
        public static float Invert(float x)
        {
            return 1.0f - x;
        }
        /// <summary>
        /// Equalizes a value relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <returns>Value</returns>
        public static float Equalize(float x, float min, float max)
        {
            float a = max - min;
            float b = x - min;
            float c = (a != 0) ? b / a : x;
            return Maths.Float(c);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Equalize(float min, float max, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Equalize(x / (float)length, min, max);
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
        public static float[] Linear(RangeFloat range, float delta, int length)
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
        public static float[] Linear(float xmax, float xmin, float delta, int length)
        {
            float[] table = new float[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Linear(x / (float)length, xmax, xmin, delta);
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
        /// <returns>Value</returns>
        public static float Linear(float x, float xmax, float xmin, float delta)
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
        public static float[] Levels(float xmin, float xmax, float ymin, float ymax, int length)
        {
            float[] table = new float[length];
            float k = 0, b = 0;
            float v;

            if (xmax != xmin)
            {
                k = (ymax - ymin) / (xmax - xmin);
                b = (ymin) - k * xmin;
            }

            for (int i = 0; i < length; i++)
            {
                v = i / (float)length;

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
        public static float[] Levels(RangeFloat input, RangeFloat output, int length)
        {
            return Intensity.Levels(input.Min, input.Max, output.Min, output.Max, length);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static float[] Quantize(int levels, int length)
        {
            if (levels > length)
                throw new Exception("Number of levels cannot be greater than length");

            int interval = length / levels + 1;
            float[] table = new float[length];
            float min = float.MaxValue;
            float max = float.MinValue;
            float v;
            int q, i;

            // calculating
            for (i = 0; i < length; i++)
            {
                q = (i / interval) * interval;
                v = q / (length - 1.0f);

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
