using System;
using SkiaDrawing;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses to work with the statistical characteristics of the image.
    /// </summary>
    public static class Statistics
    {
        #region Histogram
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Array</returns>
        public unsafe static int[] Histogram(this BitmapData bmData)
        {
            int[] rgb = new int[256];
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            byte brightness;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    brightness = (byte)RGB.Average(p[2], p[1], p[0]);
                    rgb[brightness]++;
                }
            }
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="channel">Channel of RGBA model</param>
        /// <returns>Array</returns>
        public unsafe static int[] Histogram(this BitmapData bmData, RGBA channel)
        {
            int[] rgb = new int[256];
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            if (channel == RGBA.Red)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[2])]++;
                    }
                }
            }
            else if (channel == RGBA.Green)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[1])]++;
                    }
                }
            }
            else if (channel == RGBA.Blue)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[0])]++;
                    }
                }
            }
            else if (channel == RGBA.Alpha)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[3])]++;
                    }
                }
            }
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Array</returns>
        public static int[] Histogram(this Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int[] rgb = Histogram(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="channel">Channel of RGBA model</param>
        /// <returns>Array</returns>
        public static int[] Histogram(this Bitmap Data, RGBA channel)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int[] rgb = Histogram(bmData, channel);
            BitmapFormat.Unlock(Data, bmData);
            return rgb;
        }
        #endregion

        #region Density
        /// <summary>
        /// Gets an array of values of the density function.
        /// </summary>
        /// <param name="H">Histogram</param>
        /// <returns>Array</returns>
        public static int[] CDF(int[] H)
        {
            int length = H.Length;
            int[] cdf = new int[length];
            cdf[0] = H[0];

            for (int i = 1; i < length; i++)
            {
                cdf[i] = H[i] + cdf[i - 1];
            }
            return cdf;
        }
        /// <summary>
        /// Gets an array of equalized histogram values by recalculating the brightness density function.
        /// </summary>
        /// <param name="H">Histogram</param>
        /// <returns>Array</returns>
        public static float[] Equalize(int[] H)
        {
            // CDF calculating
            int length = H.Length;
            int[] cdf = Statistics.CDF(H);
            _ = Statistics.Min(cdf, out int min);
            _ = Statistics.Max(cdf, out int max);

            // table
            float[] table = new float[length];
            float factor = cdf[max] - cdf[min];

            // scaling
            for (int i = 0; i < length; i++)
            {
                table[i] = (cdf[i] - cdf[min]) / factor;
            }

            return table;
        }
        #endregion

        #region Otsu's threshold
        /// <summary>
        /// Calculates the optimal threshold using the Otsu method for the original bitmap.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Integer number</returns>
        public static int OtsuThreshold(this BitmapData bmData)
        {
            float[] v = new float[256];
            int[] h = Statistics.Histogram(bmData);
            float diff, w1, w2;

            for (int k = 1; k < 256; k++)
            {
                w1 = Omega(0, k, h);
                w2 = Omega(k + 1, 255, h);
                diff = Mu(0, k, h) * w2 - Mu(k + 1, 255, h) * w1;
                v[k] = (diff * diff) / (w1 * w2);
            }

            _ = v.Max(out int index);
            return index;
        }
        /// <summary>
        /// Calculates the optimal threshold using the Otsu method for the original bitmap.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Integer number</returns>
        public static int OtsuThreshold(this Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int threshold = OtsuThreshold(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return threshold;
        }
        #endregion

        #region SIS threshold
        /// <summary>
        /// Calculates the optimal threshold for the original bitmap.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Integer number</returns>
        public unsafe static int SISThreshold(this BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int width1 = width - 1, height1 = height - 1, offset = stride - width;
            double ex, ey, weight, weightTotal = 0, total = 0;
            byte* p = (byte*)bmData.Scan0.ToPointer(); p += stride;
            int y, x;

            for (y = 1; y < height1; y++)
            {
                p++;

                for (x = 1; x < width1; x++, p++)
                {
                    ex = Math.Abs(p[1] - p[-1]);
                    ey = Math.Abs(p[stride] - p[-stride]);
                    weight = (ex > ey) ? ex : ey;
                    weightTotal += weight;
                    total += weight * (*p);
                }
                p += offset + 1;
            }

            return (weightTotal == 0) ? 0 : (byte)(total / weightTotal);
        }
        /// <summary>
        /// Calculates the optimal threshold for the original bitmap.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Integer number</returns>
        public static int SISThreshold(this Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            int threshold = SISThreshold(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return threshold;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Omega.
        /// </summary>
        /// <param name="init">Init</param>
        /// <param name="end">End</param>
        /// <param name="h">Histogram</param>
        /// <returns>Double precision floating point number</returns>
        private static int Omega(int init, int end, int[] h)
        {
            int sum = 0, i;

            for (i = init; i <= end; i++)
            {
                sum += h[i];
            }

            return sum;
        }
        /// <summary>
        /// Mean.
        /// </summary>
        /// <param name="init">Init</param>
        /// <param name="end">End</param>
        /// <param name="h">Histogram</param>
        /// <returns>Double precision floating point number</returns>
        private static int Mu(int init, int end, int[] h)
        {
            int sum = 0, i;

            for (i = init; i <= end; i++)
            {
                sum += i * h[i];
            }

            return sum;
        }
        #endregion

        #region AForge
        /// <summary>
        /// Gets the index of the maximum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        public static int Max(int[] data)
        {
            return Max(data, out _);
        }
        /// <summary>
        /// Gets the index of the maximum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <param name="index">Max index</param>
        /// <returns>Integer number</returns>
        public static int Max(int[] data, out int index)
        {
            index = 0;
            int length = data.Length - 1, i;
            int maximum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] > maximum)
                {
                    maximum = data[i];
                    index = i;
                }
            }

            return maximum;
        }
        /// <summary>
        /// Gets the index of the minimum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        public static int Min(int[] data)
        {
            return Min(data, out _);
        }
        /// <summary>
        /// Gets the index of the minimum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <param name="index">Min index</param>
        /// <returns>Tuple of integer numbers</returns>
        public static int Min(int[] data, out int index)
        {
            index = 0;
            int length = data.Length - 1, i;
            int minimum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] < minimum)
                {
                    minimum = data[i];
                    index = i;
                }
            }

            return minimum;
        }
        /// <summary>
        /// Returns the summary of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static float Sum(int[] values)
        {
            int length = values.Length;
            float sum = 0;

            // for all values
            for (int i = 0; i < length; i++)
            {
                sum += values[i];
            }
            return sum;
        }
        /// <summary>
        /// Returns the mean of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static float Mean(int[] values)
        {
            int hits;
            long total = 0;
            float mean = 0;
            int length = values.Length;

            // for all values
            for (int i = 0, n = length; i < n; i++)
            {
                hits = values[i];
                // accumulate mean
                mean += i * hits;
                // accumalate total
                total += hits;
            }
            return (total == 0) ? 0 : mean / total;
        }
        /// <summary>
        /// Returns the standart deviation of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static float StdDev(int[] values)
        {
            return StdDev(values, Mean(values));
        }
        /// <summary>
        /// Returns the standart deviation of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <param name="mean">Mean</param>
        /// <returns>Value</returns>
        public static float StdDev(int[] values, float mean)
        {
            float stddev = 0;
            float diff;
            int hits;
            int total = 0;

            // for all values
            for (int i = 0, n = values.Length; i < n; i++)
            {
                hits = values[i];
                diff = (float)i - mean;
                // accumulate std.dev.
                stddev += diff * diff * hits;
                // accumalate total
                total += hits;
            }

            return (total == 0) ? 0 : Maths.Sqrt(stddev / total);
        }
        /// <summary>
        /// Returns the median of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static int Median(int[] values)
        {
            int total = 0, n = values.Length;

            // for all values
            for (int i = 0; i < n; i++)
            {
                // accumalate total
                total += values[i];
            }

            int halfTotal = total / 2;
            int median = 0, v = 0;

            // find median value
            for (; median < n; median++)
            {
                v += values[median];
                if (v >= halfTotal)
                    break;
            }

            return median;
        }
        /// <summary>
        /// Returns range of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <param name="percent">Percent</param>
        /// <returns>Value</returns>
        public static RangeInt GetRange(int[] values, float percent)
        {
            int total = 0, n = values.Length;

            // for all values
            for (int i = 0; i < n; i++)
            {
                // accumalate total
                total += values[i];
            }

            int min, max, hits;
            int h = (int)(total * (percent + (1 - percent) / 2));

            // get range min value
            for (min = 0, hits = total; min < n; min++)
            {
                hits -= values[min];
                if (hits < h)
                    break;
            }
            // get range max value
            for (max = n - 1, hits = total; max >= 0; max--)
            {
                hits -= values[max];
                if (hits < h)
                    break;
            }
            return new RangeInt(min, max);
        }
        /// <summary>
        /// Returns entropy of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static float Entropy(int[] values)
        {
            int n = values.Length;
            int total = 0;
            float entropy = 0;
            float p;

            // calculate total amount of hits
            for (int i = 0; i < n; i++)
            {
                total += values[i];
            }

            if (total != 0)
            {
                // for all values
                for (int i = 0; i < n; i++)
                {
                    // get item's probability
                    p = (float)values[i] / total;
                    // calculate entropy
                    if (p != 0)
                        entropy += -p * Maths.Log2(p);
                }
            }
            return entropy;
        }
        /// <summary>
        /// Returns mode of a histogram.
        /// </summary>
        /// <param name="values">Histogram</param>
        /// <returns>Value</returns>
        public static int Mode(int[] values)
        {
            int length = values.Length;
            int mode = 0, curMax = 0;

            for (int i = 0; i < length; i++)
            {
                if (values[i] > curMax)
                {
                    curMax = values[i];
                    mode = i;
                }
            }

            return mode;
        }
        #endregion
    }
}
