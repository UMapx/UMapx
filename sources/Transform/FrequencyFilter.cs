using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the frequency filter.
    /// </summary>
    [Serializable]
    public class FrequencyFilter : IFilter
    {
        #region Private data
        private RangeInt frequencyRange = new RangeInt(0, 64);
        #endregion

        #region Filter components
        /// <summary>
        /// Gets or sets the frequency range.
        /// </summary>
        public RangeInt FrequencyRange
        {
            get
            {
                return frequencyRange;
            }
            set
            {
                frequencyRange = value;
            }
        }
        /// <summary>
        /// Initializes the frequency filter.
        /// </summary>
        public FrequencyFilter() { }
        /// <summary>
        /// Initializes the frequency filter.
        /// </summary>
        /// <param name="frequencyRange">Frequency range</param>
        public FrequencyFilter(RangeInt frequencyRange)
        {
            this.frequencyRange = frequencyRange;
        }
        /// <summary>
        /// Initializes the frequency filter.
        /// </summary>
        /// <param name="min">Minimum frequency</param>
        /// <param name="max">Maximum frequency</param>
        public FrequencyFilter(int min, int max)
        {
            this.frequencyRange = new RangeInt(min, max);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            int height = data.GetLength(0);
            int width = data.GetLength(1);

            int hh = height >> 1;
            int hw = width >> 1;

            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, j, x, y, d;

            for (i = 0; i < height; i++)
            {
                y = i - hh;

                for (j = 0; j < width; j++)
                {
                    x = j - hw;
                    d = (int)Math.Sqrt(x * x + y * y);

                    if ((d > max) || (d < min))
                    {
                        data[i, j] = 0;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            int length = data.Length;
            int hh = length >> 1;

            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, d;

            for (i = 0; i < length; i++)
            {
                d = i - hh;

                if ((d > max) || (d < min))
                {
                    data[i] = 0;
                }
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex[,] data)
        {
            int height = data.GetLength(0);
            int width = data.GetLength(1);

            int hh = height >> 1;
            int hw = width >> 1;

            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, j, x, y, d;

            for (i = 0; i < height; i++)
            {
                y = i - hh;

                for (j = 0; j < width; j++)
                {
                    x = j - hw;
                    d = (int)Math.Sqrt(x * x + y * y);

                    if ((d > max) || (d < min))
                    {
                        data[i, j].Real = 0;
                        data[i, j].Imag = 0;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex[] data)
        {
            int length = data.Length;
            int hh = length >> 1;

            int min = frequencyRange.Min;
            int max = frequencyRange.Max;

            int i, d;

            for (i = 0; i < length; i++)
            {
                d = i - hh;

                if ((d > max) || (d < min))
                {
                    data[i].Real = 0;
                    data[i].Imag = 0;
                }
            }

            return;
        }
        #endregion
    }
}
