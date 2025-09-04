using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the threshold filter.
    /// </summary>
    [Serializable]
    public class ThresholdFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Threshold value.
        /// </summary>
        private float threshold = 0;
        /// <summary>
        /// Threshold mode.
        /// </summary>
        private ThresholdMode type = ThresholdMode.Abs;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        public ThresholdFilter() { }
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        /// <param name="threshold">Threshold value</param>
        /// <param name="mode">Threshold mode</param>
        public ThresholdFilter(float threshold, ThresholdMode mode = ThresholdMode.Abs)
        {
            this.Threshold = threshold;
            this.Type = mode;
        }
        /// <summary>
        /// Gets or sets the threshold mode.
        /// </summary>
        public ThresholdMode Type
        {
            get
            {
                return this.type;
            }
            set
            {
                this.type = value;
            }
        }
        /// <summary>
        /// Gets or sets the threshold value.
        /// </summary>
        public float Threshold
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(float[] data)
        {
            int length = data.Length;
            int i;

            if (this.type == ThresholdMode.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Math.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.type == ThresholdMode.Over)
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i] > threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i] < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex32[] data)
        {
            int length = data.Length;
            int i;

            if (this.type == ThresholdMode.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Maths.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.type == ThresholdMode.Over)
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i].Real > threshold)
                    {
                        data[i].Real = 0;
                    }
                    if (data[i].Imag > threshold)
                    {
                        data[i].Imag = 0;
                    }
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    if (data[i].Real < threshold)
                    {
                        data[i].Real = 0;
                    }
                    if (data[i].Imag < threshold)
                    {
                        data[i].Imag = 0;
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(float[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.type == ThresholdMode.Abs)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (Math.Abs(data[i, j]) < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else if (this.type == ThresholdMode.Over)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j] > threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j] < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex32[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.type == ThresholdMode.Abs)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (Maths.Abs(data[i, j]) < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else if (this.type == ThresholdMode.Over)
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j].Real > threshold)
                        {
                            data[i, j] = 0;
                        }
                        if (data[i, j].Imag > threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < height; i++)
                {
                    for (j = 0; j < width; j++)
                    {
                        if (data[i, j].Real < threshold)
                        {
                            data[i, j] = 0;
                        }
                        if (data[i, j].Imag < threshold)
                        {
                            data[i, j] = 0;
                        }
                    }
                }
            }
        }
        #endregion
    }
}
