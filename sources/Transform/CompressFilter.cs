using System;
using UMapx.Core;

namespace UMapx.Transform
{
    /// <summary>
    /// Defines the compression filter by threshold value.
    /// </summary>
    [Serializable]
    public class CompressFilter : IFilter
    {
        #region Private data
        /// <summary>
        /// Threshold value.
        /// </summary>
        private double threshold = 0;
        /// <summary>
        /// Compress type.
        /// </summary>
        private Compress compresstype = Compress.Abs;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the compression filter by threshold value.
        /// </summary>
        public CompressFilter() { }
        /// <summary>
        /// Initializes the compression filter by threshold value.
        /// </summary>
        /// <param name="threshold">Threshold value</param>
        /// <param name="compresstype">Compress type</param>
        public CompressFilter(double threshold, Compress compresstype = Compress.Abs)
        {
            this.threshold = threshold;
            this.compresstype = compresstype;
        }
        /// <summary>
        /// Gets or sets the compress type.
        /// </summary>
        public Compress CompressType
        {
            get
            {
                return this.compresstype;
            }
            set
            {
                this.compresstype = value;
            }
        }
        /// <summary>
        /// Gets or sets the threshold value.
        /// </summary>
        public double Threshold
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
        public void Apply(double[] data)
        {
            int length = data.Length;
            int i;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Math.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
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
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Array</param>
        public void Apply(Complex[] data)
        {
            int length = data.Length;
            int i;

            if (this.compresstype == Compress.Abs)
            {
                for (i = 0; i < length; i++)
                {
                    if (Maths.Abs(data[i]) < threshold)
                    {
                        data[i] = 0;
                    }
                }
            }
            else if (this.compresstype == Compress.Over)
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

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(double[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.compresstype == Compress.Abs)
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
            else if (this.compresstype == Compress.Over)
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
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="data">Matrix</param>
        public void Apply(Complex[,] data)
        {
            int width = data.GetLength(1);
            int height = data.GetLength(0);
            int i, j;

            if (this.compresstype == Compress.Abs)
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
            else if (this.compresstype == Compress.Over)
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
            return;
        }
        #endregion

        #region Compress type
        /// <summary>
        /// Defines the compress type.
        /// </summary>
        public enum Compress
        {
            #region Types
            /// <summary>
            /// Absolute compression.
            /// </summary>
            Abs,
            /// <summary>
            /// ompression of values is less than threshold.
            /// </summary>
            Under,
            /// <summary>
            /// Compression of values is greater than the threshold.
            /// </summary>
            Over,
            #endregion
        }
        #endregion
    }
}
