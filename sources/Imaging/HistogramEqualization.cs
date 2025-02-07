using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the global histogram equalization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cromwell-intl.com/3d/histogram/
    /// </remarks>
    /// </summary>
    [Serializable]
    public class HistogramEqualization : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the global histogram equalization filter.
        /// </summary>
        public HistogramEqualization() { }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int[] H = Statistics.Histogram(bmData);
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            float[] table = Statistics.Equalize(H);
            float length = table.Length - 1;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(table[p[2]] * length);
                    p[1] = Maths.Byte(table[p[1]] * length);
                    p[0] = Maths.Byte(table[p[0]] * length);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
}
