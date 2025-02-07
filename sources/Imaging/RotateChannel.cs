using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the channel rotation filter.
    /// </summary>
    [Serializable]
    public class RotateChannel : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the channel rotation filter.
        /// </summary>
        public RotateChannel() { }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte red, green, blue;
            int y, x, height = bmData.Height, width = bmData.Width;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    blue = p[0];
                    green = p[1];
                    red = p[2];

                    p[0] = red;
                    p[1] = blue;
                    p[2] = green;
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
        }
        #endregion
    }
}
