using System;
using SkiaDrawing;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the mask color filter.
    /// </summary>
    [Serializable]
    public class MaskColorFilter : IBitmapFilter2
    {
        #region Filter components

        /// <summary>
        /// Initializes the mask color filter.
        /// </summary>
        /// <param name="color">Mask color</param>
        public MaskColorFilter(Color color)
        {
            Color = color;
        }

        /// <summary>
        /// Gets or sets color.
        /// </summary>
        public Color Color { get; set; }

        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (x = 0; x < width; x++)
            {
                for (y = 0; y < height; y++, p += 4, pSrc += 4)
                {
                    var average = (byte)RGB.Average(pSrc[2], pSrc[1], pSrc[0]);

                    if (average > 0)
                    {
                        p[2] = Merge.merge(Color.R, p[2], Color.A, 255 - Color.A);
                        p[1] = Merge.merge(Color.G, p[1], Color.A, 255 - Color.A);
                        p[0] = Merge.merge(Color.B, p[0], Color.A, 255 - Color.A);
                    }
                }
            }
            return;
        }

        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
        }

        #endregion
    }
}
