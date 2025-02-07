using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the stereo effect filter for a pair of images.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
    /// </remarks>
    /// </summary>
    [Serializable]
    public class StereoAnaglyph : IBitmapFilter2
    {
        #region Private data
        private Anaglyph algorithm = Anaglyph.Gray;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the stereo effect filter for a pair of images.
        /// </summary>
        /// <param name="algorithm">Algorithm</param>
        public StereoAnaglyph(Anaglyph algorithm)
        {
            this.algorithm = algorithm;
        }
        /// <summary>
        /// Gets or sets the algorithm.
        /// </summary>
        public Anaglyph Algorithm
        {
            get { return algorithm; }
            set { algorithm = value; }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width, height = bmData.Height;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int x, y;

            switch (algorithm)
            {
                case Anaglyph.True:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = 0;
                            p[0] = (byte)(pSrc[2] * 0.299 + pSrc[1] * 0.587 + pSrc[0] * 0.114);
                        }
                    }
                    break;

                case Anaglyph.Gray:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = (byte)(pSrc[2] * 0.299 + pSrc[1] * 0.587 + pSrc[0] * 0.114);
                            p[0] = p[1];
                        }
                    }
                    break;

                case Anaglyph.Color:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            // keep Red as it is and take only Green and Blue from the second image
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;

                case Anaglyph.HalfColor:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;

                case Anaglyph.Optimized:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[1] * 0.7 + p[0] * 0.3);
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;
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

        #region Enums
        /// <summary>
        /// Defines the stereo effect creation algorithm.
        /// </summary>
        /// <remarks>
        /// More information can be found on the website:
        /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
        /// </remarks>
        public enum Anaglyph
        {
            /// <summary>
            /// Creates a stereo effect for a pair of images UMapxing to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            True,
            /// <summary>
            /// Creates a stereo effect for a pair of images UMapxing to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Gray,
            /// <summary>
            /// Creates a stereo effect for a pair of images UMapxing to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=R<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Color,
            /// <summary>
            /// Creates a stereo effect for a pair of images UMapxing to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            HalfColor,
            /// <summary>
            /// Creates a stereo effect for a pair of images UMapxing to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.7*G<sub>l</sub>+0.3*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Optimized
        }
        #endregion
    }
}
