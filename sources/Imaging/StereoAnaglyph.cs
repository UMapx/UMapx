using System;
using System.Drawing;
using System.Drawing.Imaging;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the stereo effect filter for a pair of images.
    /// </summary>
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
    /// </remarks>
    [Serializable]
    public class StereoAnaglyph : IBitmapFilter2
    {
        #region Private data
        private AnaglyphMode algorithm = AnaglyphMode.Gray;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the stereo effect filter for a pair of images.
        /// </summary>
        /// <param name="algorithm">Algorithm</param>
        public StereoAnaglyph(AnaglyphMode algorithm)
        {
            this.algorithm = algorithm;
        }
        /// <summary>
        /// Gets or sets the algorithm.
        /// </summary>
        public AnaglyphMode Algorithm
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
                case AnaglyphMode.True:
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

                case AnaglyphMode.Gray:
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

                case AnaglyphMode.Color:
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

                case AnaglyphMode.HalfColor:
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

                case AnaglyphMode.Optimized:
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
