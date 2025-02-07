using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the color filter based on the YUV structure.
    /// </summary>
    [Serializable]
    public class YUVPhotoFilter : IBitmapFilter
    {
        #region Private data
        private float s;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color filter based on the YUV structure.
        /// </summary>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public YUVPhotoFilter(Color color, float strength = 0.5f)
        {
            Color = color; Strength = strength;
        }
        /// <summary>
        /// Initializes the color filter based on the YUV structure.
        /// </summary>
        public YUVPhotoFilter()
        {
            Color = Color.White; Strength = 0.5f;
        }
        /// <summary>
        /// Gets or sets the filter color.
        /// </summary>
        public Color Color
        {
            get
            {
                return color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets the filter strength [0, 1].
        /// </summary>
        public float Strength
        {
            get
            {
                return this.s;
            }
            set
            {
                this.s = Maths.Float(value);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            float z = 1 - s;

            Parallel.For(0, height, y =>
            {
                YUV nYUV, iYUV; RGB rgb;
                int nR, nG, nB, iR, iG, iB;
                int x, ystride, k, luma;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    iR = p[k + 2]; iG = p[k + 1]; iB = p[k];

                    luma = RGB.HDTV(iR, iG, iB);
                    nYUV = YUV.FromRGB(luma, luma, luma);
                    iYUV = YUV.FromRGB(color.R, color.G, color.B);
                    rgb = AddColor(nYUV, iYUV).ToRGB;
                    nR = rgb.Red; nG = rgb.Green; nB = rgb.Blue;

                    p[k + 2] = Maths.Byte(nR * s + iR * z);
                    p[k + 1] = Maths.Byte(nG * s + iG * z);
                    p[k] = Maths.Byte(nB * s + iB * z);
                }
            }
            );

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
        /// <summary>
        /// Checks if the color is a shade of gray.
        /// </summary>
        /// <param name="color">Color</param>
        /// <returns>Boolean</returns>
        public static bool IsGrayColor(Color color)
        {
            if (color.R == color.G && color.G == color.B)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Blend two colors in YUV space.
        /// </summary>
        /// <param name="yuv1">First color</param>
        /// <param name="yuv2">Second color</param>
        /// <returns>YUV</returns>
        public static YUV AddColor(YUV yuv1, YUV yuv2)
        {
            return new YUV(yuv1.Y, yuv1.U + yuv2.U, yuv1.V + yuv2.V);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Initializes the sepia filter.
        /// </summary>
        public static YUVPhotoFilter Sepia
        {
            get
            {
                return new YUVPhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Initializes the orange filter.
        /// </summary>
        public static YUVPhotoFilter Orange
        {
            get
            {
                return new YUVPhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Initializes the yellow filter.
        /// </summary>
        public static YUVPhotoFilter Yellow
        {
            get
            {
                return new YUVPhotoFilter(Color.Yellow);
            }
        }
        #endregion
    }
}
