using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the CMYK filter.
    /// </summary>
    [Serializable]
    public class CMYKFilter : IBitmapFilter
    {
        #region Private data
        private float cyan;
        private float magenta;
        private float yellow;
        private float keycolor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the CMYK filter.
        /// </summary>
        /// <param name="cyan">Cyan [-1, 1]</param>
        /// <param name="magenta">Magenta [-1, 1]</param>
        /// <param name="yellow">Yellow [-1, 1]</param>
        /// <param name="keycolor">Keycolor [-1, 1]</param>
        public CMYKFilter(float cyan, float magenta, float yellow, float keycolor)
        {
            Cyan = cyan;
            Magenta = magenta;
            Yellow = yellow;
            Keycolor = keycolor;
        }
        /// <summary>
        /// Initializes the CMYK filter.
        /// </summary>
        public CMYKFilter()
        {
            new CMYKFilter(0, 0, 0, 0);
        }
        /// <summary>
        /// Cyan [-1, 1].
        /// </summary>
        public float Cyan
        {
            get
            {
                return this.cyan;
            }
            set
            {
                this.cyan = value;
            }
        }
        /// <summary>
        /// Magenta [-1, 1].
        /// </summary>
        public float Magenta
        {
            get
            {
                return this.magenta;
            }
            set
            {
                this.magenta = value;
            }
        }
        /// <summary>
        /// Yellow [-1, 1].
        /// </summary>
        public float Yellow
        {
            get
            {
                return this.yellow;
            }
            set
            {
                this.yellow = value;
            }
        }
        /// <summary>
        /// Keycolor [-1, 1].
        /// </summary>
        public float Keycolor
        {
            get
            {
                return this.keycolor;
            }
            set
            {
                this.keycolor = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                CMYK cmyk; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    cmyk = CMYK.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    cmyk.Cyan += cyan;
                    cmyk.Magenta += magenta;
                    cmyk.Yellow += yellow;
                    cmyk.Keycolor += keycolor;
                    rgb = cmyk.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
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
            return;
        }
        #endregion
    }
}
