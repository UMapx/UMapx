using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the RGB filter.
    /// </summary>
    [Serializable]
    public class RGBFilter : IBitmapFilter
    {
        #region Private data
        private int red;
        private int green;
        private int blue;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the RGB filter.
        /// </summary>
        /// <param name="red">Red [-255, 255]</param>
        /// <param name="green">Green [-255, 255]</param>
        /// <param name="blue">Blue [-255, 255]</param>
        public RGBFilter(int red, int green, int blue)
        {
            Red = red;
            Green = green;
            Blue = blue;
        }
        /// <summary>
        /// Initializes the RGB filter.
        /// </summary>
        public RGBFilter()
        {
            new RGBFilter(0, 0, 0);
        }
        /// <summary>
        /// Red [-255, 255].
        /// </summary>
        public int Red
        {
            get
            {
                return this.red;
            }
            set
            {
                this.red = value;
            }
        }
        /// <summary>
        /// Green [-255, 255].
        /// </summary>
        public int Green
        {
            get
            {
                return this.green;
            }
            set
            {
                this.green = value;
            }
        }
        /// <summary>
        /// Blue [-255, 255].
        /// </summary>
        public int Blue
        {
            get
            {
                return this.blue;
            }
            set
            {
                this.blue = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(p[2] + red);
                    p[1] = Maths.Byte(p[1] + green);
                    p[0] = Maths.Byte(p[0] + blue);
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
