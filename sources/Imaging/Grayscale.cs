using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the grayscale filter.
    /// </summary>
    [Serializable]
    public class Grayscale : IBitmapFilter
    {
        #region Private data
        private double cr;
        private double cg;
        private double cb;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grayscale filter.
        /// </summary>
        /// <param name="cr">Red</param>
        /// <param name="cg">Green</param>
        /// <param name="cb">Blue</param>
        public Grayscale(double cr, double cg, double cb)
        {
            Cr = cr;
            Cg = cg;
            Cb = cb;
        }
        /// <summary>
        /// Initializes the grayscale filter.
        /// </summary>
        public Grayscale()
        {
            Cr = 0.333f;
            Cg = 0.333f;
            Cb = 0.333f;
        }
        /// <summary>
        /// Gets or sets the red channel value.
        /// </summary>
        public double Cr
        {
            get
            {
                return this.cr;
            }
            set
            {
                this.cr = value;
            }
        }
        /// <summary>
        /// Gets or sets the green channel value.
        /// </summary>
        public double Cg
        {
            get
            {
                return this.cg;
            }
            set
            {
                this.cg = value;
            }
        }
        /// <summary>
        /// Gets or sets the blue channel value.
        /// </summary>
        public double Cb
        {
            get
            {
                return this.cb;
            }
            set
            {
                this.cb = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[0] = p[1] = p[2] = Maths.Byte(cr * p[2] + cg * p[1] + cb * p[0]);
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
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Initializes the grayscale filter (BT709).
        /// </summary>
        public static Grayscale BT709
        {
            get
            {
                return new Grayscale(0.212f, 0.715f, 0.072f);
            }
        }
        /// <summary>
        /// Initializes the grayscale filter (R-Y).
        /// </summary>
        public static Grayscale RY
        {
            get
            {
                return new Grayscale(0.5f, 0.419f, 0.081f);
            }
        }
        /// <summary>
        /// Initializes the grayscale filter Y.
        /// </summary>
        public static Grayscale Y
        {
            get
            {
                return new Grayscale(0.299f, 0.587f, 0.114f);
            }
        }
        #endregion
    }
}
