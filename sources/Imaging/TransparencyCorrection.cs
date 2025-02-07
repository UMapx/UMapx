using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the transparency correction filter.
    /// </summary>
    [Serializable]
    public class TransparencyCorrection : Rebuilder, IBitmapFilter
    {
        #region Private data
        private float[] values;
        private float transparency;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the transparency correction filter.
        /// </summary>
        /// <param name="transparency">Transparency [-1, 1]</param>
        public TransparencyCorrection(float transparency)
        {
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the transparency correction filter.
        /// </summary>
        public TransparencyCorrection()
        {
            Transparency = 0;
        }
        /// <summary>
        /// Gets or sets the transparency value [-1, 1].
        /// </summary>
        public float Transparency
        {
            get
            {
                return this.transparency;
            }
            set
            {
                this.transparency = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.transparency, 256);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, height = bmData.Height, width = bmData.Width;
            float length = values.Length - 1;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[3] = Maths.Byte(values[p[3]] * length);
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
