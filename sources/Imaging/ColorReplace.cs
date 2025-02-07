using System;
using SkiaDrawing;
using System.Threading.Tasks;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the color replacement filter.
    /// </summary>
    [Serializable]
    public class ColorReplace : IBitmapFilter
    {
        #region Private data
        private Color input = Color.Transparent;
        private Color output = Color.Transparent;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color replacement filter.
        /// </summary>
        /// <param name="input">Input color</param>
        /// <param name="output">Output color</param>
        public ColorReplace(Color input, Color output)
        {
            Input = input;
            Output = output;
        }
        /// <summary>
        /// Initializes the color replacement filter.
        /// </summary>
        public ColorReplace() { }
        /// <summary>
        /// Input color.
        /// </summary>
        public Color Input
        {
            get
            {
                return this.input;
            }
            set
            {
                this.input = value;
            }
        }
        /// <summary>
        /// Output color.
        /// </summary>
        public Color Output
        {
            get
            {
                return this.output;
            }
            set
            {
                this.output = value;
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

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    if ((p[k + 3] == input.A) && (p[k + 2] == input.R) && (p[k + 1] == input.G) && (p[k] == input.B))
                    {
                        p[k + 3] = output.A;
                        p[k + 2] = output.R;
                        p[k + 1] = output.G;
                        p[k + 0] = output.B;
                    }
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
        #endregion
    }
}
