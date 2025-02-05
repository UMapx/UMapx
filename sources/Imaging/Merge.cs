using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the merge filter.
    /// </summary>
    [Serializable]
    public class Merge : IBitmapFilter2
    {
        #region Private data
        private int transparency = 255;
        private PointInt point;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        /// <param name="transparency">Transparency [0, 255]</param>
        public Merge(int transparency = 255)
        {
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        /// <param name="point">A pair of integers representing an ordered pair of X and Y coordinates</param>
        /// <param name="transparency">Transparency [0, 255]</param>
        public Merge(PointInt point, int transparency = 255)
        {
            Point = point;
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="transparency">Transparency [0, 255]</param>
        public Merge(int x, int y, int transparency = 255)
        {
            Point = new PointInt(x, y);
            Transparency = transparency;
        }
        /// <summary>
        /// Gets or sets the transparency value [0, 255].
        /// </summary>
        public int Transparency
        {
            get
            {
                return this.transparency;
            }
            set
            {
                this.transparency = value;
            }
        }
        /// <summary>
        /// Gets or sets a pair of integers representing an ordered pair of X and Y coordinates.
        /// </summary>
        public PointInt Point
        {
            get
            {
                return point;
            }
            set
            {
                point = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image:
            int srcStride = bmSrc.Stride, dstStride = bmData.Stride;

            // source pixel's coordinates from rectangle:
            int startX = (int)Maths.Max(point.X, 0);
            int startY = (int)Maths.Max(point.Y, 0);
            int endX = (int)Maths.Min(bmSrc.Width, bmData.Width - startX);
            int endY = (int)Maths.Min(bmSrc.Height, bmData.Height - startY);

            // do the job:
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // transparency:
            double t = transparency / 255.0;

            // do job
            Parallel.For(0, endY, y =>
            {
                // pixel offsets:
                int x, k, l;
                int a0, a1;

                // for each pixel:
                for (x = 0; x < endX; x++)
                {
                    // local data:
                    k = dstStride * (y + startY) + 4 * (x + startX);
                    l = srcStride * (y) + 4 * (x);

                    // calculating transparency:
                    a1 = (int)(src[l + 3] * t);
                    a0 = 255 - a1;

                    // applying filter:
                    dst[k + 0] = merge(dst[k + 0], src[l + 0], a0, a1);
                    dst[k + 1] = merge(dst[k + 1], src[l + 1], a0, a1);
                    dst[k + 2] = merge(dst[k + 2], src[l + 2], a0, a1);
                }
            }
            );

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

        #region Merging function components

        /// <summary>
        /// Merge function.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="a0"></param>
        /// <param name="a1"></param>
        /// <returns></returns>
        internal static byte merge(byte x, byte y, int a0, int a1)
        {
            return Maths.Byte((x * a0 + y * a1) / 255);
        }

        #endregion
    }
}
