using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the point addition filter.
    /// </summary>
    [Serializable]
    public class PointAddition : Rebuilder, IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Point matrix.
        /// </summary>
        protected PointInt[,] points = new PointInt[0, 0];
        /// <summary>
        /// Image width.
        /// </summary>
        protected int width;
        /// <summary>
        /// Image height.
        /// </summary>
        protected int height;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the point addition filter.
        /// </summary>
        /// <param name="points">Array of ordered pairs of X and Y</param>
        public PointAddition(PointInt[,] points)
        {
            Points = points;
        }
        /// <summary>
        /// Initializes the point addition filter.
        /// </summary>
        public PointAddition() { }
        /// <summary>
        /// Gets or sets point matrix.
        /// </summary>
        public PointInt[,] Points
        {
            get
            {
                return this.points;
            }
            set
            {
                this.points = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get image sizes
            this.width = bmData.Width;
            this.height = bmData.Height;

            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            // do job
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                int x, y, i, c, x1, y1, jstride, k;
                jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    x = points[i, j].X;
                    y = points[i, j].Y;

                    x1 = i + x; y1 = j + y;

                    if (y1 >= 0 && y1 < height && x1 >= 0 && x1 < width)
                    {
                        c = (j + y) * stride + (i + x) * 4;

                        p[k] = pSrc[c]; // Blue
                        p[k + 1] = pSrc[c + 1]; // Green
                        p[k + 2] = pSrc[c + 2]; // Red
                        p[k + 3] = pSrc[c + 3]; // Alpha
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
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap current = BitmapFormat.Bitmap(bmData);
            Bitmap Src = (Bitmap)current.Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
            current.Dispose();
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
            return;
        }
        #endregion
    }
}
