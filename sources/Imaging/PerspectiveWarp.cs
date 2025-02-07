using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines a perspective warp filter.
    /// </summary>
    [Serializable]
    public class PerspectiveWarp : IBitmapFilter2, IBitmapFilter
    {
        #region Constructor

        /// <summary>
        /// Initializes a perspective warp filter.
        /// </summary>
        public PerspectiveWarp() { }
        /// <summary>
        /// Initializes a perspective warp filter.
        /// </summary>
        /// <param name="topLeft">Top left point</param>
        /// <param name="topRight">Top right point</param>
        /// <param name="bottomLeft">Bottom left point</param>
        /// <param name="bottomRight">Bottom right point</param>
        public PerspectiveWarp(PointFloat topLeft, PointFloat topRight, PointFloat bottomLeft, PointFloat bottomRight) 
        {
            TopLeft = topLeft;
            TopRight = topRight;
            BottomLeft = bottomLeft;
            BottomRight = bottomRight;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Gets or sets top left point.
        /// </summary>
        public PointFloat TopLeft { get; set; }
        /// <summary>
        /// Gets or sets top right point.
        /// </summary>
        public PointFloat TopRight { get; set; }
        /// <summary>
        /// Gets or sets bottom left point.
        /// </summary>
        public PointFloat BottomLeft { get;set; }
        /// <summary>
        /// Gets or sets bottom right point.
        /// </summary>
        public PointFloat BottomRight { get; set; }

        #endregion

        #region Methods

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
            return;
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
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            var width = bmData.Width;
            var height = bmData.Height;
            var rectangle = new Rectangle(0, 0, width, height);

            var transform = Float3x3.Perspective(rectangle, TopLeft, TopRight, BottomLeft, BottomRight);
            transform = transform.Invert();

            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            var stride = bmData.Stride;

            for (var x = 0; x < width; x++)
            {
                for (var y = 0; y < height; y++)
                {
                    var tileCoord = transform.TransformPoint(new PointFloat(x, y));

                    tileCoord.X = tileCoord.X < 0 || tileCoord.X >= width ? 0 : tileCoord.X;
                    tileCoord.Y = tileCoord.Y < 0 || tileCoord.Y >= height ? 0 : tileCoord.Y;

                    var k = (int)tileCoord.Y * stride + (int)tileCoord.X * 4;
                    var j = y * stride + x * 4;

                    p[j + 0] = pSrc[k + 0];
                    p[j + 1] = pSrc[k + 1];
                    p[j + 2] = pSrc[k + 2];
                    p[j + 3] = pSrc[k + 3];
                }
            }
        }

        #endregion
    }
}
