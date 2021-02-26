using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses for editing and transforming images. 
    /// </summary>
    public static class BitmapTransform
    {
        #region Rotate
        /// <summary>
        /// Rotates the bitmap by the specified angle. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="angle">Angle in degrees</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Rotate(this Bitmap b, float angle)
        {
            return Rotate(b, new PointFloat(b.Width / 2.0f, b.Height / 2.0f), angle);
        }
        /// <summary>
        /// Rotates the bitmap by the specified angle. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="x">X</param>
        /// <param name="y">Y</param>
        /// <param name="angle">Angle in degrees</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Rotate(this Bitmap b, int x, int y, float angle)
        {
            return Rotate(b, new PointFloat(x, y), angle);
        }
        /// <summary>
        /// Rotates the bitmap by the specified angle. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="point">Float point</param>
        /// <param name="angle">Angle in degrees</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Rotate(this Bitmap b, PointFloat point, float angle)
        {
            Bitmap bmp = new Bitmap(b.Width, b.Height);
            bmp.SetResolution(b.HorizontalResolution, b.VerticalResolution);
            Graphics graphics = Graphics.FromImage(bmp);
            graphics.InterpolationMode = InterpolationMode.HighQualityBilinear;
            graphics.TranslateTransform(point.X, point.Y);
            graphics.RotateTransform(angle);
            graphics.TranslateTransform(-point.X, -point.Y);
            graphics.DrawImage(b, new PointF(0, 0));
            graphics.Dispose();
            return bmp;
        }
        #endregion

        #region Flip
        /// <summary>
        /// Flips bitmap by X axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap FlipX(this Bitmap b)
        {
            var clone = (Bitmap)b.Clone();
            clone.RotateFlip(RotateFlipType.RotateNoneFlipX);
            return clone;
        }
        /// <summary>
        /// Flips bitmap by Y axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap FlipY(this Bitmap b)
        {
            var clone = (Bitmap)b.Clone();
            clone.RotateFlip(RotateFlipType.RotateNoneFlipY);
            return clone;
        }
        #endregion

        #region Crop
        /// <summary>
        /// Crops bitmap.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="rectangle">Rectangle</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Crop(this Bitmap b, Rectangle rectangle)
        {
            return b.Clone(rectangle, b.PixelFormat);
        }
        #endregion

        #region Resize
        /// <summary>
        /// Resizes bitmap.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Resize(this Bitmap b, int width, int height)
        {
            return new Bitmap(b, width, height);
        }
        /// <summary>
        /// Resizes bitmap.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Scale value</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Resize(this Bitmap b, float value)
        {
            return Resize(b, (int)(b.Width * value), (int)(b.Height * value));
        }
        #endregion

        #region Shift
        /// <summary>
        /// Shifts bitmap by X axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Shift value</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ShiftX(this Bitmap b, int value)
        {
            Bitmap bmp = (Bitmap)b.Clone();
            Graphics graphics = Graphics.FromImage(bmp);
            graphics.DrawImage(b, value, 0, b.Width, b.Height);
            graphics.DrawImage(b, value - b.Width, 0, b.Width, b.Height);
            graphics.Dispose();
            return bmp;
        }
        /// <summary>
        /// Shifts bitmap by Y axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Shift value</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ShiftY(this Bitmap b, int value)
        {
            Bitmap bmp = (Bitmap)b.Clone();
            Graphics graphics = Graphics.FromImage(bmp);
            graphics.DrawImage(b, 0, value, b.Width, b.Height);
            graphics.DrawImage(b, 0, value - b.Height, b.Width, b.Height);
            graphics.Dispose();
            return bmp;
        }
        #endregion

        #region Transparecy
        /// <summary>
        /// Sets the transparency of the bitmap. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Transparency [0, 255]</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Transparency(this Bitmap b, int value)
        {
            int width = b.Width, height = b.Height;
            Bitmap bmp = new Bitmap(width, height);
            Graphics graphics = Graphics.FromImage(bmp);
            ColorMatrix cmx = new ColorMatrix();
            cmx.Matrix33 = value / 255.0f;
            ImageAttributes attributes = new ImageAttributes();
            attributes.SetColorMatrix(cmx);
            graphics.DrawImage(b, new Rectangle(0, 0, width, height), 0, 0, width, height, GraphicsUnit.Pixel, attributes);
            graphics.Dispose();
            attributes.Dispose();
            return bmp;
        }
        #endregion

        #region Merge
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>s
        /// <param name="b">Background bitmap</param>
        /// <param name="fb">Foreground bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Merge(this Bitmap b, Bitmap fb)
        {
            var rectangle = new Rectangle(0, 0, fb.Width, fb.Height);
            return Merge(b, fb, rectangle, 255);
        }
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>s
        /// <param name="b">Background bitmap</param>
        /// <param name="fb">Foreground bitmap</param>
        /// <param name="rectangle">Rectangle</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Merge(this Bitmap b, Bitmap fb, Rectangle rectangle)
        {
            return Merge(b, fb, rectangle, 255);
        }
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>s
        /// <param name="b">Background bitmap</param>
        /// <param name="fb">Foreground bitmap</param>
        /// <param name="transparency">Transparency value [0, 255]</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Merge(this Bitmap b, Bitmap fb, int transparency)
        {
            var rectangle = new Rectangle(0, 0, fb.Width, fb.Height);
            return Merge(b, fb, rectangle, transparency);
        }
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>s
        /// <param name="b">Background bitmap</param>
        /// <param name="fb">Foreground bitmap</param>
        /// <param name="rectangle">Rectangle</param>
        /// <param name="transparency">Transparency value [0, 255]</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Merge(this Bitmap b, Bitmap fb, Rectangle rectangle, int transparency)
        {
            var bmp = (Bitmap)b.Clone();
            fb = Resize(Transparency(fb, transparency), rectangle.Width, rectangle.Height);
            Graphics graphics = Graphics.FromImage(b);
            graphics.DrawImage(fb, rectangle);
            graphics.Dispose();
            return b;
        }
        #endregion
    }
}
