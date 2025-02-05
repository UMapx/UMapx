using System;
using UMapx.Core;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Uses for editing and transforming images using GDI+.
    /// </summary>
    public static class BitmapTransform
    {
        #region Rotate
        /// <summary>
        /// Rotates bitmap by rotation value.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="rotation">Rotation</param>
        /// <returns></returns>
        public static Bitmap Rotate(this Bitmap b, RotationMode rotation)
        {
            var clone = (Bitmap)b.Clone();
            switch (rotation)
            {
                case RotationMode.R0:
                    break;
                case RotationMode.R90:
                    clone.RotateFlip(RotateFlipType.Rotate90FlipNone);
                    break;
                case RotationMode.R180:
                    clone.RotateFlip(RotateFlipType.Rotate180FlipNone);
                    break;
                case RotationMode.R270:
                    clone.RotateFlip(RotateFlipType.Rotate270FlipNone);
                    break;
                default:
                    break;
            }
            return clone;
        }
        /// <summary>
        /// Rotates the bitmap by the specified angle. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="angle">Angle</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Rotate(this Bitmap b, float angle)
        {
            return Rotate(b, new PointFloat(b.Width / 2.0f, b.Height / 2.0f), angle);
        }
        
        #region Private
        /// <summary>
        /// Rotates the bitmap by the specified angle. 
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="point">Float point</param>
        /// <param name="angle">Angle in degrees</param>
        /// <returns>Bitmap</returns>
        private static Bitmap Rotate(this Bitmap b, PointFloat point, float angle)
        {
            Bitmap bmp = new Bitmap(b.Width, b.Height);
            bmp.SetResolution(b.HorizontalResolution, b.VerticalResolution);
            Graphics graphics = Graphics.FromImage(bmp);
            graphics.InterpolationMode = SkiaDrawing.InterpolationMode.HighQualityBilinear;
            graphics.TranslateTransform(point.X, point.Y);
            graphics.RotateTransform(angle);
            graphics.TranslateTransform(-point.X, -point.Y);
            graphics.DrawImage(b, new PointF(0, 0));
            graphics.Dispose();
            return bmp;
        }
        #endregion

        /// <summary>
        /// Returns rotated image.
        /// </summary>
        /// <param name="image">Bitmap</param>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Rotate(this Bitmap image, float angle, Color color)
        {
            // create an empty Bitmap image
            Bitmap bmp = new Bitmap(image.Width, image.Height);

            // turn the Bitmap into a Graphics object
            using var g = Graphics.FromImage(bmp);
            g.Clear(color);

            // now we set the rotation point to the center of our image
            g.TranslateTransform((float)bmp.Width / 2, (float)bmp.Height / 2);

            // now rotate the image
            g.RotateTransform(angle);
            g.TranslateTransform(-(float)bmp.Width / 2, -(float)bmp.Height / 2);

            // set the InterpolationMode to HighQualityBicubic so to ensure a high
            // quality image once it is transformed to the specified size
            g.InterpolationMode = SkiaDrawing.InterpolationMode.HighQualityBicubic;

            // now draw our new image onto the graphics object
            g.DrawImage(image, new Point(0, 0));

            //return the image
            return bmp;
        }
        #endregion

        #region Flip
        /// <summary>
        /// Flips bitmap by direction.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="direction">Direction</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Flip(this Bitmap b, Direction direction)
        {
            var clone = (Bitmap)b.Clone();
            switch (direction)
            {
                case Direction.Horizontal:
                    clone.RotateFlip(RotateFlipType.RotateNoneFlipX);
                    break;
                case Direction.Vertical:
                    clone.RotateFlip(RotateFlipType.RotateNoneFlipY);
                    break;
                case Direction.Both:
                    clone.RotateFlip(RotateFlipType.RotateNoneFlipXY);
                    break;
                default:
                    break;
            }
            
            return clone;
        }
        #endregion

        #region Crop
        /// <summary>
        /// Returns cropped image.
        /// </summary>
        /// <param name="image">Bitmap</param>
        /// <param name="rectangle">Rectangle</param>
        /// <param name="clamp">Clamp crop or not</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Crop(this Bitmap image, Rectangle rectangle, bool clamp = true)
        {
            // image params
            int width = image.Width;
            int height = image.Height;

            // check section params
            int x = clamp ? Maths.Range(rectangle.X, 0, width) : rectangle.X;
            int y = clamp ? Maths.Range(rectangle.Y, 0, height) : rectangle.Y;
            int w = clamp ? Maths.Range(rectangle.Width, 0, width - x) : rectangle.Width;
            int h = clamp ? Maths.Range(rectangle.Height, 0, height - y) : rectangle.Height;

            // exception
            if (x == 0 &&
                y == 0 &&
                w == 0 &&
                h == 0) return image;

            // crop image to rectangle section
            var src = new Rectangle(x, y, w, h);
            var bitmap = new Bitmap(src.Width, src.Height);
            var dest = new Rectangle(0, 0, bitmap.Width, bitmap.Height);
            using var g = Graphics.FromImage(bitmap);
            g.DrawImage(image, dest, src, GraphicsUnit.Pixel);

            return bitmap;
        }
        #endregion

        #region Resize
        /// <summary>
        /// Returns resized image.
        /// </summary>
        /// <param name="image">Bitmap</param>
        /// <param name="size">Size</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Resize(this Bitmap image, Size size)
        {
            return new Bitmap(image, size.Width, size.Height);
        }

        /// <summary>
        /// Returns resized image with preserved proportions.
        /// </summary>
        /// <param name="image">Bitmap</param>
        /// <param name="size">Size</param>
        /// <param name="color">Border color</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ResizePreserved(this Bitmap image, Size size, Color color)
        {
            // size
            int width = image.Width;
            int height = image.Height;
            int max = Math.Max(width, height);

            //  borders
            var rectangle = new Rectangle(
                (max - width) / 2,
                (max - height) / 2,
                width,
                height);

            // drawing
            using var background = new Bitmap(max, max);

            using var g = Graphics.FromImage(background);
            g.Clear(color);
            g.DrawImage(image, rectangle);

            return BitmapTransform.Resize(background, size);
        }

        /// <summary>
        /// Returns resized image with preserved proportions.
        /// </summary>
        /// <param name="image">Bitmap</param>
        /// <param name="size">Size</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ResizePreserved(this Bitmap image, Size size)
        {
            // size
            int width = size.Width;
            int height = size.Height;
            int max = Math.Max(width, height);

            //  borders
            var rectangle = new Rectangle(
                (max - width) / 2,
                (max - height) / 2,
                width,
                height);

            // drawing
            using var resized = BitmapTransform.Resize(image, new Size(max, max));
            return BitmapTransform.Crop(resized, rectangle);
        }

        #endregion

        #region Shift
        /// <summary>
        /// Shifts the bitmap.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="h">The number of positions to which a shift in height occurs</param>
        /// <param name="w">The number of positions by which the shift occurs in width</param>
        /// <returns></returns>
        public static Bitmap Shift(this Bitmap b, int w, int h)
        {
            if (w == 0)
            {
                if (h == 0)
                {
                    return b;
                }
                return ShiftY(b, h);
            }
            else
            {
                if (h == 0)
                {
                    return ShiftX(b, w);
                }
                else
                {
                    using var c = ShiftY(b, h);
                    return ShiftX(c, w);
                }
            }
        }

        #region Private
        /// <summary>
        /// Shifts bitmap by X axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Shift value</param>
        /// <returns>Bitmap</returns>
        private static Bitmap ShiftX(Bitmap b, int value)
        {
            var bmp = (Bitmap)b.Clone();
            using var graphics = Graphics.FromImage(bmp);
            graphics.DrawImage(b, value, 0, b.Width, b.Height);
            graphics.DrawImage(b, value - b.Width, 0, b.Width, b.Height);
            return bmp;
        }
        /// <summary>
        /// Shifts bitmap by Y axis.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="value">Shift value</param>
        /// <returns>Bitmap</returns>
        private static Bitmap ShiftY(Bitmap b, int value)
        {
            var bmp = (Bitmap)b.Clone();
            using var graphics = Graphics.FromImage(bmp);
            graphics.DrawImage(b, 0, value, b.Width, b.Height);
            graphics.DrawImage(b, 0, value - b.Height, b.Width, b.Height);
            return bmp;
        }
        #endregion
        #endregion

        #region Merge
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>s
        /// <param name="background">Background bitmap</param>
        /// <param name="foreground">Foreground bitmap</param>
        /// <returns>Bitmap</returns>
        public static void Merge(this Bitmap background, Bitmap foreground)
        {
            var rectangle = new Rectangle(0, 0, foreground.Width, foreground.Height);
            Merge(background, foreground, rectangle);
        }
        /// <summary>
        /// Merges two bitmaps.
        /// </summary>
        /// <param name="background">Background image</param>
        /// <param name="foreground">Foreground image</param>
        /// <param name="rectangle">Rectangle</param>
        public static void Merge(this Bitmap background, Bitmap foreground, Rectangle rectangle)
        {
            using var graphics = Graphics.FromImage(background);
            graphics.DrawImage(foreground, rectangle);
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
            ColorMatrix cmx = new ColorMatrix
            {
                Matrix33 = value / 255.0f
            };
            ImageAttributes attributes = new ImageAttributes();
            attributes.SetColorMatrix(cmx);
            graphics.DrawImage(b, new Rectangle(0, 0, width, height), 0, 0, width, height, GraphicsUnit.Pixel, attributes);
            graphics.Dispose();
            attributes.Dispose();
            return bmp;
        }
        #endregion
    }
}
