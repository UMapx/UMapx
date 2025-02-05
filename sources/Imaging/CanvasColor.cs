using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the color canvas class.
    /// </summary>
    [Serializable]
    public class CanvasColor : ICanvas
    {
        #region Private data
        private Color color = Color.White;
        private Bitmap bitmap = new Bitmap(256, 256);
        private int width = 256, height = 256;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the color canvas class.
        /// </summary>
        public CanvasColor() { }
        /// <summary>
        /// Initializes the color canvas class.
        /// </summary>
        /// <param name="width">Canvas width</param>
        /// <param name="height">Canvas height</param>
        /// <param name="color">Color</param>
        public CanvasColor(int width, int height, Color color)
        {
            Width = width; Height = height; Color = color;
        }
        /// <summary>
        /// Gets or sets canvas color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets the width of the canvas.
        /// </summary>
        public int Width
        {
            get
            {
                return this.width;
            }
            set
            {
                this.width = value;
            }
        }
        /// <summary>
        /// Gets or sets the height of the canvas.
        /// </summary>
        public int Height
        {
            get
            {
                return this.height;
            }
            set
            {
                this.height = value;
            }
        }
        /// <summary>
        /// Creates canvas.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap Create()
        {
            bitmap = new Bitmap(width, height);
            Graphics graphics = Graphics.FromImage(bitmap);
            SolidBrush brush = new SolidBrush(color);
            graphics.FillRectangle(brush, 0, 0, width, height);
            graphics.Dispose();
            brush.Dispose();
            return bitmap;
        }
        #endregion
    }
}
