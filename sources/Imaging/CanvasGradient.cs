using System;
using SkiaDrawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the gradient canvas class.
    /// </summary>
    [Serializable]
    public class CanvasGradient : ICanvas
    {
        #region Private data
        private Color color1 = Color.White;
        private Color color2 = Color.Black;
        private Bitmap bitmap = new Bitmap(256, 256);
        private int width = 256, height = 256;
        private double angle = 0;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the gradient canvas class.
        /// </summary>
        public CanvasGradient() { }
        /// <summary>
        /// Initializes the gradient canvas class.
        /// </summary>
        /// <param name="width">Canvas width</param>
        /// <param name="height">Canvas height</param>
        /// <param name="angle">Angle</param>
        /// <param name="color1">First color</param>
        /// <param name="color2">Second color</param>
        public CanvasGradient(int width, int height, double angle, Color color1, Color color2)
        {
            Width = width; Height = height; Angle = angle; Color1 = color1; Color2 = color2;
        }
        /// <summary>
        /// Gets or sets the first color.
        /// </summary>
        public Color Color1
        {
            get
            {
                return this.color1;
            }
            set
            {
                this.color1 = value;
            }
        }
        /// <summary>
        /// Gets or sets the second color.
        /// </summary>
        public Color Color2
        {
            get
            {
                return this.color2;
            }
            set
            {
                this.color2 = value;
            }
        }
        /// <summary>
        /// Gets or sets the angle value.
        /// </summary>
        public double Angle
        {
            get
            {
                return this.angle;
            }
            set
            {
                this.angle = value;
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
            Rectangle rectangle = new Rectangle(0, 0, width, height);
            LinearGradientBrush brush = new LinearGradientBrush(rectangle, color1, color2, (float)angle);
            graphics.FillRectangle(brush, rectangle);
            graphics.Dispose();
            brush.Dispose();
            return bitmap;
        }
        #endregion
    }
}
