using System;
using SkiaDrawing;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines inference painter.
    /// </summary>
    [Serializable]
    public class Painter : IDisposable
    {
        #region Private data

        private Pen _boxPen = new Pen(Color.Red, 5);
        private Pen _pointPen = new Pen(Color.Yellow, 10);
        private Font _font = new Font("Arial", 24);

        #endregion

        #region Class components
        /// <summary>
        /// Gets or sets box pen.
        /// </summary>
        public Pen BoxPen
        {
            get
            {
                return _boxPen;
            }
            set
            {
                if ((_boxPen != null) && (!ReferenceEquals(_boxPen, value)))
                {
                    _boxPen.Dispose();
                }
                _boxPen = value;
            }
        }
        /// <summary>
        /// Gets or sets point pen.
        /// </summary>
        public Pen PointPen
        {
            get
            {
                return _pointPen;
            }
            set
            {
                if ((_pointPen != null) && (!ReferenceEquals(_pointPen, value)))
                {
                    _pointPen.Dispose();
                }
                _pointPen = value;
            }
        }
        /// <summary>
        /// Gets or sets text font.
        /// </summary>
        public Font TextFont
        {
            get
            {
                return _font;
            }
            set
            {
                if ((_font != null) && (!ReferenceEquals(_font, value)))
                {
                    _font.Dispose();
                }
                _font = value;
            }
        }
        /// <summary>
        /// Gets or sets text color.
        /// </summary>
        public Color TextColor { get; set; } = Color.White;
        /// <summary>
        /// Gets or sets box transparency.
        /// </summary>
        public byte Transparency { get; set; } = 0;
        /// <summary>
        /// Draw labels inside the box or not.
        /// </summary>
        public bool InsideBox { get; set; }
        /// <summary>
        /// Initializes inference painter.
        /// </summary>
        public Painter() { }
        #endregion

        #region Public draw methods
        /// <summary>
        /// Draws tracking and recognition results.
        /// </summary>
        /// <param name="graphics">Graphics</param>
        /// <param name="paintData">Paint data</param>
        public void Draw(Graphics graphics, params PaintData[] paintData)
        {
            int length = paintData.Length;

            for (int i = 0; i < length; i++)
            {
                if (!paintData[i].Rectangle.IsEmpty)
                {
                    Draw(
                        graphics,
                        paintData[i].Title,
                        paintData[i].Rectangle);
                }
            }

            for (int i = 0; i < length; i++)
            {
                if (paintData[i].Points?.Length > 0)
                {
                    Draw(
                        graphics,
                        paintData[i].Points);
                }
            }

            for (int i = 0; i < length; i++)
            {
                if (paintData[i].Labels?.Length > 0)
                    Draw(
                        graphics,
                        new Rectangle[] { paintData[i].Rectangle },
                        new string[][] { paintData[i].Labels });
            }

            return;
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Draws tracking and recognition results.
        /// </summary>
        /// <param name="graphics">Graphics</param>
        /// <param name="title">Title</param>
        /// <param name="rectangles">Rectangles</param>
        private void Draw(Graphics graphics, string title, params Rectangle[] rectangles)
        {
            using var textBrush = new SolidBrush(TextColor);
            using var inferenceBrush = new SolidBrush(BoxPen.Color);
            using var inferenceTransparentBrush = new SolidBrush(Color.FromArgb(Transparency, BoxPen.Color));
            var length = rectangles.Length;
            var depth = BoxPen.Width;

            for (int i = 0; i < length; i++)
            {
                var rectangle = rectangles[i];

                if (rectangle.IsEmpty)
                {
                    continue;
                }
                else
                {
                    graphics.FillRectangle(inferenceTransparentBrush, rectangle);
                    graphics.DrawRectangle(BoxPen, rectangle);

                    if (!string.IsNullOrEmpty(title))
                    {
                        var faceLabel = GetLabel(graphics, TextFont, rectangle, depth, title);
                        var w = graphics.MeasureString(faceLabel, TextFont);
                        var t = new RectangleF(
                            rectangle.X - depth / 2,
                            rectangle.Y - w.Height - depth / 2,
                            w.Width + depth,
                            w.Height);

                        graphics.FillRectangle(inferenceBrush, t);
                        graphics.DrawString(faceLabel, TextFont, textBrush, t.X + depth / 2, t.Y + depth / 2);
                    }
                }
            }

            return;
        }
        /// <summary>
        /// Draws tracking and recognition results.
        /// </summary>
        /// <param name="graphics">Graphics</param>
        /// <param name="points">Points</param>
        private void Draw(Graphics graphics, params Point[] points)
        {
            using var b = new SolidBrush(PointPen.Color);
            var length = points.Length;
            var depth = PointPen.Width;

            for (int i = 0; i < length; i++)
            {
                var point = points[i];

                if (point.IsEmpty)
                {
                    continue;
                }
                else
                {
                    graphics.FillEllipse(b, point.X - depth, point.Y - depth, depth, depth);
                }
            }
        }
        /// <summary>
        /// Draws tracking and recognition results.
        /// </summary>
        /// <param name="graphics">GraphicsBitmap</param>
        /// <param name="rectangles">Rectangles</param>
        /// <param name="labels">Labels</param>
        private void Draw(Graphics graphics, Rectangle[] rectangles, params string[][] labels)
        {
            if (rectangles.Length != labels.Length)
            {
                return;
            }
            else
            {
                using var textBrush = new SolidBrush(TextColor);
                using var inferenceBrush = new SolidBrush(BoxPen.Color);
                using var inferenceTransparentBrush = new SolidBrush(Color.FromArgb(Transparency, BoxPen.Color));
                var length = rectangles.Length;
                var depth = BoxPen.Width;

                for (int i = 0; i < length; i++)
                {
                    var rectangle = rectangles[i];

                    if (rectangle.IsEmpty)
                    {
                        continue;
                    }
                    else
                    {
                        var label = GetLabel(graphics, TextFont, rectangle, depth, labels[i]);

                        if (!string.IsNullOrEmpty(label))
                        {
                            var s = graphics.MeasureString(label, TextFont);
                            RectangleF r;

                            if (InsideBox)
                            {
                                r = new RectangleF(
                                    rectangle.X + depth / 2,
                                    rectangle.Y + depth / 2,
                                    rectangle.Width - depth,
                                    s.Height + depth);
                            }
                            else
                            {
                                r = new RectangleF(
                                    rectangle.X - depth / 2,
                                    rectangle.Y + rectangle.Height + depth / 2,
                                    rectangle.Width + depth,
                                    s.Height + depth);
                            }

                            graphics.FillRectangle(inferenceBrush, r);
                            graphics.DrawString(label, TextFont, textBrush, r.X + depth / 2, r.Y + depth / 2);
                        }
                    }
                }
            }
        }
        /// <summary>
        /// Returns label string.
        /// </summary>
        /// <param name="g">Graphics</param>
        /// <param name="font">Font</param>
        /// <param name="rectangle">Rectangle</param>
        /// <param name="depth">Depth</param>
        /// <param name="unit">Unit string</param>
        /// <returns>String</returns>
        private static string GetLabel(Graphics g, Font font, Rectangle rectangle, float depth, params string[] unit)
        {
            var label = string.Empty;
            var count = unit.Length;

            for (int j = 0; j < count; j++)
            {
                var line = string.Empty;
                var subline = unit[j];
                var length = subline.Length;

                for (int k = 0; k < length; k++)
                {
                    line += subline[k];

                    if (g.MeasureString(line + "..", font).Width >= rectangle.Width)
                    {
                        label += $"..{Environment.NewLine}";
                        line = string.Empty;
                    }

                    label += subline[k];
                }

                label += (j < count - 1) ? Environment.NewLine : string.Empty;
            }

            return label;
        }
        #endregion

        #region IDisposable

        private bool _disposed;

        /// <inheritdoc/>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <inheritdoc/>
        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                if (disposing)
                {
                    _boxPen?.Dispose();
                    _pointPen?.Dispose();
                    _font?.Dispose();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~Painter()
        {
            Dispose(false);
        }

        #endregion
    }
}
