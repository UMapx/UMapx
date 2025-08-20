using System;
using System.Collections.Generic;
using System.Drawing;
using UMapx.Core;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the figure to plotting in a Cartesian coordinate system.
    /// </summary>
    [Serializable]
    public class Figure : IDisposable
    {
        #region Private data
        private int figure_width, figure_height;
        private int canvas_width, canvas_height;
        private float xmin = -5, xmax = 5, ymin = -5, ymax = 5;
        private int xscale = 10, yscale = 10;
        private Style _style = new Style();
        private readonly List<GraphPane> _panes = new List<GraphPane>();
        private Bitmap _imagePane;
        #endregion

        #region Control voids and overrides
        /// <summary>
        /// Initializes the figure.
        /// </summary> 
        public Figure() { }
        /// <summary>
        /// Initializes the figure.
        /// </summary>
        /// <param name="style">Figure style</param>
        public Figure(Style style)
        {
            Style = style;
        }
        #endregion

        #region Figure properties
        /// <summary>
        /// Gets or sets figure style.
        /// </summary>
        public Style Style
        {
            get
            { 
                return _style; 
            }
            set
            {
                _style?.Dispose();
                _style = value;
            }
        }
        /// <summary>
        /// Gets or sets X label.
        /// </summary>
        public string LabelX { get; set; } = "Label X";
        /// <summary>
        /// Gets or sets Y label.
        /// </summary>
        public string LabelY { get; set; } = "Label Y";
        /// <summary>
        /// Gets or sets X title.
        /// </summary>
        public string Title { get; set; } = "Title";
        /// <summary>
        /// Gets or sets the legend.
        /// </summary>
        public Legend Legend { get; set; } = new Legend();
        /// <summary>
        /// Gets or sets X range [min, max].
        /// </summary>
        public RangeFloat RangeX
        {
            get
            {
                return new RangeFloat(xmin, xmax);
            }
            set
            {
                if (value.Min == value.Max)
                    throw new ArgumentOutOfRangeException("Start and end points cannot be the same");

                xmin = value.Min; xmax = value.Max;
            }
        }
        /// <summary>
        /// Gets or sets Y range [min, max].
        /// </summary>
        public RangeFloat RangeY
        {
            get
            {
                return new RangeFloat(ymin, ymax);
            }
            set
            {
                if (value.Min == value.Max)
                    throw new ArgumentOutOfRangeException("Start and end points cannot be the same");

                ymin = value.Min; ymax = value.Max;
            }
        }
        /// <summary>
        /// Gets or sets the scale of a range of digital elevations along the X and Y axes.
        /// </summary>
        public PointInt Marks
        {
            get
            {
                return new PointInt(xscale, yscale);
            }
            set
            {
                if (value.X < 1 || value.Y < 1)
                    throw new ArgumentException("The range of marks cannot be less than 1");

                xscale = value.X; yscale = value.Y;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor for the canvas [0.5, 0.8]. 
        /// </summary>
        public float Scaling { get; set; } = 0.7f;
        /// <summary>
        /// Gets or sets shapes.
        /// </summary>
        public bool Shapes { get; set; } = true;
        /// <summary>
        /// Gets or sets grid.
        /// </summary>
        public bool Grid { get; set; } = false;
        /// <summary>
        /// Gets or sets property of auto range axes.
        /// </summary>
        public bool AutoRange { get; set; } = true;
        #endregion

        #region Figure methods
        /// <summary>
        /// Draws figure to graphics object.
        /// </summary>
        /// <param name="graphics">Graphics</param>
        public void Draw(Graphics graphics)
        {
            #region Figure data
            // figure and canvas options
            var size = graphics.VisibleClipBounds.Size;
            figure_width  = (int)size.Width;
            figure_height = (int)size.Height;
            canvas_width  = (int)(figure_width * Scaling);
            canvas_height = (int)(figure_height * Scaling);

            // encapsulation of figure and canvas graphics
            using var figure = new Bitmap(figure_width, figure_height);
            using Graphics figure_graphics = Graphics.FromImage(figure);
            figure_graphics.Clear(Style.ColorFrame);

            using var canvas = new Bitmap(canvas_width, canvas_height);
            using Graphics canvas_graphics = Graphics.FromImage(canvas);
            canvas_graphics.Clear(Style.ColorBack);

            // offsets
            int dx = canvas_width / xscale;
            int dy = canvas_height / yscale;
            int dw = (figure_width - canvas_width) / 2;
            int dh = (figure_height - canvas_height) / 2;
            #endregion

            #region Autorange
            //autorange or not?
            if (_imagePane is object)
            {
                xmin = ymin = 0;
                xmax = size.Width;
                ymax = size.Height;
            }
            else if (AutoRange)
            {
                xmin = ymin = float.MaxValue;
                xmax = ymax = float.MinValue;

                foreach (var pane in _panes)
                {
                    xmin = Math.Min(xmin, pane.X.Min());
                    xmax = Math.Max(xmax, pane.X.Max());
                    ymin = Math.Min(ymin, pane.Y.Min());
                    ymax = Math.Max(ymax, pane.Y.Max());
                }
            }
            else
            {
                // user defined rangeX and rangeY
            }
            #endregion

            #region Numeric marks
            // points:
            float[] X = Points.GetPoints(xmin, xmax, xscale);
            float[] Y = Points.GetPoints(ymin, ymax, yscale);

            using SolidBrush br = new SolidBrush(Style.ColorMarks);
            using Pen pen1 = new Pen(Style.ColorGrid, Style.DepthShapes);
            using Pen pen2 = new Pen(Style.ColorShapes, Style.DepthShapes);
            int min = Math.Min(dx, dy), s = min / 8;
            string numerics;
            int xlength, ylength, i;
            int xpoint, ypoint;
            SizeF numsize;

            xlength = X.Length;
            float numX;

            for (i = 0; i < xlength; i++)
            {
                numX = X[i];
                numerics = GetNumString(numX);
                xpoint = (int)Points.X2Point(numX, xmin, xmax, canvas_width);
                figure_graphics.DrawString(numerics, Style.FontMarks, br, xpoint + dw - 5, canvas_height + dh + 5);
            }

            ylength = Y.Length;
            float numY;

            for (i = 0; i < ylength; i++)
            {
                numY = Y[i];
                numerics = GetNumString(numY);
                numsize = figure_graphics.MeasureString(numerics, Style.FontMarks);
                ypoint = (int)Points.Y2Point(numY, ymin, ymax, canvas_height);
                figure_graphics.DrawString(numerics, Style.FontMarks, br, dw - numsize.Width - 5, ypoint + dh - 10);
            }
            #endregion

            #region Grid painting
            if (Grid == true)
            {
                xlength = X.Length;

                for (i = 0; i < xlength; i++)
                {
                    numX = X[i];
                    xpoint = (int)Points.X2Point(numX, xmin, xmax, canvas_width);
                    canvas_graphics.DrawLine(pen1, xpoint, 0, xpoint, canvas_height);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, ymin, ymax, canvas_height);
                    canvas_graphics.DrawLine(pen1, 0, ypoint, canvas_width, ypoint);
                }
            }
            #endregion

            #region Title and labels
            Paint_Title(figure_graphics, Title, figure_width, figure_height, dw, dh);
            Paint_LabelX(figure_graphics, LabelX, figure_width, figure_height, dw, dh);
            Paint_LabelY(figure_graphics, LabelY, figure_width, figure_height, dw, dh);
            #endregion

            #region Graph painting
            if (_imagePane is object)
            {
                // 2D plotting
                var rectangle = new Rectangle(0, 0, canvas_width, canvas_height);
                canvas_graphics.DrawImage(_imagePane, rectangle);
            }
            //else
            {
                // 1D plotting
                float r0 = 4;
                float radius;

                // Plotting:
                foreach (var current in _panes)
                {
                    radius = (current.Depth + r0) * 2;

                    if (current.Pane == Pane.Plot)
                    {
                        switch (current.Symbol)
                        {
                            case Symbol.None:
                                PlotLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case Symbol.Circle:
                                PlotCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Ball:
                                PlotCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case Symbol.Rectangle:
                                PlotRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Polygon:
                                PlotRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                    else if (current.Pane == Pane.Stem)
                    {
                        switch (current.Symbol)
                        {
                            case Symbol.None:
                                StemLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case Symbol.Circle:
                                StemCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Ball:
                                StemCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case Symbol.Rectangle:
                                StemRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Polygon:
                                StemRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                    else
                    {
                        switch (current.Symbol)
                        {
                            case Symbol.None:
                                ScatterLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case Symbol.Circle:
                                ScatterCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Ball:
                                ScatterCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case Symbol.Rectangle:
                                ScatterRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case Symbol.Polygon:
                                ScatterRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                }
            }
            #endregion

            #region Shapes
            if (Shapes == true)
            {
                xlength = X.Length;

                for (i = 0; i < xlength; i++)
                {
                    numX = X[i];
                    xpoint = (int)Points.X2Point(numX, xmin, xmax, canvas_width);
                    canvas_graphics.DrawLine(pen2, xpoint, canvas_height, xpoint, canvas_height - s);
                    canvas_graphics.DrawLine(pen2, xpoint, 0, xpoint, s);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, ymin, ymax, canvas_height);
                    canvas_graphics.DrawLine(pen2, 0, ypoint, s, ypoint);
                    canvas_graphics.DrawLine(pen2, canvas_width, ypoint, canvas_width - s, ypoint);
                }
            }

            canvas_graphics.DrawRectangle(pen2, 0, 0, canvas_width - 1, canvas_height - 1);
            #endregion

            #region Legend
            if (Legend != null && Legend.Show)
            {
                Paint_Legend(canvas_graphics);
            }
            #endregion

            #region Merging
            graphics.DrawImage(figure, new Point(0, 0));
            graphics.DrawImage(canvas, new Point(dw, dh));
            #endregion
        }
        /// <summary>
        /// Show image at the figure.
        /// </summary>
        /// <param name="bitmap">Bitmap</param>
        public void Image(Bitmap bitmap)
        {
            _imagePane?.Dispose();
            _imagePane = (Bitmap)bitmap.Clone();
        }
        /// <summary>
        /// Adds a graph pane.
        /// </summary>
        /// <param name="pane">Graph pane</param>
        public void Graph(GraphPane pane)
        {
            if (pane.X.Length != pane.Y.Length)
                throw new ArgumentException("Vectors must be of the same length");

            _panes.Add(pane);
        }
        /// <summary>
        /// Removes all graphs from the figure. 
        /// </summary>
        public void Clear()
        {
            _panes.Clear();
            _imagePane?.Dispose();
        }
        #endregion

        #region Private voids

        #region Legend voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="g"></param>
        /// <param name="cx"></param>
        /// <param name="cy"></param>
        /// <param name="size"></param>
        /// <param name="color"></param>
        /// <param name="type"></param>
        /// <param name="depth"></param>
        private void DrawLegendMarker(Graphics g, int cx, int cy, int size, Color color, Symbol type, float depth)
        {
            int half = size / 2;
            var rect = new Rectangle(cx, cy - half, size, size);

            using var pen = new Pen(color, Math.Max(1f, depth));
            using var br = new SolidBrush(color);

            switch (type)
            {
                case Symbol.None:
                    g.DrawLine(pen, cx, cy, cx + size, cy);
                    var r2 = new Rectangle(cx + size / 2 - 2, cy - 2, 4, 4);
                    g.FillRectangle(br, r2);
                    break;

                case Symbol.Circle:
                case Symbol.Ball:
                    if (IsFilled(type)) g.FillEllipse(br, rect);
                    else g.DrawEllipse(pen, rect);
                    break;

                case Symbol.Rectangle:
                case Symbol.Polygon:
                    if (IsFilled(type)) g.FillRectangle(br, rect);
                    else g.DrawRectangle(pen, rect);
                    break;

                default:
                    g.DrawLine(pen, cx, cy, cx + size, cy);
                    break;
            }

            bool IsFilled(Symbol type)
            {
                return type == Symbol.Ball || type == Symbol.Polygon;
            }
        }
        #endregion

        #region Plot function voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        private void PlotLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            int i, length = y.Length;
            Point[] points = new Point[x.Length];
            using Pen pen = new Pen(color, depth);
            float xi, yi;

            for (i = 0; i < length; i++)
            {
                xi = Points.ClipPoint(x[i], xmin, xmax);
                yi = Points.ClipPoint(y[i], ymin, ymax);

                if (Points.IsSingularPoint(xi)) continue;
                if (Points.IsSingularPoint(yi)) continue;

                points[i].X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                points[i].Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);
            }
            graphics.DrawLines(pen, points);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void PlotCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            int i, length = y.Length;
            Point[] points = new Point[x.Length];
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);
                    points[i] = point;

                    graphics.FillEllipse(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);
                    points[i] = point;

                    graphics.DrawEllipse(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            graphics.DrawLines(pen, points);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void PlotRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            int i, length = y.Length;
            Point[] points = new Point[x.Length];
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);
                    points[i] = point;

                    graphics.FillRectangle(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);
                    points[i] = point;

                    graphics.DrawRectangle(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            graphics.DrawLines(pen, points);
        }
        #endregion

        #region Stem function voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        private void StemLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            Point point = new Point();
            using Pen pen = new Pen(color, depth);
            int i, length = y.Length;
            float dy = (ymax - ymin) / canvas_height;
            int zero = (int)(canvas_height + ymin / dy);
            float xi, yi;

            for (i = 0; i < length; i++)
            {
                xi = Points.ClipPoint(x[i], xmin, xmax);
                yi = Points.ClipPoint(y[i], ymin, ymax);

                if (Points.IsSingularPoint(xi)) continue;
                if (Points.IsSingularPoint(yi)) continue;

                point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                graphics.DrawLine(pen, point.X, zero, point.X, point.Y);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void StemCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            int i, length = y.Length;
            float dy = (ymax - ymin) / canvas_height;
            int zero = (int)(canvas_height + ymin / dy);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.FillEllipse(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                    graphics.DrawLine(pen, point.X, zero, point.X, point.Y);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.DrawEllipse(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                    graphics.DrawLine(pen, point.X, zero, point.X, point.Y);
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void StemRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill = false)
        {
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            int i, length = y.Length;
            float dy = (ymax - ymin) / canvas_height;
            int zero = (int)(canvas_height + ymin / dy);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.FillRectangle(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                    graphics.DrawLine(pen, point.X, zero, point.X, point.Y);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.DrawRectangle(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                    graphics.DrawLine(pen, point.X, zero, point.X, point.Y);
                }
            }
        }
        #endregion

        #region Scatter function voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        private void ScatterLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void ScatterCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            int i, length = y.Length;
            float dy = (ymax - ymin) / canvas_height;
            int zero = (int)(canvas_height + ymin / dy);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.FillEllipse(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.DrawEllipse(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="depth"></param>
        /// <param name="color"></param>
        /// <param name="radius"></param>
        /// <param name="fill"></param>
        private void ScatterRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill = false)
        {
            Point point = new Point();
            using SolidBrush br = new SolidBrush(color);
            using Pen pen = new Pen(color, depth);
            int i, length = y.Length;
            float dy = (ymax - ymin) / canvas_height;
            int zero = (int)(canvas_height + ymin / dy);
            float xi, yi;

            if (fill == true)
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.FillRectangle(br, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    xi = Points.ClipPoint(x[i], xmin, xmax);
                    yi = Points.ClipPoint(y[i], ymin, ymax);

                    if (Points.IsSingularPoint(xi)) continue;
                    if (Points.IsSingularPoint(yi)) continue;

                    point.X = (int)Points.X2Point(xi, xmin, xmax, canvas_width);
                    point.Y = (int)Points.Y2Point(yi, ymin, ymax, canvas_height);

                    graphics.DrawRectangle(pen, point.X - radius / 2, point.Y - radius / 2, radius, radius);
                }
            }
        }
        #endregion

        #region Painter voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="g"></param>
        private void Paint_Legend(Graphics g)
        {
            if (_panes.Count == 0) return;

            using var textBrush = new SolidBrush(Style.ColorText);
            var font = Style.FontMarks;
            int marker = Legend.MarkerSize;
            int rowH = Math.Max(Legend.RowHeight, marker);
            int maxTextW = 0;

            foreach (var it in _panes)
            {
                var sz = g.MeasureString(it.Label, font);
                if (sz.Width > maxTextW) maxTextW = (int)Math.Ceiling(sz.Width);
            }

            int innerPad = 8;
            int contentW = marker + Legend.MarkerGap + maxTextW;
            int contentH = _panes.Count * rowH;
            int boxW = contentW + innerPad * 2;
            int boxH = contentH + innerPad * 2;

            int x, y;

            switch (Legend.Anchor)
            {
                default:
                case LegendAnchor.TopRight:
                    x = canvas_width - Legend.Padding - boxW;
                    y = Legend.Padding;
                    break;
                case LegendAnchor.TopLeft:
                    x = Legend.Padding;
                    y = Legend.Padding;
                    break;
                case LegendAnchor.BottomRight:
                    x = canvas_width - Legend.Padding - boxW;
                    y = canvas_height - Legend.Padding - boxH;
                    break;
                case LegendAnchor.BottomLeft:
                    x = Legend.Padding;
                    y = canvas_height - Legend.Padding - boxH;
                    break;
            }

            var back = Color.FromArgb((int)(Legend.Opacity * 255), Style.ColorBack);
            using var backBrush = new SolidBrush(back);
            using var borderPen = new Pen(Style.ColorShapes, 1f);

            g.FillRectangle(backBrush, x, y, boxW, boxH);

            if (Legend.Border)
                g.DrawRectangle(borderPen, x, y, boxW, boxH);

            int cx = x + innerPad;
            int cy = y + innerPad;

            for (int idx = 0; idx < _panes.Count; idx++)
            {
                var it = _panes[idx];
                int my = cy + idx * rowH + rowH / 2;

                DrawLegendMarker(g, cx, my, marker, it.Color, it.Symbol, it.Depth);
                g.DrawString(it.Label, font, textBrush, cx + marker + Legend.MarkerGap, my - rowH / 2 + 1);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="title"></param>
        /// <param name="sizeX"></param>
        /// <param name="sizeY"></param>
        /// <param name="dw"></param>
        /// <param name="dh"></param>
        private void Paint_Title(Graphics graphics, string title, int sizeX, int sizeY, int dw, int dh)
        {
            using var format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            using var br = new SolidBrush(Style.ColorText);
            graphics.DrawString(title, Style.FontText, br, new PointF(sizeX / 2, dh / 2), format);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="xlabel"></param>
        /// <param name="sizeX"></param>
        /// <param name="sizeY"></param>
        /// <param name="dw"></param>
        /// <param name="dh"></param>
        private void Paint_LabelX(Graphics graphics, string xlabel, int sizeX, int sizeY, int dw, int dh)
        {
            using var format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            var size = graphics.MeasureString(xlabel, Style.FontText);
            using var br = new SolidBrush(Style.ColorText);
            graphics.DrawString(xlabel, Style.FontText, br, new PointF(sizeX / 2, sizeY - dh / 2 + size.Height / 4), format);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="graphics"></param>
        /// <param name="ylabel"></param>
        /// <param name="sizeX"></param>
        /// <param name="sizeY"></param>
        /// <param name="dw"></param>
        /// <param name="dh"></param>
        private void Paint_LabelY(Graphics graphics, string ylabel, int sizeX, int sizeY, int dw, int dh)
        {
            using var fmt = new StringFormat
            {
                Alignment = StringAlignment.Center,
                LineAlignment = StringAlignment.Center
            };

            var state = graphics.Save();

            float cx = dw / 2f;
            float cy = sizeY / 2f;

            using var br = new SolidBrush(Style.ColorText);
            graphics.TranslateTransform(cx, cy);
            graphics.RotateTransform(-90f);
            graphics.DrawString(ylabel, Style.FontText, br, PointF.Empty, fmt);
            graphics.Restore(state);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="num"></param>
        /// <returns></returns>
        private string GetNumString(float num)
        {
            string numerics = num.ToString("0.00");
            int l1 = numerics.Length - 1;

            if (numerics[l1] == '0')
            {
                numerics = numerics.Substring(0, l1);
                int l2 = l1 - 1;

                if (numerics[l2] == '0')
                {
                    numerics = numerics.Substring(0, l2 - 1);
                }
                return numerics;
            }

            return numerics;
        }
        #endregion

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
                    _panes.Clear();
                    _imagePane?.Dispose();
                    _style?.Dispose();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~Figure()
        {
            Dispose(false);
        }

        #endregion
    }
}
