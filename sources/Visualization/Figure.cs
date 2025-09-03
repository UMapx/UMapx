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
    public class Figure
    {
        #region Private data
        private int _figure_width, _figure_height;
        private int _canvas_width, _canvas_height;
        private int _xscale = 10, _yscale = 10;
        private float _xmin = -5, _xmax = 5;
        private float _ymin = -5, _ymax = 5;
        private float _scaling = 0.65f;
        private readonly FigureStyle _style;
        private readonly List<GraphPane> _panes = new List<GraphPane>();
        private Bitmap _imagePane;
        #endregion

        #region Figure constructor
        /// <summary>
        /// Initializes the figure.
        /// </summary>
        /// <param name="style">Figure style</param>
        public Figure(FigureStyle style)
        {
            _style = style;
        }
        #endregion

        #region Figure properties
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
        /// Gets or sets the grid.
        /// </summary>
        public Grid Grid { get; set; } = new Grid();
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
                return new RangeFloat(_xmin, _xmax);
            }
            set
            {
                if (value.Min == value.Max)
                    throw new ArgumentOutOfRangeException("Start and end points cannot be the same");

                if (MathF.IsSingular(value.Min) || MathF.IsSingular(value.Max))
                    throw new ArgumentOutOfRangeException("Start of end points cannot be singular");

                _xmin = value.Min; _xmax = value.Max;
            }
        }
        /// <summary>
        /// Gets or sets Y range [min, max].
        /// </summary>
        public RangeFloat RangeY
        {
            get
            {
                return new RangeFloat(_ymin, _ymax);
            }
            set
            {
                if (value.Min == value.Max)
                    throw new ArgumentOutOfRangeException("Start and end points cannot be the same");

                if (MathF.IsSingular(value.Min) || MathF.IsSingular(value.Max))
                    throw new ArgumentOutOfRangeException("Start of end points cannot be singular");

                _ymin = value.Min; _ymax = value.Max;
            }
        }
        /// <summary>
        /// Gets or sets the scale of a range of digital elevations along the X and Y axes.
        /// </summary>
        public PointInt Marks
        {
            get
            {
                return new PointInt(_xscale, _yscale);
            }
            set
            {
                if (value.X < 1 || value.Y < 1)
                    throw new ArgumentOutOfRangeException("The range of marks cannot be less than 1");

                _xscale = value.X; _yscale = value.Y;
            }
        }
        /// <summary>
        /// Gets or sets the scale factor for the canvas [0.5, 0.8]. 
        /// </summary>
        public float Scaling
        {
            get 
            { 
                return _scaling; 
            }
            set 
            {
                if (value <= 0 || value > 1)
                    throw new ArgumentOutOfRangeException("Scale factor cannot be less than 0 or more than 1");

                _scaling = value; 
            }
        }
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
            _figure_width  = (int)size.Width;
            _figure_height = (int)size.Height;
            _canvas_width  = (int)(_figure_width * _scaling);
            _canvas_height = (int)(_figure_height * _scaling);

            // encapsulation of figure and canvas graphics
            using var figure = new Bitmap(_figure_width, _figure_height);
            using Graphics figure_graphics = Graphics.FromImage(figure);
            figure_graphics.Clear(_style.ColorFrame);
            figure_graphics.CompositingQuality = System.Drawing.Drawing2D.CompositingQuality.HighQuality;
            figure_graphics.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBilinear;
            figure_graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            using var canvas = new Bitmap(_canvas_width, _canvas_height);
            using Graphics canvas_graphics = Graphics.FromImage(canvas);
            canvas_graphics.Clear(_style.ColorBack);
            canvas_graphics.CompositingQuality = System.Drawing.Drawing2D.CompositingQuality.HighQuality;
            canvas_graphics.InterpolationMode = System.Drawing.Drawing2D.InterpolationMode.HighQualityBilinear;
            canvas_graphics.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.HighQuality;

            // offsets
            int dx = _canvas_width / _xscale;
            int dy = _canvas_height / _yscale;
            int dw = (_figure_width - _canvas_width) / 2;
            int dh = (_figure_height - _canvas_height) / 2;
            #endregion

            #region Autorange
            //autorange or not?
            if (_imagePane is object)
            {
                _xmin = _ymin = 0;
                _xmax = _imagePane.Width;
                _ymax = _imagePane.Height;
            }
            else if (AutoRange && _panes.Count != 0)
            {
                var xmin = float.PositiveInfinity; var xmax = float.NegativeInfinity;
                var ymin = float.PositiveInfinity; var ymax = float.NegativeInfinity;

                foreach (var pane in _panes)
                {
                    xmin = Math.Min(xmin, pane.X.GetMin() ?? xmin);
                    xmax = Math.Max(xmax, pane.X.GetMax() ?? xmax);
                    ymin = Math.Min(ymin, pane.Y.GetMin() ?? ymin);
                    ymax = Math.Max(ymax, pane.Y.GetMax() ?? ymax);
                }

                _xmin = MathF.IsSingular(xmin) ? _xmin : xmin;
                _xmax = MathF.IsSingular(xmax) ? _xmax : xmax;
                _ymin = MathF.IsSingular(ymin) ? _ymin : ymin;
                _ymax = MathF.IsSingular(ymax) ? _ymax : ymax;
            }
            else
            {
                // user defined rangeX and rangeY
            }
            #endregion

            #region Numeric marks
            // points:
            float[] X = Points.GetPoints(_xmin, _xmax, _xscale);
            float[] Y = Points.GetPoints(_ymin, _ymax, _yscale);

            using SolidBrush br = new SolidBrush(_style.ColorMarks);
            using Pen pen1 = new Pen(_style.ColorGrid, _style.DepthShapes);
            using Pen pen2 = new Pen(_style.ColorShapes, _style.DepthShapes);
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
                xpoint = (int)Points.X2Point(numX, _xmin, _xmax, _canvas_width);
                figure_graphics.DrawString(numerics, _style.FontMarks, br, xpoint + dw - 5, _canvas_height + dh + 5);
            }

            ylength = Y.Length;
            float numY;

            for (i = 0; i < ylength; i++)
            {
                numY = Y[i];
                numerics = GetNumString(numY);
                numsize = figure_graphics.MeasureString(numerics, _style.FontMarks);
                ypoint = (int)Points.Y2Point(numY, _ymin, _ymax, _canvas_height);
                figure_graphics.DrawString(numerics, _style.FontMarks, br, dw - numsize.Width - 5, ypoint + dh - 10);
            }
            #endregion

            #region Grid painting
            if (Grid.Show)
            {
                xlength = X.Length;

                for (i = 0; i < xlength; i++)
                {
                    numX = X[i];
                    xpoint = (int)Points.X2Point(numX, _xmin, _xmax, _canvas_width);
                    DrawGridVertical(canvas_graphics, xpoint, _canvas_height, pen1);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, _ymin, _ymax, _canvas_height);
                    DrawGridHorizontal(canvas_graphics, ypoint, _canvas_width, pen1);
                }
            }
            #endregion

            #region Title and labels
            Paint_Title(figure_graphics, Title, _figure_width, _figure_height, dw, dh);
            Paint_LabelX(figure_graphics, LabelX, _figure_width, _figure_height, dw, dh);
            Paint_LabelY(figure_graphics, LabelY, _figure_width, _figure_height, dw, dh);
            #endregion

            #region Graph painting
            if (_imagePane is object)
            {
                // 2D plotting
                var rectangle = new Rectangle(0, 0, _canvas_width, _canvas_height);
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

                    if (current.PaneType == PaneType.Plot)
                    {
                        switch (current.GraphType)
                        {
                            case GraphType.None:
                                PlotLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case GraphType.Circle:
                                PlotCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Ball:
                                PlotCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case GraphType.Rectangle:
                                PlotRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Polygon:
                                PlotRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                    else if (current.PaneType == PaneType.Stem)
                    {
                        switch (current.GraphType)
                        {
                            case GraphType.None:
                                StemLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case GraphType.Circle:
                                StemCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Ball:
                                StemCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case GraphType.Rectangle:
                                StemRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Polygon:
                                StemRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                    else
                    {
                        switch (current.GraphType)
                        {
                            case GraphType.None:
                                ScatterLine(canvas_graphics, current.X, current.Y, current.Depth, current.Color);
                                break;

                            case GraphType.Circle:
                                ScatterCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Ball:
                                ScatterCircle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;

                            case GraphType.Rectangle:
                                ScatterRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, false);
                                break;

                            case GraphType.Polygon:
                                ScatterRectangle(canvas_graphics, current.X, current.Y, current.Depth, current.Color, radius, true);
                                break;
                        }
                    }
                }
            }
            #endregion

            #region Shapes
            if (Grid.Shapes)
            {
                xlength = X.Length;

                for (i = 0; i < xlength; i++)
                {
                    numX = X[i];
                    xpoint = (int)Points.X2Point(numX, _xmin, _xmax, _canvas_width);
                    canvas_graphics.DrawLine(pen2, xpoint, _canvas_height, xpoint, _canvas_height - s);
                    canvas_graphics.DrawLine(pen2, xpoint, 0, xpoint, s);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, _ymin, _ymax, _canvas_height);
                    canvas_graphics.DrawLine(pen2, 0, ypoint, s, ypoint);
                    canvas_graphics.DrawLine(pen2, _canvas_width, ypoint, _canvas_width - s, ypoint);
                }
            }

            canvas_graphics.DrawRectangle(pen2, 0, 0, _canvas_width - 1, _canvas_height - 1);
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
            _imagePane = bitmap;
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
            _imagePane = null;
            _xmin = -5; _xmax = 5;
            _ymin = -5; _ymax = 5;
        }
        #endregion

        #region Private voids

        #region Grid voids
        /// <summary>
        /// Draws a vertical grid line at the specified x pixel using the configured grid style.
        /// </summary>
        /// <remarks>
        /// - Uses <see cref="Grid.Style"/> to choose between solid, dashed (custom pattern), or dotted lines.<br/>
        /// - Clones <paramref name="basePen"/> when style customization is required to avoid mutating the caller's pen.
        /// </remarks>
        /// <param name="g">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X pixel coordinate where the grid line is drawn</param>
        /// <param name="height">Total canvas height in pixels</param>
        /// <param name="basePen">Base pen (color/width) to use; will be cloned for style-specific tweaks</param>
        private void DrawGridVertical(Graphics g, int x, int height, Pen basePen)
        {
            switch (Grid.Style)
            {
                case GridStyle.Solid:
                    g.DrawLine(basePen, x, 0, x, height);
                    break;

                case GridStyle.Dashed:
                    using (var pen = (Pen)basePen.Clone())
                    {
                        pen.DashStyle = System.Drawing.Drawing2D.DashStyle.Custom;
                        pen.DashPattern = new float[] { Math.Max(1f, Grid.DashLength), Math.Max(1f, Grid.GapLength) };
                        g.DrawLine(pen, x, 0, x, height);
                    }
                    break;

                case GridStyle.Dot:
                    using (var pen = (Pen)basePen.Clone())
                    {
                        pen.DashStyle = System.Drawing.Drawing2D.DashStyle.Dot;
                        pen.DashCap = System.Drawing.Drawing2D.DashCap.Round;
                        g.DrawLine(pen, x, 0, x, height);
                    }
                    break;
            }
        }
        /// <summary>
        /// Draws a horizontal grid line at the specified y pixel using the configured grid style.
        /// </summary>
        /// <remarks>
        /// - Uses <see cref="Grid.Style"/> to choose between solid, dashed (custom pattern), or dotted lines.<br/>
        /// - Clones <paramref name="basePen"/> when style customization is required to avoid mutating the caller's pen.
        /// </remarks>
        /// <param name="g">Target <see cref="Graphics"/> surface</param>
        /// <param name="y">Y pixel coordinate where the grid line is drawn</param>
        /// <param name="width">Total canvas width in pixels</param>
        /// <param name="basePen">Base pen (color/width) to use; will be cloned for style-specific tweaks</param>
        private void DrawGridHorizontal(Graphics g, int y, int width, Pen basePen)
        {
            switch (Grid.Style)
            {
                case GridStyle.Solid:
                    g.DrawLine(basePen, 0, y, width, y);
                    break;

                case GridStyle.Dashed:
                    using (var pen = (Pen)basePen.Clone())
                    {
                        pen.DashStyle = System.Drawing.Drawing2D.DashStyle.Custom;
                        pen.DashPattern = new float[] { Math.Max(1f, Grid.DashLength), Math.Max(1f, Grid.GapLength) };
                        g.DrawLine(pen, 0, y, width, y);
                    }
                    break;

                case GridStyle.Dot:
                    using (var pen = (Pen)basePen.Clone())
                    {
                        pen.DashStyle = System.Drawing.Drawing2D.DashStyle.Dot;
                        pen.DashCap = System.Drawing.Drawing2D.DashCap.Round;
                        g.DrawLine(pen, 0, y, width, y);
                    }
                    break;
            }
        }

        #endregion

        #region Legend voids
        /// <summary>
        /// Draws a legend marker (shape/line sample) at the specified position.
        /// </summary>
        /// <remarks>
        /// - Marker appearance depends on <paramref name="type"/>; filled variants use <paramref name="depth"/> as stroke width.<br/>
        /// - For <see cref="GraphType.None"/> a short line sample with a small square is drawn.
        /// </remarks>
        /// <param name="g">Target <see cref="Graphics"/> surface</param>
        /// <param name="cx">Left X pixel of the marker box</param>
        /// <param name="cy">Vertical center Y pixel of the marker box</param>
        /// <param name="size">Marker box size (width and height)</param>
        /// <param name="color">Marker color</param>
        /// <param name="type">Graph marker type</param>
        /// <param name="depth">Stroke thickness for outline</param>
        private void DrawLegendMarker(Graphics g, int cx, int cy, int size, Color color, GraphType type, float depth)
        {
            int half = size / 2;
            var rect = new Rectangle(cx, cy - half, size, size);

            using var pen = new Pen(color, Math.Max(1f, depth));
            using var br = new SolidBrush(color);

            switch (type)
            {
                case GraphType.None:
                    g.DrawLine(pen, cx, cy, cx + size, cy);
                    var r2 = new Rectangle(cx + size / 2 - 2, cy - 2, 4, 4);
                    g.FillRectangle(br, r2);
                    break;

                case GraphType.Circle:
                case GraphType.Ball:
                    if (IsFilled(type)) g.FillEllipse(br, rect);
                    else g.DrawEllipse(pen, rect);
                    break;

                case GraphType.Rectangle:
                case GraphType.Polygon:
                    if (IsFilled(type)) g.FillRectangle(br, rect);
                    else g.DrawRectangle(pen, rect);
                    break;

                default:
                    g.DrawLine(pen, cx, cy, cx + size, cy);
                    break;
            }

            bool IsFilled(GraphType type)
            {
                return type == GraphType.Ball || type == GraphType.Polygon;
            }
        }
        #endregion

        #region Plot function voids
        /// <summary>
        /// Plots a polyline through (x, y) using the given stroke.
        /// </summary>
        /// <remarks>
        /// - Invalid or out-of-range points (NaN/Inf/±∞ or clipped by axes) are skipped and break the polyline into segments.<br/>
        /// - Uses <see cref="DrawPolylineSkipInvalid"/> for robust rendering.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stroke thickness</param>
        /// <param name="color">Stroke color</param>
        private void PlotLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            using var pen = new Pen(color, depth);
            DrawPolylineSkipInvalid(graphics, pen, x, y);
        }
        /// <summary>
        /// Plots circular markers at (x, y) and connects valid points with a polyline.
        /// </summary>
        /// <remarks>
        /// - Each point is clipped to the current axes; invalid points are skipped.<br/>
        /// - When <paramref name="fill"/> is true, filled discs are drawn; otherwise only outlines are drawn.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stroke thickness for outlines and connecting line</param>
        /// <param name="color">Marker and line color</param>
        /// <param name="radius">Marker diameter in pixels</param>
        /// <param name="fill">Whether to fill the markers</param>
        private void PlotCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);

            int length = y.Length;
            for (int i = 0; i < length; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi)) continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillEllipse(br, px - radius / 2, py - radius / 2, radius, radius);
                else graphics.DrawEllipse(pen, px - radius / 2, py - radius / 2, radius, radius);
            }

            DrawPolylineSkipInvalid(graphics, pen, x, y);
        }
        /// <summary>
        /// Plots square markers at (x, y) and connects valid points with a polyline.
        /// </summary>
        /// <remarks>
        /// - Each point is clipped to the current axes; invalid points are skipped.<br/>
        /// - When <paramref name="fill"/> is true, filled squares are drawn; otherwise only outlines are drawn.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stroke thickness for outlines and connecting line</param>
        /// <param name="color">Marker and line color</param>
        /// <param name="radius">Marker side length in pixels</param>
        /// <param name="fill">Whether to fill the markers</param>
        private void PlotRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);

            int length = y.Length;
            for (int i = 0; i < length; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi)) continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillRectangle(br, px - radius / 2, py - radius / 2, radius, radius);
                else graphics.DrawRectangle(pen, px - radius / 2, py - radius / 2, radius, radius);
            }

            DrawPolylineSkipInvalid(graphics, pen, x, y);
        }
        /// <summary>
        /// Draws a polyline through (x, y) while skipping invalid or clipped points, splitting into segments.
        /// </summary>
        /// <remarks>
        /// - Converts world coordinates to device pixels via <c>Points.X2Point</c> and <c>Points.Y2Point</c>.<br/>
        /// - Accumulates a segment until an invalid point is encountered, then draws and starts a new segment.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="pen">Pen to draw with</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        private void DrawPolylineSkipInvalid(Graphics graphics, Pen pen, float[] x, float[] y)
        {
            var seg = new List<Point>(Math.Min(x.Length, y.Length));

            for (int i = 0; i < y.Length; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi))
                {
                    if (seg.Count > 1) graphics.DrawLines(pen, seg.ToArray());
                    seg.Clear();
                    continue;
                }

                seg.Add(new Point(
                    (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width),
                    (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height)));
            }

            if (seg.Count > 1) graphics.DrawLines(pen, seg.ToArray());
        }
        #endregion

        #region Stem function voids
        /// <summary>
        /// Renders a classic stem plot: vertical lines from y = 0 to each data point.
        /// </summary>
        /// <remarks>
        /// - Uses the current Y-axis transform to locate the zero baseline.<br/>
        /// - Invalid/clipped points are skipped.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stem thickness</param>
        /// <param name="color">Stem color</param>
        private void StemLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            using var pen = new Pen(color, depth);
            int zero = (int)Points.Y2Point(0f, _ymin, _ymax, _canvas_height);

            int n = y.Length;
            for (int i = 0; i < n; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi))
                    continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                graphics.DrawLine(pen, px, zero, px, py);
            }
        }
        /// <summary>
        /// Renders a stem plot with circular markers at the stem tips.
        /// </summary>
        /// <remarks>
        /// - Draws each marker (filled or outlined) and a vertical stem to y = 0.<br/>
        /// - Invalid/clipped points are skipped.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stem/outline thickness</param>
        /// <param name="color">Marker and stem color</param>
        /// <param name="radius">Marker diameter in pixels</param>
        /// <param name="fill">Whether to fill the marker</param>
        private void StemCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);
            int zero = (int)Points.Y2Point(0f, _ymin, _ymax, _canvas_height);
            float r2 = radius / 2f;

            int n = y.Length;
            for (int i = 0; i < n; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi))
                    continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillEllipse(br, px - r2, py - r2, radius, radius);
                else graphics.DrawEllipse(pen, px - r2, py - r2, radius, radius);

                graphics.DrawLine(pen, px, zero, px, py);
            }
        }
        /// <summary>
        /// Renders a stem plot with square markers at the stem tips.
        /// </summary>
        /// <remarks>
        /// - Draws each marker (filled or outlined) and a vertical stem to y = 0.<br/>
        /// - Invalid/clipped points are skipped.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stem/outline thickness</param>
        /// <param name="color">Marker and stem color</param>
        /// <param name="radius">Marker side length in pixels</param>
        /// <param name="fill">Whether to fill the marker</param>
        private void StemRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill = false)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);
            int zero = (int)Points.Y2Point(0f, _ymin, _ymax, _canvas_height);
            float r2 = radius / 2f;

            int n = y.Length;
            for (int i = 0; i < n; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi))
                    continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillRectangle(br, px - r2, py - r2, radius, radius);
                else graphics.DrawRectangle(pen, px - r2, py - r2, radius, radius);

                graphics.DrawLine(pen, px, zero, px, py);
            }
        }
        #endregion

        #region Scatter function voids
        /// <summary>
        /// Draws a polyline for scatter data that should be connected.
        /// </summary>
        /// <remarks>
        /// - Delegates to <see cref="DrawPolylineSkipInvalid"/> to handle invalid/clipped points.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Stroke thickness</param>
        /// <param name="color">Stroke color</param>
        private void ScatterLine(Graphics graphics, float[] x, float[] y, float depth, Color color)
        {
            using var pen = new Pen(color, depth);
            DrawPolylineSkipInvalid(graphics, pen, x, y);
        }
        /// <summary>
        /// Draws unconnected circular markers at (x, y).
        /// </summary>
        /// <remarks>
        /// - Invalid/clipped points are skipped. No connecting line is drawn.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Outline thickness when <paramref name="fill"/> is false</param>
        /// <param name="color">Marker color</param>
        /// <param name="radius">Marker diameter in pixels</param>
        /// <param name="fill">Whether to fill the markers</param>
        private void ScatterCircle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);
            float r2 = radius / 2f;

            int n = y.Length;
            for (int i = 0; i < n; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi)) continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillEllipse(br, px - r2, py - r2, radius, radius);
                else graphics.DrawEllipse(pen, px - r2, py - r2, radius, radius);
            }
        }
        /// <summary>
        /// Draws unconnected square markers at (x, y).
        /// </summary>
        /// <remarks>
        /// - Invalid/clipped points are skipped. No connecting line is drawn.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="x">X data in world coordinates</param>
        /// <param name="y">Y data in world coordinates</param>
        /// <param name="depth">Outline thickness when <paramref name="fill"/> is false</param>
        /// <param name="color">Marker color</param>
        /// <param name="radius">Marker side length in pixels</param>
        /// <param name="fill">Whether to fill the markers</param>
        private void ScatterRectangle(Graphics graphics, float[] x, float[] y, float depth, Color color, float radius, bool fill = false)
        {
            using var br = new SolidBrush(color);
            using var pen = new Pen(color, depth);
            float r2 = radius / 2f;

            int n = y.Length;
            for (int i = 0; i < n; i++)
            {
                float xi = Points.ClipPoint(x[i], _xmin, _xmax);
                float yi = Points.ClipPoint(y[i], _ymin, _ymax);

                if (Points.IsSingularPoint(xi) || Points.IsSingularPoint(yi)) continue;

                int px = (int)Points.X2Point(xi, _xmin, _xmax, _canvas_width);
                int py = (int)Points.Y2Point(yi, _ymin, _ymax, _canvas_height);

                if (fill) graphics.FillRectangle(br, px - r2, py - r2, radius, radius);
                else graphics.DrawRectangle(pen, px - r2, py - r2, radius, radius);
            }
        }
        #endregion

        #region Painter voids
        /// <summary>
        /// Paints the plot legend box with marker samples and labels for each pane/series.
        /// </summary>
        /// <remarks>
        /// - Computes content size from labels and marker size, positions the box by <see cref="Legend.Anchor"/>.<br/>
        /// - Applies background opacity and optional border based on legend style settings.
        /// </remarks>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        private void Paint_Legend(Graphics graphics)
        {
            if (_panes.Count == 0) return;

            using var textBrush = new SolidBrush(_style.ColorText);
            var font = _style.FontMarks;
            int marker = Legend.MarkerSize;
            int rowH = Math.Max(Legend.RowHeight, marker);
            int maxTextW = 0;

            foreach (var it in _panes)
            {
                var sz = graphics.MeasureString(it.Label, font);
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
                    x = _canvas_width - Legend.Padding - boxW;
                    y = Legend.Padding;
                    break;
                case LegendAnchor.TopLeft:
                    x = Legend.Padding;
                    y = Legend.Padding;
                    break;
                case LegendAnchor.BottomRight:
                    x = _canvas_width - Legend.Padding - boxW;
                    y = _canvas_height - Legend.Padding - boxH;
                    break;
                case LegendAnchor.BottomLeft:
                    x = Legend.Padding;
                    y = _canvas_height - Legend.Padding - boxH;
                    break;
            }

            var back = Color.FromArgb((int)(Legend.Opacity * 255), _style.ColorBack);
            using var backBrush = new SolidBrush(back);
            using var borderPen = new Pen(_style.ColorShapes, 1f);

            graphics.FillRectangle(backBrush, x, y, boxW, boxH);

            if (Legend.Border)
                graphics.DrawRectangle(borderPen, x, y, boxW, boxH);

            int cx = x + innerPad;
            int cy = y + innerPad;

            for (int idx = 0; idx < _panes.Count; idx++)
            {
                var it = _panes[idx];
                int my = cy + idx * rowH + rowH / 2;

                DrawLegendMarker(graphics, cx, my, marker, it.Color, it.GraphType, it.Depth);
                graphics.DrawString(it.Label, font, textBrush, cx + marker + Legend.MarkerGap, my - rowH / 2 + 1);
            }
        }
        /// <summary>
        /// Draws the plot title centered at the top inside the drawable area.
        /// </summary>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="title">Title text</param>
        /// <param name="sizeX">Total canvas width in pixels</param>
        /// <param name="sizeY">Total canvas height in pixels</param>
        /// <param name="dw">Horizontal padding/margin used by the layout</param>
        /// <param name="dh">Vertical padding/margin used by the layout</param>
        private void Paint_Title(Graphics graphics, string title, int sizeX, int sizeY, int dw, int dh)
        {
            using var format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            using var br = new SolidBrush(_style.ColorText);
            graphics.DrawString(title, _style.FontText, br, new PointF(sizeX / 2, dh / 2), format);
        }
        /// <summary>
        /// Draws the X-axis label centered below the plot area.
        /// </summary>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="xlabel">X-axis label text</param>
        /// <param name="sizeX">Total canvas width in pixels</param>
        /// <param name="sizeY">Total canvas height in pixels</param>
        /// <param name="dw">Horizontal padding/margin used by the layout</param>
        /// <param name="dh">Vertical padding/margin used by the layout</param>
        private void Paint_LabelX(Graphics graphics, string xlabel, int sizeX, int sizeY, int dw, int dh)
        {
            using var format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            var size = graphics.MeasureString(xlabel, _style.FontText);
            using var br = new SolidBrush(_style.ColorText);
            graphics.DrawString(xlabel, _style.FontText, br, new PointF(sizeX / 2, sizeY - dh / 2 + size.Height / 4), format);
        }
        /// <summary>
        /// Draws the Y-axis label centered at the left, rotated 90°.
        /// </summary>
        /// <param name="graphics">Target <see cref="Graphics"/> surface</param>
        /// <param name="ylabel">Y-axis label text</param>
        /// <param name="sizeX">Total canvas width in pixels</param>
        /// <param name="sizeY">Total canvas height in pixels</param>
        /// <param name="dw">Horizontal padding/margin used by the layout</param>
        /// <param name="dh">Vertical padding/margin used by the layout</param>
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

            using var br = new SolidBrush(_style.ColorText);
            graphics.TranslateTransform(cx, cy);
            graphics.RotateTransform(-90f);
            graphics.DrawString(ylabel, _style.FontText, br, PointF.Empty, fmt);
            graphics.Restore(state);
        }
        #endregion

        #region Helper voids
        /// <summary>
        /// Formats a numeric tick/label value with sane defaults for scientific vs fixed notation.
        /// </summary>
        /// <remarks>
        /// - Returns empty string for singular values (NaN/Inf).<br/>
        /// - Uses scientific notation for large magnitudes (≥ 1e4) or tiny nonzero magnitudes (&lt; 1e-3).<br/>
        /// - Otherwise prints up to three decimals.
        /// </remarks>
        /// <param name="v">Value to format</param>
        /// <returns>Formatted string for <paramref name="v"/>; empty for singular values</returns>
        private string GetNumString(float v)
        {
            if (Points.IsSingularPoint(v))
                return string.Empty;

            if (v == 0 /*|| Maths.Abs(v) < 1e-8f*/)
                return "0";

            float a = MathF.Abs(v);

            if (a >= 1e4f || (a > 0f && a < 1e-3f))
                return v.ToString("0.##E+0");

            return v.ToString("0.###");
        }
        #endregion

        #endregion
    }
}
