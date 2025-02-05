using System;
using System.Collections.Generic;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the figure to plotting in a Cartesian coordinate system.
    /// </summary>
    [Serializable]
    public class Figure : IDisposable
    {
        #region Panes
        /// <summary>
        /// Plot panes.
        /// </summary>
        protected List<GraphPane> PlotPanes = new List<GraphPane>();
        /// <summary>
        /// Stem panes.
        /// </summary>
        protected List<GraphPane> StemPanes = new List<GraphPane>();
        /// <summary>
        /// Scatter panes.
        /// </summary>
        protected List<GraphPane> ScatterPanes = new List<GraphPane>();
        /// <summary>
        /// Image pane.
        /// </summary>
        protected Bitmap ImagePane;
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

        #region Private data
        private Style _style = new Style();
        private int figure_width, figure_height;
        private int canvas_width, canvas_height;
        private float xmin = -5, xmax = 5, ymin = -5, ymax = 5;
        private int xscale = 10, yscale = 10;
        #endregion

        #region Figure properties
        /// <summary>
        /// Gets or sets figure style.
        /// </summary>
        public Style Style
        {
            get { return _style; }
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
                    throw new Exception("Start and end points cannot be the same");

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
                    throw new Exception("Start and end points cannot be the same");

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

        #region Figure options and graph options
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
            using var canvas = new Bitmap(canvas_width, canvas_height);
            using Graphics figure_graphics = Graphics.FromImage(figure); figure_graphics.Clear(Style.ColorFrame);
            using Graphics canvas_grpahics = Graphics.FromImage(canvas); canvas_grpahics.Clear(Style.ColorBack);

            // offsets
            int dx = (canvas_width / xscale),           dy = (canvas_height / yscale);
            int dw = (figure_width - canvas_width) / 2, dh = (figure_height - canvas_height) / 2;
            #endregion

            #region Autorange
            //autorange or not?
            if (ImagePane is object)
            {
                xmin = ymin = 0;
                xmax = size.Width;
                ymax = size.Height;
            }
            else if (AutoRange)
            {
                xmin = ymin = float.MaxValue;
                xmax = ymax = float.MinValue;

                foreach (var pane in PlotPanes)
                {
                    xmin = Math.Min(xmin, pane.X.Min());
                    xmax = Math.Max(xmax, pane.X.Max());
                    ymin = Math.Min(ymin, pane.Y.Min());
                    ymax = Math.Max(ymax, pane.Y.Max());
                }

                foreach (var pane in StemPanes)
                {
                    xmin = Math.Min(xmin, pane.X.Min());
                    xmax = Math.Max(xmax, pane.X.Max());
                    ymin = Math.Min(ymin, pane.Y.Min());
                    ymax = Math.Max(ymax, pane.Y.Max());
                }

                foreach (var pane in ScatterPanes)
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
                    canvas_grpahics.DrawLine(pen1, xpoint, 0, xpoint, canvas_height);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, ymin, ymax, canvas_height);
                    canvas_grpahics.DrawLine(pen1, 0, ypoint, canvas_width, ypoint);
                }
            }
            #endregion

            #region Title and labels
            Paint_Title(figure_graphics, Title, figure_width, figure_height, dw, dh);
            Paint_LabelX(figure_graphics, LabelX, figure_width, figure_height, dw, dh);
            Paint_LabelY(figure_graphics, LabelY, figure_width, figure_height, dw, dh);
            #endregion

            #region Graph painting
            if (ImagePane is object)
            {
                // 2D plotting
                Rectangle rectangle = new Rectangle(0, 0, canvas_width, canvas_height);
                canvas_grpahics.DrawImage(ImagePane, rectangle);
            }
            //else
            {
                // 1D plotting
                float r0 = 4;
                float radius;

                // Plotting:
                foreach (var current in PlotPanes)
                {
                    radius = (current.Depth + r0) * 2;

                    switch (current.Type)
                    {
                        case Symbol.None:
                            PlotLine(canvas_grpahics, current.X, current.Y, current.Depth, current.Color);
                            break;

                        case Symbol.Circle:
                            PlotCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Ball:
                            PlotCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;

                        case Symbol.Rectangle:
                            PlotRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Polygon:
                            PlotRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;
                    }
                }
                // Stemming:
                foreach (var current in StemPanes)
                {
                    radius = (current.Depth + r0) * 2;

                    switch (current.Type)
                    {
                        case Symbol.None:
                            StemLine(canvas_grpahics, current.X, current.Y, current.Depth, current.Color);
                            break;

                        case Symbol.Circle:
                            StemCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Ball:
                            StemCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;

                        case Symbol.Rectangle:
                            StemRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Polygon:
                            StemRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;
                    }
                }
                // Scattering:
                foreach (var current in ScatterPanes)
                {
                    radius = (current.Depth + r0) * 2;

                    switch (current.Type)
                    {
                        case Symbol.None:
                            ScatterLine(canvas_grpahics, current.X, current.Y, current.Depth, current.Color);
                            break;

                        case Symbol.Circle:
                            ScatterCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Ball:
                            ScatterCircle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;

                        case Symbol.Rectangle:
                            ScatterRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, false);
                            break;

                        case Symbol.Polygon:
                            ScatterRectangle(canvas_grpahics, current.X, current.Y, current.Depth, current.Color, radius, true);
                            break;
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
                    canvas_grpahics.DrawLine(pen2, xpoint, canvas_height, xpoint, canvas_height - s);
                    canvas_grpahics.DrawLine(pen2, xpoint, 0, xpoint, s);
                }

                ylength = Y.Length;

                for (i = 0; i < ylength; i++)
                {
                    numY = Y[i];
                    ypoint = (int)Points.Y2Point(numY, ymin, ymax, canvas_height);
                    canvas_grpahics.DrawLine(pen2, 0, ypoint, s, ypoint);
                    canvas_grpahics.DrawLine(pen2, canvas_width, ypoint, canvas_width - s, ypoint);
                }
            }

            canvas_grpahics.DrawRectangle(pen2, 0, 0, canvas_width - 1, canvas_height - 1);
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
            ImagePane = (Bitmap)bitmap.Clone();
        }
        /// <summary>
        /// Adds graph pane to continuous plot.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="depth">Depth</param>
        /// <param name="color">Color</param>
        /// <param name="type">Symbol type</param>
        public void Plot(float[] x, float[] y, float depth, Color color, Symbol type)
        {
            if (x.Length != y.Length)
                throw new Exception("Vectors must be of the same length");

            var pane = new GraphPane()
            {
                X = x,
                Y = y,
                Depth = depth,
                Color = color,
                Type = type,
            };

            PlotPanes.Add(pane);
        }
        /// <summary>
        /// Adds graph pane to continuous plot.
        /// </summary>
        /// <param name="pane">Graph pane</param>
        public void Plot(GraphPane pane)
        {
            if (pane.X.Length != pane.Y.Length)
                throw new Exception("Vectors must be of the same length");

            PlotPanes.Add(pane);
        }
        /// <summary>
        /// Adds graph pane to stem plot.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="depth">Depth</param>
        /// <param name="color">Color</param>
        /// <param name="type">Symbol type</param>
        public void Stem(float[] x, float[] y, float depth, Color color, Symbol type)
        {
            if (x.Length != y.Length)
                throw new Exception("Vectors must be of the same length");

            var pane = new GraphPane()
            {
                X = x,
                Y = y,
                Depth = depth,
                Color = color,
                Type = type,
            };

            StemPanes.Add(pane);
        }
        /// <summary>
        /// Adds graph pane to stem plot.
        /// </summary>
        /// <param name="pane">Graph pane</param>
        public void Stem(GraphPane pane)
        {
            if (pane.X.Length != pane.Y.Length)
                throw new Exception("Vectors must be of the same length");

            StemPanes.Add(pane);
        }
        /// <summary>
        /// Add graph pane to scatter plot.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="depth">Depth</param>
        /// <param name="color">Color</param>
        /// <param name="type">Symbol type</param>
        public void Scatter(float[] x, float[] y, float depth, Color color, Symbol type)
        {
            if (x.Length != y.Length)
                throw new Exception("Vectors must be of the same length");

            var pane = new GraphPane()
            {
                X = x,
                Y = y,
                Depth = depth,
                Color = color,
                Type = type,
            };

            ScatterPanes.Add(pane);
        }
        /// <summary>
        /// Adds graph pane to scatter plot.
        /// </summary>
        /// <param name="pane">Graph pane</param>
        public void Scatter(GraphPane pane)
        {
            if (pane.X.Length != pane.Y.Length)
                throw new Exception("Vectors must be of the same length");

            ScatterPanes.Add(pane);
        }
        /// <summary>
        /// Removes all graphs from the figure. 
        /// </summary>
        public void Clear()
        {
            PlotPanes.Clear();
            StemPanes.Clear();
            ScatterPanes.Clear();
            ImagePane?.Dispose();
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
            return;
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
            return;
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
            return;
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

            return;
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

            return;
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

            return;
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

            return;
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

            return;
        }
        #endregion

        #region Painter voids
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
            StringFormat format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            using var br = new SolidBrush(Style.ColorText);
            graphics.DrawString(title, Style.FontText, br, new PointF(sizeX / 2, dh / 2), format);
            format.Dispose();
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
            StringFormat format = new StringFormat
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            SizeF size = graphics.MeasureString(xlabel, Style.FontText);
            using var br = new SolidBrush(Style.ColorText);
            graphics.DrawString(xlabel, Style.FontText, br, new PointF(sizeX / 2, sizeY - dh / 2 + size.Height / 4), format);
            format.Dispose();
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
            StringFormat format = new StringFormat(StringFormatFlags.DirectionVertical)
            {
                LineAlignment = StringAlignment.Center,
                Alignment = StringAlignment.Center
            };

            graphics.TranslateTransform(sizeX / 2, sizeY / 2);
            graphics.RotateTransform(180);
            graphics.TranslateTransform(-sizeX / 2, -sizeY / 2);

            SizeF size = graphics.MeasureString(ylabel, Style.FontText);
            using var br = new SolidBrush(Style.ColorText);
            graphics.DrawString(ylabel, Style.FontText, br, new PointF(sizeX - dw / 2 + size.Height / 4, sizeY / 2), format);
            format.Dispose();
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
                    ImagePane?.Dispose();
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
