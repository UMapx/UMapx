using System;
using System.Drawing;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the plot series.
    /// </summary>
    [Serializable]
    public class PlotSeries
    {
        #region Class components
        /// <summary>
        /// Initializes the plot series.
        /// </summary>
        /// <param name="x">Value</param>
        /// <param name="y">Function</param>
        /// <param name="depth">Depth</param>
        /// <param name="color">Color</param>
        /// <param name="shapeType">Shape type</param>
        /// <param name="label">Label</param>
        /// <param name="seriesType">Series type</param>
        public PlotSeries(float[] x, float[] y, float depth, Color color, SeriesType seriesType, ShapeType shapeType, string label)
        {
            X = x;
            Y = y;
            Depth = depth;
            Color = color;
            ShapeType = shapeType;
            Label = label;
            SeriesType = seriesType;
        }
        /// <summary>
        /// Gets or sets argument array.
        /// </summary>
        public float[] X { get; set; }
        /// <summary>
        /// Gets or sets function array.
        /// </summary>
        public float[] Y { get; set; }
        /// <summary>
        /// Gets or sets depth.
        /// </summary>
        public float Depth { get; set; }
        /// <summary>
        /// Gets or sets color.
        /// </summary>
        public Color Color { get; set; }
        /// <summary>
        /// Gets or sets shape type.
        /// </summary>
        public ShapeType ShapeType { get; set; }
        /// <summary>
        /// Gets or sets label.
        /// </summary>
        public string Label { get; set; }
        /// <summary>
        /// Gets or sets series type.
        /// </summary>
        public SeriesType SeriesType { get; set; }
        #endregion
    }
}
