using System;
using System.Drawing;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the graph pane.
    /// </summary>
    [Serializable]
    public class GraphPane
    {
        #region Class components
        /// <summary>
        /// Initializes the graph pane.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Function</param>
        /// <param name="depth">Depth</param>
        /// <param name="color">Color</param>
        /// <param name="symbol">Symbol type</param>
        /// <param name="label">Label</param>
        /// <param name="pane">Pane type</param>
        public GraphPane(float[] x, float[] y, float depth, Color color, Symbol symbol, string label, Pane pane)
        {
            X = x;
            Y = y;
            Depth = depth;
            Color = color;
            Symbol = symbol;
            Label = label;
            Pane = pane;
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
        /// Gets or sets symbol type.
        /// </summary>
        public Symbol Symbol { get; set; }
        /// <summary>
        /// Gets or sets label.
        /// </summary>
        public string Label { get; set; }
        /// <summary>
        /// Gets or sets pane type.
        /// </summary>
        public Pane Pane { get; set; }
        #endregion
    }
}
