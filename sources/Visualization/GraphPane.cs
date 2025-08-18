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
        public GraphPane() { }
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
        /// Gets or sets type.
        /// </summary>
        public Symbol Type { get; set; }
        #endregion
    }
}
