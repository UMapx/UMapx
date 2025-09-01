using System;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines a legend.
    /// </summary>
    [Serializable]
    public class Legend
    {
        #region Class components
        /// <summary>
        /// Initializes the legend.
        /// </summary>
        public Legend() { }
        /// <summary>
        /// Show the legend panel or not.
        /// </summary>
        public bool Show { get; set; } = true;
        /// <summary>
        /// Gets or sets legend anchor inside the canvas.
        /// </summary>
        public LegendAnchor Anchor { get; set; } = LegendAnchor.BottomRight;
        /// <summary>
        /// Gets or sets outer padding inside canvas (px).
        /// </summary>
        public int Padding { get; set; } = 10;
        /// <summary>
        /// Gets or sets spacing between marker and text (px).
        /// </summary>
        public int MarkerGap { get; set; } = 8;
        /// <summary>
        /// Gets or sets row height for each legend item (px).
        /// </summary>
        public int RowHeight { get; set; } = 18;
        /// <summary>
        /// Gets or sets marker nominal size (px).
        /// </summary>
        public int MarkerSize { get; set; } = 12;
        /// <summary>
        /// Gets or sets background opacity [0, 1].
        /// </summary>
        public float Opacity { get; set; } = 0.85f;
        /// <summary>
        /// Gets or sets draw legend border.
        /// </summary>
        public bool Border { get; set; } = true;
        #endregion
    }
}
