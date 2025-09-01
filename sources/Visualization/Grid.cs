using System;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines a grid.
    /// </summary>
    [Serializable]
    public class Grid
    {
        #region Class components
        /// <summary>
        /// Initializes the grid.
        /// </summary>
        public Grid() { }
        /// <summary>
        /// Show the grid or not.
        /// </summary>
        public bool Show { get; set; } = false;
        /// <summary>
        /// Gets or sets grid style.
        /// </summary>
        public GridStyle Style { get; set; } = GridStyle.Dashed;
        /// <summary>
        /// Gets or sets grid dash length.
        /// </summary>
        public float DashLength { get; set; } = 6f;
        /// <summary>
        /// Gets or sets grid gap length.
        /// </summary>
        public float GapLength { get; set; } = 4f;
        #endregion
    }
}
