using System;
using System.Drawing;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the figure style.
    /// </summary>
    [Serializable]
    public class Style
    {
        #region Initialize
        /// <summary>
        /// Initializes the figure style.
        /// </summary>
        public Style() { }
        #endregion

        #region Parameters
        /// <summary>
        /// Gets or sets frame color.
        /// </summary>
        public Color ColorFrame { get; set; } = Color.FromArgb(240, 240, 240);
        /// <summary>
        /// Gets or sets background color.
        /// </summary>
        public Color ColorBack { get; set; } = Color.White;
        /// <summary>
        /// Gets or sets grid color.
        /// </summary>
        public Color ColorGrid { get; set; } = Color.FromArgb(180, 180, 180);
        /// <summary>
        /// Gets or sets shapes color.
        /// </summary>
        public Color ColorShapes { get; set; } = Color.Black;
        /// <summary>
        /// Gets or sets text color.
        /// </summary>
        public Color ColorText { get; set; } = Color.Black;
        /// <summary>
        /// Gets or sets marks color.
        /// </summary>
        public Color ColorMarks { get; set; } = Color.Black;
        /// <summary>
        /// Gets or sets marks font.
        /// </summary>
        public Font FontMarks { get; set; } = new Font("Trojan Pro", 10, FontStyle.Regular);
        /// <summary>
        /// Gets or sets text font.
        /// </summary>
        public Font FontText { get; set; } = new Font("Trojan Pro", 12, FontStyle.Regular);
        /// <summary>
        /// Gets or sets shapes depth.
        /// </summary>
        public float DepthShapes { get; set; } = 2f;
        /// <summary>
        /// Gets or sets X grid.
        /// </summary>
        public bool GridX { get; set; } = true;
        /// <summary>
        /// Gets or sets Y grid.
        /// </summary>
        public bool GridY { get; set; } = true;
        #endregion

        #region Completed styles
        /// <summary>
        /// Returns MATLAB style.
        /// </summary>
        public static Style MATLAB
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Helvetica", 10, FontStyle.Regular);
                style.FontText = new Font("Helvetica", 12, FontStyle.Regular);
                style.ColorFrame = Color.FromArgb(204, 204, 204);
                style.ColorBack = Color.White;
                style.ColorGrid = Color.FromArgb(230, 230, 230);
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns MathCad style.
        /// </summary>
        public static Style MathCad
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Arial", 10, FontStyle.Regular);
                style.FontText = new Font("Arial", 12, FontStyle.Regular);
                style.ColorFrame = Color.White;
                style.ColorBack = Color.White;
                style.ColorGrid = Color.FromArgb(77, 201, 77);
                style.ColorShapes = Color.FromArgb(77, 201, 77);
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns MicroCap style.
        /// </summary>
        public static Style MicroCap
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Arial", 10, FontStyle.Regular);
                style.FontText = new Font("Arial", 12, FontStyle.Regular);
                style.ColorFrame = Color.White;
                style.ColorBack = Color.White;
                style.ColorGrid = Color.FromArgb(128, 128, 128);
                style.ColorShapes = Color.FromArgb(128, 128, 128);
                style.ColorText = Color.FromArgb(72, 0, 255);
                style.ColorMarks = Color.FromArgb(72, 0, 255);
                style.DepthShapes = 0.5f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns Excel style.
        /// </summary>
        public static Style Excel
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Arial Cyr", 10, FontStyle.Regular);
                style.FontText = new Font("Arial Cyr", 12, FontStyle.Regular);
                style.ColorFrame = Color.White;
                style.ColorBack = Color.FromArgb(219, 238, 244);
                style.ColorGrid = Color.FromArgb(127, 137, 139);
                style.ColorShapes = Color.FromArgb(127, 137, 139);
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns standart style.
        /// </summary>
        public static Style Standart
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.FromArgb(240, 240, 240);
                style.ColorBack = Color.White;
                style.ColorGrid = Color.FromArgb(180, 180, 180);
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns beige style.
        /// </summary>
        public static Style Beige
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.Beige;
                style.ColorBack = Color.Bisque;
                style.ColorGrid = Color.LightGray;
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns cyan style.
        /// </summary>
        public static Style Cyan
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.Lavender;
                style.ColorBack = Color.LightCyan;
                style.ColorGrid = Color.LightGray;
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns rose style.
        /// </summary>
        public static Style Rose
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.MistyRose;
                style.ColorBack = Color.White;
                style.ColorGrid = Color.LightGray;
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns coral style.
        /// </summary>
        public static Style Coral
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.Coral;
                style.ColorBack = Color.LightCoral;
                style.ColorGrid = Color.LightGray;
                style.ColorShapes = Color.Black;
                style.ColorText = Color.Black;
                style.ColorMarks = Color.Black;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        /// <summary>
        /// Returns black style.
        /// </summary>
        public static Style Black
        {
            get
            {
                Style style = new Style();
                style.FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
                style.FontText = new Font("Trojan Pro", 12, FontStyle.Regular);
                style.ColorFrame = Color.FromArgb(19, 19, 19);
                style.ColorBack = Color.FromArgb(35, 35, 35);
                style.ColorGrid = Color.FromArgb(75, 75, 70);
                style.ColorShapes = Color.White;
                style.ColorText = Color.White;
                style.ColorMarks = Color.White;
                style.DepthShapes = 1.9f;
                style.GridX = true;
                style.GridY = true;
                return style;
            }
        }
        #endregion
    }
}
