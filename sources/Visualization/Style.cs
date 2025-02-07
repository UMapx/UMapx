using System;
using SkiaDrawing;

namespace UMapx.Visualization
{
    /// <summary>
    /// Defines the figure style.
    /// </summary>
    [Serializable]
    public class Style : IDisposable
    {
        #region Private data

        private Font _fontMarks = new Font("Trojan Pro", 10, FontStyle.Regular);
        private Font _fontText = new Font("Trojan Pro", 12, FontStyle.Regular);

        #endregion

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
        public Font FontMarks
        {
            get
            {
                return _fontMarks;
            }
            set
            {
                if ((_fontMarks != null) && (!ReferenceEquals(_fontMarks, value)))
                {
                    _fontMarks.Dispose();
                }
                _fontMarks = value;
            }
        }
        /// <summary>
        /// Gets or sets text font.
        /// </summary>
        public Font FontText
        {
            get
            {
                return _fontText;
            }
            set
            {
                if ((_fontText != null) && (!ReferenceEquals(_fontText, value)))
                {
                    _fontText.Dispose();
                }
                _fontText = value;
            }
        }
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
                    _fontText?.Dispose();
                    _fontMarks?.Dispose();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~Style()
        {
            Dispose(false);
        }
        #endregion

        #region Completed styles
        /// <summary>
        /// Returns MATLAB style.
        /// </summary>
        public static Style MATLAB
        {
            get
            {
                Style style = new Style
                {
                    FontMarks = new Font("Helvetica", 10, FontStyle.Regular),
                    FontText = new Font("Helvetica", 12, FontStyle.Regular),
                    ColorFrame = Color.FromArgb(204, 204, 204),
                    ColorBack = Color.White,
                    ColorGrid = Color.FromArgb(230, 230, 230),
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Arial", 10, FontStyle.Regular),
                    FontText = new Font("Arial", 12, FontStyle.Regular),
                    ColorFrame = Color.White,
                    ColorBack = Color.White,
                    ColorGrid = Color.FromArgb(77, 201, 77),
                    ColorShapes = Color.FromArgb(77, 201, 77),
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Arial", 10, FontStyle.Regular),
                    FontText = new Font("Arial", 12, FontStyle.Regular),
                    ColorFrame = Color.White,
                    ColorBack = Color.White,
                    ColorGrid = Color.FromArgb(128, 128, 128),
                    ColorShapes = Color.FromArgb(128, 128, 128),
                    ColorText = Color.FromArgb(72, 0, 255),
                    ColorMarks = Color.FromArgb(72, 0, 255),
                    DepthShapes = 0.5f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Arial Cyr", 10, FontStyle.Regular),
                    FontText = new Font("Arial Cyr", 12, FontStyle.Regular),
                    ColorFrame = Color.White,
                    ColorBack = Color.FromArgb(219, 238, 244),
                    ColorGrid = Color.FromArgb(127, 137, 139),
                    ColorShapes = Color.FromArgb(127, 137, 139),
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.FromArgb(240, 240, 240),
                    ColorBack = Color.White,
                    ColorGrid = Color.FromArgb(180, 180, 180),
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.Beige,
                    ColorBack = Color.Bisque,
                    ColorGrid = Color.LightGray,
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.Lavender,
                    ColorBack = Color.LightCyan,
                    ColorGrid = Color.LightGray,
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.MistyRose,
                    ColorBack = Color.White,
                    ColorGrid = Color.LightGray,
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.Coral,
                    ColorBack = Color.LightCoral,
                    ColorGrid = Color.LightGray,
                    ColorShapes = Color.Black,
                    ColorText = Color.Black,
                    ColorMarks = Color.Black,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
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
                Style style = new Style
                {
                    FontMarks = new Font("Trojan Pro", 10, FontStyle.Regular),
                    FontText = new Font("Trojan Pro", 12, FontStyle.Regular),
                    ColorFrame = Color.FromArgb(19, 19, 19),
                    ColorBack = Color.FromArgb(35, 35, 35),
                    ColorGrid = Color.FromArgb(75, 75, 70),
                    ColorShapes = Color.White,
                    ColorText = Color.White,
                    ColorMarks = Color.White,
                    DepthShapes = 1.9f,
                    GridX = true,
                    GridY = true
                };
                return style;
            }
        }
        #endregion
    }
}
