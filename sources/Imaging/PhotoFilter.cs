using System;
using SkiaDrawing;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the color photo filter.
    /// </summary>
    [Serializable]
    public class PhotoFilter : IBitmapFilter
    {
        #region Private data
        private float s;
        private Color color;
        private IFloatMesh blendf;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        /// <param name="blendf">Blend function</param>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public PhotoFilter(IFloatMesh blendf, Color color, float strength = 0.5f)
        {
            BlendFunction = blendf;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public PhotoFilter(Color color, float strength = 0.5f)
        {
            BlendFunction = BlendMode.Pegtop;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        public PhotoFilter()
        {
            BlendFunction = BlendMode.Pegtop;
            Color = Color.White;
            Strength = 0.5f;
        }
        /// <summary>
        /// Gets or sets the blend function.
        /// </summary>
        public IFloatMesh BlendFunction
        {
            get
            {
                return this.blendf;
            }
            set
            {
                this.blendf = value;
            }
        }
        /// <summary>
        /// gets or sets the filter color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets filter strenght [0, 1].
        /// </summary>
        public float Strength
        {
            get
            {
                return this.s;
            }
            set
            {
                this.s = Maths.Float(value);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;
            int stride = bmData.Stride;

            // filter color
            float rf = color.R / 255.0f;
            float gf = color.G / 255.0f;
            float bf = color.B / 255.0f;

            // do job
            Parallel.For(0, height, y =>
            {
                // bitmap color
                float lm;
                float rb;
                float gb;
                float bb;

                int x, ystride, k;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    // luma
                    lm = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0f;

                    // blending
                    rb = this.blendf(lm, rf) * 255.0f;
                    gb = this.blendf(lm, gf) * 255.0f;
                    bb = this.blendf(lm, bf) * 255.0f;

                    // recording
                    p[k + 2] = Maths.Byte(rb * s + p[k + 2] * (1.0f - s));
                    p[k + 1] = Maths.Byte(gb * s + p[k + 1] * (1.0f - s));
                    p[k    ] = Maths.Byte(bb * s + p[k    ] * (1.0f - s));
                }
            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Initializes the cold filter (82).
        /// </summary>
        public static PhotoFilter Cold82
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 181, 255));
            }
        }
        /// <summary>
        /// Initializes the cold filter LBB.
        /// </summary>
        public static PhotoFilter ColdLBB
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 93, 255));
            }
        }
        /// <summary>
        /// Initializes the hot filter (81).
        /// </summary>
        public static PhotoFilter Warm81
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(235, 177, 19));
            }
        }
        /// <summary>
        /// Initializes the hot filter LBA.
        /// </summary>
        public static PhotoFilter WarmLBA
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(250, 150, 0));
            }
        }
        /// <summary>
        /// Initializes the sepia filter.
        /// </summary>
        public static PhotoFilter Sepia
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Initializes the red filter.
        /// </summary>
        public static PhotoFilter Red
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(234, 26, 26));
            }
        }
        /// <summary>
        /// Initializes the blue filter.
        /// </summary>
        public static PhotoFilter Blue
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(29, 53, 234));
            }
        }
        /// <summary>
        /// Initializes the green filter.
        /// </summary>
        public static PhotoFilter Green
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(25, 201, 25));
            }
        }
        /// <summary>
        /// Initializes the "underwater" filter.
        /// </summary>
        public static PhotoFilter Underwater
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 194, 177));
            }
        }
        /// <summary>
        /// Initializes the purple filter.
        /// </summary>
        public static PhotoFilter Purple
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(227, 24, 227));
            }
        }
        /// <summary>
        /// Initializes the orange filter.
        /// </summary>
        public static PhotoFilter Orange
        {
            get
            {
                return new PhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Initializes the yellow filter.
        /// </summary>
        public static PhotoFilter Yellow
        {
            get
            {
                return new PhotoFilter(Color.Yellow);
            }
        }
        #endregion
    }
}
