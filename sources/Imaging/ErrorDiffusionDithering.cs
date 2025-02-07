using System;
using SkiaDrawing;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the error diffusion dithering filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ErrorDiffusionDithering : Rebuilder, IBitmapFilter
    {
        #region Private data
        private int x;
        private int y;
        private int width;
        private int height;
        private int stride;
        private float[][] matrix;
        private float summary;
        private int levels;
        private float[] table;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the error diffusion dithering filter.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="matrix">Matrix</param>
        public ErrorDiffusionDithering(int levels, float[][] matrix)
        {
            this.Levels = levels;
            this.Matrix = matrix;
        }
        /// <summary>
        /// Gets or sets the number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return this.levels;
            }
            set
            {
                this.levels = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the matrix.
        /// </summary>
        public float[][] Matrix
        {
            get
            {
                return this.matrix;
            }
            set
            {
                this.matrix = value;
                this.summary = 0;
                int n = matrix.Length;
                float[] row;
                int i, j, k;

                for (i = 0; i < n; i++)
                {
                    row = matrix[i];
                    k = row.Length;

                    for (j = 0; j < k; j++)
                    {
                        summary += row[j];
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Bitmap</returns>
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            // params
            this.width = bmData.Width;
            this.height = bmData.Height;
            this.stride = bmData.Stride;
            int length = table.Length;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int r, g, b;
            Color color;

            // for each line
            for (y = 0; y < height; y++)
            {
                // for each pixels
                for (x = 0; x < width; x++, p += 4)
                {
                    // current
                    r = p[2];
                    g = p[1];
                    b = p[0];

                    // get color from palette, which is the closest to current pixel's value
                    color = GetColor(r, g, b, table);
                    p[2] = color.R;
                    p[1] = color.G;
                    p[0] = color.B;

                    // do error diffusion
                    Diffuse(r - color.R, g - color.G, b - color.B, p);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Bitmap</returns>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            Apply(bmData);
            BitmapFormat.Unlock(Data, bmData);
            return;
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.table = Intensity.Quantize(this.levels, 256).Mul(255);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="rError"></param>
        /// <param name="gError"></param>
        /// <param name="bError"></param>
        /// <param name="ptr"></param>
        protected unsafe void Diffuse(int rError, int gError, int bError, byte* ptr)
        {
            float edR;	// error diffusion
            float edG;	// error diffusion
            float edB;	// error diffusion

            // do error diffusion to right-standing neighbors
            float[] row = matrix[0];
            int jI, jP, i, k, jC, n;
            int length = matrix.Length;


            for (jI = 1, jP = 4, jC = 0, k = row.Length; jC < k; jI++, jC++, jP += 4)
            {
                if (x + jI >= width)
                    break;

                edR = ptr[jP + 2] + (rError * row[jC]) / summary;
                ptr[jP + 2] = Maths.Byte(edR);

                edG = ptr[jP + 1] + (gError * row[jC]) / summary;
                ptr[jP + 1] = Maths.Byte(edG);

                edB = ptr[jP + 0] + (bError * row[jC]) / summary;
                ptr[jP + 0] = Maths.Byte(edB);
            }

            // do error diffusion to bottom neigbors
            for (i = 1, n = length; i < n; i++)
            {
                if (y + i >= height)
                    break;

                // move pointer to next image line
                ptr += stride;

                // get coefficients of the row
                row = matrix[i];

                // process the row
                for (jC = 0, k = row.Length, jI = -(k >> 1), jP = -(k >> 1) * 4; jC < k; jI++, jC++, jP += 4)
                {
                    if (x + jI >= width)
                        break;
                    if (x + jI < 0)
                        continue;

                    edR = ptr[jP + 2] + (rError * row[jC]) / summary;
                    ptr[jP + 2] = Maths.Byte(edR);

                    edG = ptr[jP + 1] + (gError * row[jC]) / summary;
                    ptr[jP + 1] = Maths.Byte(edG);

                    edB = ptr[jP + 0] + (bError * row[jC]) / summary;
                    ptr[jP + 0] = Maths.Byte(edB);

                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="red"></param>
        /// <param name="green"></param>
        /// <param name="blue"></param>
        /// <param name="table"></param>
        /// <returns></returns>
        private Color GetColor(int red, int green, int blue, float[] table)
        {
            byte r = Maths.Byte(table[red]);
            byte g = Maths.Byte(table[green]);
            byte b = Maths.Byte(table[blue]);
            return Color.FromArgb(r, g, b);
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Initializes the Atkinson dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Atkinson
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[3][] {
                new float[2]   {          1, 1 },
                new float[5]   { 0, 1, 1, 1, 0 },
                new float[5]   { 0, 0, 1, 0, 0 } });
            }
        }
        /// <summary>
        /// Initializes the Burkes dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Burkes
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                    new float[2] { 8, 4 },
                    new float[5] { 2, 4, 8, 4, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Fan dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Fan
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                    new float[1] { 8 },
                    new float[4] { 1, 1, 2, 4 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra lite dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering SierraLite
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                new float[1] { 2 },
                new float[2] { 1, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Sierra
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[3][] {
                new float[2] { 5, 3 },
                new float[5] { 2, 4, 5, 4, 2 },
                new float[3] { 2, 3, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra lite dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering SierraTowsRows
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                new float[2] { 4, 3 },
                new float[5] { 1, 2, 3, 2, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Flyd-Steinberg dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering FloydSteinberg
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                new float[1] {       7 },
                new float[3] { 3, 5, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Jarvis-Judice-Ninke dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering JarvisJudiceNinke
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[3][] {
                new float[2] {          7, 5 },
                new float[5] { 3, 5, 7, 5, 3 },
                new float[5] { 1, 3, 5, 3, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Stevenson dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Stevenson
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[3][] {
                new float[4] { 12, 26, 30, 16 },
                new float[3] { 12, 26, 12   },
                new float[4] { 5, 12, 12, 5 } });
            }
        }
        /// <summary>
        /// Initializes the Shiau dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Shiau
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[2][] {
                new float[1] { 4 },
                new float[3] { 1, 1, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Stucki dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Stucki
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new float[3][] {
                new float[2] { 8, 4 },
                new float[5] { 2, 4, 8, 4, 2 },
                new float[5] { 1, 2, 4, 2, 1 } });
            }
        }
        #endregion
    }
}
