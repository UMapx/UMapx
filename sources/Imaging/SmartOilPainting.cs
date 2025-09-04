using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the smart oil painting filter.
    /// </summary>
    [Serializable]
    public unsafe class SmartOilPainting : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int r0;
        private int r1;
        private int levels;
        private float sharpening;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the smart oil painting filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="levels">Levels [2, 256]</param>
        /// <param name="sharpening">Sharpening (> 0)</param>
        public SmartOilPainting(int radius = 4, int levels = 64, float sharpening = 0.5f)
        {
            this.Size = new SizeInt(radius, radius);
            this.levels = Math.Max(1, Math.Min(256, levels));
            this.sharpening = sharpening;
        }
        /// <summary>
        /// Initializes the smart oil painting filter.
        /// </summary>
        /// <param name="height">Filter height</param>
        /// <param name="width">Filter width</param>
        /// <param name="levels">Levels [2, 256]</param>
        /// <param name="sharpening">Sharpening (> 0)</param>
        public SmartOilPainting(int width, int height, int levels = 64, float sharpening = 0.5f)
        {
            Size = new SizeInt(width, height);
            this.levels = Math.Max(1, Math.Min(256, levels));
            this.sharpening = sharpening;
        }
        /// <summary>
        /// Initializes the smart oil painting filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="levels">Levels [2, 256]</param>
        /// <param name="sharpening">Sharpening (> 0)</param>
        public SmartOilPainting(SizeInt size, int levels = 64, float sharpening = 0.5f)
        {
            Size = size;
            this.levels = Math.Max(1, Math.Min(256, levels));
            this.sharpening = sharpening;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(this.r0, this.r1);
            }
            set
            {
                this.r0 = value.Width;
                this.r1 = value.Height;
            }
        }
        /// <summary>
        /// Gets or sets number of levels [2, 256].
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
            }
        }
        /// <summary>
        /// Gets or sets sharpening value (> 0).
        /// </summary>
        public float Sharpening
        {
            get
            {
                return this.sharpening;
            }
            set
            {
                this.sharpening = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width;
            int height = bmData.Height;
            int stride = bmData.Stride;

            byte* src = (byte*)bmSrc.Scan0;
            byte* dst = (byte*)bmData.Scan0;

            Parallel.For(0, height, y =>
            {
                int[] intensityCount = new int[levels];
                int[] rSum = new int[levels];
                int[] gSum = new int[levels];
                int[] bSum = new int[levels];

                for (int x = 0; x < width; x++)
                {
                    Array.Clear(intensityCount, 0, levels);
                    Array.Clear(rSum, 0, levels);
                    Array.Clear(gSum, 0, levels);
                    Array.Clear(bSum, 0, levels);

                    int maxCount = 0;
                    int maxIndex = 0;

                    for (int dy = -r1; dy <= r1; dy++)
                    {
                        int ny = MathsF.Range(y + dy, 0, height - 1);
                        for (int dx = -r0; dx <= r0; dx++)
                        {
                            int nx = MathsF.Range(x + dx, 0, width - 1);

                            byte* p = &src[ny * stride + nx * 4];

                            int b = p[0];
                            int g = p[1];
                            int r = p[2];

                            int intensity = (r + g + b) * levels / (3 * 256);
                            intensity = MathsF.Range(intensity, 0, levels - 1);

                            intensityCount[intensity]++;
                            rSum[intensity] += r;
                            gSum[intensity] += g;
                            bSum[intensity] += b;

                            if (intensityCount[intensity] > maxCount)
                            {
                                maxCount = intensityCount[intensity];
                                maxIndex = intensity;
                            }
                        }
                    }

                    maxCount = maxCount == 0 ? 1 : maxCount;

                    int rr = rSum[maxIndex] / maxCount;
                    int gg = gSum[maxIndex] / maxCount;
                    int bb = bSum[maxIndex] / maxCount;

                    // Optional sharpening boost to enhance brush stroke edges
                    byte* orig = &src[y * stride + x * 4];
                    rr = (int)(rr + sharpening * (rr - orig[2]));
                    gg = (int)(gg + sharpening * (gg - orig[1]));
                    bb = (int)(bb + sharpening * (bb - orig[0]));

                    byte* d = &dst[y * stride + x * 4];
                    d[2] = (byte)MathsF.Range(rr, 0, 255);
                    d[1] = (byte)MathsF.Range(gg, 0, 255);
                    d[0] = (byte)MathsF.Range(bb, 0, 255);
                    d[3] = 255;
                }
            });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap current = BitmapFormat.Bitmap(bmData);
            Bitmap Src = (Bitmap)current.Clone();
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Src, bmSrc);
            Src.Dispose();
            current.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
            Src.Dispose();
        }
        #endregion
    }

}
