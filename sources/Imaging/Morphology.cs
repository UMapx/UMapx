using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the morphology filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Mathematical_morphology
    /// </remarks>
    /// </summary>
    [Serializable]
    public class Morphology : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="mode">Morphology mode</param>
        public Morphology(int radius = 3, MorphologyMode mode = MorphologyMode.Median)
        {
            Size = new SizeInt(radius, radius);
            Mode = mode;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="mode">Morphology mode</param>
        public Morphology(int width, int height, MorphologyMode mode = MorphologyMode.Median)
        {
            Size = new SizeInt(width, height);
            Mode = mode;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="mode">Morphology mode</param>
        public Morphology(SizeInt size, MorphologyMode mode = MorphologyMode.Median)
        {
            Size = size;
            Mode = mode;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(rw, rh);
            }
            set
            {
                this.rw = value.Width;
                this.rh = value.Height;
            }
        }
        /// <summary>
        /// Gets or sets thresholds.
        /// </summary>
        public MorphologyMode Mode
        { 
            get; 
            set; 
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int stride = bmSrc.Stride;

            var ry = rh / 2;
            var rx = rw / 2;
            int windowHeight = 2 * ry + 1;
            int windowWidth = 2 * rx + 1;
            int windowSize = windowHeight * windowWidth;

            int rank = LinealgOptions.MorphologyHistogramFastFilter.GetFilterRank(Mode, windowSize);
            int range = byte.MaxValue + 1;

            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, y =>
            {
                int[] histB = new int[range];
                int[] histG = new int[range];
                int[] histR = new int[range];

                for (int dy = -ry; dy <= ry; dy++)
                {
                    int yi = MathF.Range(y + dy, 0, height - 1);

                    for (int dx = -rx; dx <= rx; dx++)
                    {
                        int xi = MathF.Range(0 + dx, 0, width - 1);
                        byte* pixel = src + yi * stride + xi * 4;

                        histB[pixel[0]]++;
                        histG[pixel[1]]++;
                        histR[pixel[2]]++;
                    }
                }

                byte* dstPixel = dst + y * stride;
                dstPixel[0] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histB, rank); // B
                dstPixel[1] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histG, rank); // G
                dstPixel[2] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histR, rank); // R

                for (int x = 1; x < width; x++)
                {
                    int outX = MathF.Range(x - rx - 1, 0, width - 1);
                    int inX = MathF.Range(x + rx, 0, width - 1);

                    for (int dy = -ry; dy <= ry; dy++)
                    {
                        int yi = MathF.Range(y + dy, 0, height - 1);

                        byte* outPixel = src + yi * stride + outX * 4;
                        byte* inPixel = src + yi * stride + inX * 4;

                        histB[outPixel[0]]--;
                        histG[outPixel[1]]--;
                        histR[outPixel[2]]--;

                        histB[inPixel[0]]++;
                        histG[inPixel[1]]++;
                        histR[inPixel[2]]++;
                    }

                    byte* pDst = dst + y * stride + x * 4;
                    pDst[0] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histB, rank); // B
                    pDst[1] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histG, rank); // G
                    pDst[2] = LinealgOptions.MorphologyHistogramFastFilter.GetHistogramRank(histR, rank); // R
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

        #region Public static
        /// <summary>
        /// Returns erosion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Erosion(int width, int height)
        {
            return new Morphology(width, height, MorphologyMode.Erosion);
        }
        /// <summary>
        /// Returns dilatation filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Dilatation(int width, int height)
        {
            return new Morphology(width, height, MorphologyMode.Dilatation);
        }
        /// <summary>
        /// Returns median filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <returns>Morphology</returns>
        public static Morphology Median(int width, int height)
        {
            return new Morphology(width, height, MorphologyMode.Median);
        }
        #endregion
    }
}
