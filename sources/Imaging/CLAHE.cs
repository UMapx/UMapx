using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Contrast Limited Adaptive Histogram Equalization (CLAHE).
    /// Operates in RGB by equalizing luminance (Rec.709) with bilinear interpolation between tile LUTs.
    /// <remarks>
    /// More information can be found on the website:
    /// https://uk.mathworks.com/help/visionhdl/ug/contrast-adaptive-histogram-equalization.html
    /// </remarks>
    /// </summary>
    [Serializable]
    public class CLAHE : IBitmapFilter
    {
        #region Private data
        private int tilesX;
        private int tilesY;
        private float clipLimit; // relative factor: CL = clipLimit * (tileArea / bins)
        private const int Bins = 256;
        private const float EPS = 1e-6f;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes CLAHE filter.
        /// </summary>
        /// <param name="tilesX">Number of tiles horizontally (>=1)</param>
        /// <param name="tilesY">Number of tiles vertically (>=1)</param>
        /// <param name="clipLimit">Contrast limit factor (>=0). 0 disables clipping</param>
        public CLAHE(int tilesX = 8, int tilesY = 8, float clipLimit = 2.0f)
        {
            TilesX = tilesX;
            TilesY = tilesY;
            ClipLimit = clipLimit;
        }
        /// <summary>
        /// Gets or sets tiles horizontally.
        /// </summary>
        public int TilesX
        {
            get => tilesX;
            set => tilesX = Math.Max(1, value);
        }
        /// <summary>
        /// Gets or sets tiles vertically.
        /// </summary>
        public int TilesY
        {
            get => tilesY;
            set => tilesY = Math.Max(1, value);
        }
        /// <summary>
        /// Contrast limit factor. Absolute clip is: clipLimit * (tileArea / Bins).
        /// </summary>
        public float ClipLimit
        {
            get => clipLimit;
            set => clipLimit = Math.Max(0f, value);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width;
            int height = bmData.Height;
            int stride = bmData.Stride;
            byte* basePtr = (byte*)bmData.Scan0.ToPointer();

            // Compute tile boundaries and centers
            var xBreaks = MakeBreaks(width, tilesX);
            var yBreaks = MakeBreaks(height, tilesY);
            var xCenters = MakeCenters(xBreaks);
            var yCenters = MakeCenters(yBreaks);

            // Build LUT per tile (luminance-based)
            // luts[tileIndex][bin] => byte
            byte[][] luts = new byte[tilesX * tilesY][];
            Parallel.For(0, tilesY, ty =>
            {
                for (int tx = 0; tx < tilesX; tx++)
                {
                    int x0 = xBreaks[tx];
                    int x1 = xBreaks[tx + 1]; // exclusive
                    int y0 = yBreaks[ty];
                    int y1 = yBreaks[ty + 1]; // exclusive
                    luts[ty * tilesX + tx] = BuildTileLUT(basePtr, stride, x0, x1, y0, y1, clipLimit);
                }
            });

            // Precompute interpolation indices/weights per column/row
            var (ix0, ix1, wx0, wx1) = PrecomputeAxis(width, xCenters);
            var (iy0, iy1, wy0, wy1) = PrecomputeAxis(height, yCenters);

            // Apply with bilinear interpolation in luminance space
            for (int y = 0; y < height; y++)
            {
                byte* row = basePtr + y * stride;
                int ty0 = iy0[y], ty1 = iy1[y];
                float wy0v = wy0[y], wy1v = wy1[y];

                for (int x = 0; x < width; x++)
                {
                    byte* px = row + (x << 2);

                    // Read B,G,R
                    float b = px[0];
                    float g = px[1];
                    float r = px[2];

                    // Luminance (Rec.709)
                    float Y = 0.2126f * r + 0.7152f * g + 0.0722f * b;
                    int bin = (int)(Y + 0.5f);
                    if (bin < 0) bin = 0; else if (bin > 255) bin = 255;

                    int tx0 = ix0[x], tx1 = ix1[x];
                    float wx0v = wx0[x], wx1v = wx1[x];

                    // gather 4 LUT values
                    byte[] lut00 = luts[ty0 * tilesX + tx0];
                    byte[] lut10 = luts[ty0 * tilesX + tx1];
                    byte[] lut01 = luts[ty1 * tilesX + tx0];
                    byte[] lut11 = luts[ty1 * tilesX + tx1];

                    float y00 = lut00[bin];
                    float y10 = lut10[bin];
                    float y01 = lut01[bin];
                    float y11 = lut11[bin];

                    // bilinear interpolation
                    float top = wx0v * y00 + wx1v * y10;
                    float bot = wx0v * y01 + wx1v * y11;
                    float Yeq = wy0v * top + wy1v * bot;

                    // Apply gain in RGB to preserve hue
                    if (Y < EPS)
                    {
                        // If near black, push towards grayscale mapped luminance
                        byte yb = (byte)Maths.Byte(Yeq);
                        px[0] = yb; px[1] = yb; px[2] = yb;
                    }
                    else
                    {
                        float gain = Yeq / Y;
                        px[0] = (byte)Maths.Byte(b * gain);
                        px[1] = (byte)Maths.Byte(g * gain);
                        px[2] = (byte)Maths.Byte(r * gain);
                    }
                    // alpha unchanged (px[3])
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            var bmData = BitmapFormat.Lock32bpp(Data);
            try { Apply(bmData); }
            finally { BitmapFormat.Unlock(Data, bmData); }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Computes tile breakpoints (inclusive/exclusive borders) along one axis,
        /// dividing <paramref name="size"/> pixels into <paramref name="tiles"/> nearly equal tiles.
        /// Remainder pixels are distributed by adding +1 to the first <c>rem</c> tiles.
        /// Example: size=10, tiles=3 → breaks: [0, 4, 7, 10].
        /// </summary>
        /// <param name="size">Axis length in pixels (width or height).</param>
        /// <param name="tiles">Number of tiles along this axis (>= 1).</param>
        /// <returns>Array of length <c>tiles+1</c> with cumulative borders; segment i is [br[i], br[i+1]).</returns>
        private static int[] MakeBreaks(int size, int tiles)
        {
            var br = new int[tiles + 1];
            int baseLen = size / tiles;   // minimum tile size
            int rem = size % tiles;       // how many tiles get +1
            int acc = 0;
            br[0] = 0;

            for (int i = 0; i < tiles; i++)
            {
                // Distribute remainder pixels to the first 'rem' tiles
                int w = baseLen + (i < rem ? 1 : 0);
                acc += w;
                br[i + 1] = acc;
            }
            return br;
        }
        /// <summary>
        /// Converts tile breakpoints to tile center coordinates on the pixel grid.
        /// Center is computed as the midpoint between inclusive pixel indices:
        /// mid([a, b)) ≈ (a + b - 1) / 2.
        /// </summary>
        /// <param name="br">Breaks as produced by <see cref="MakeBreaks"/>.</param>
        /// <returns>Array of tile center positions (float) aligned to pixel coordinates.</returns>
        private static float[] MakeCenters(int[] br)
        {
            int tiles = br.Length - 1;
            var c = new float[tiles];

            for (int i = 0; i < tiles; i++)
                // Example: segment [2, 6) covers pixels 2..5, center is (2+5)/2 = 3.5 → (2 + 6 - 1) * 0.5
                c[i] = 0.5f * (br[i] + br[i + 1] - 1);

            return c;
        }
        /// <summary>
        /// Precomputes, for each pixel position along an axis, the indices of the left/right neighboring
        /// tile centers and their interpolation weights. This is used for bilinear interpolation across tiles.
        /// Complexity: O(n), no per-pixel binary search required.
        /// </summary>
        /// <param name="n">Axis length in pixels (width or height).</param>
        /// <param name="centers">Tile center coordinates from <see cref="MakeCenters"/> (sorted ascending).</param>
        /// <returns>
        /// Tuple of:
        /// <list type="bullet">
        /// <item><description><c>i0[p]</c> – index of left tile center for pixel p</description></item>
        /// <item><description><c>i1[p]</c> – index of right tile center for pixel p (i0 or i0+1)</description></item>
        /// <item><description><c>w0[p]</c> – weight for left center</description></item>
        /// <item><description><c>w1[p]</c> – weight for right center (1 - w0)</description></item>
        /// </list>
        /// </returns>
        private static (int[] i0, int[] i1, float[] w0, float[] w1) PrecomputeAxis(int n, float[] centers)
        {
            int T = centers.Length;
            var i0 = new int[n];
            var i1 = new int[n];
            var w0 = new float[n];
            var w1 = new float[n];

            for (int p = 0; p < n; p++)
            {
                // Coarse guess of the segment index by relative position to avoid full search
                float pos = (float)p / Math.Max(1, n - 1);
                int guess = (int)(pos * (T - 1) + 0.5f);
                if (guess < 0) guess = 0; else if (guess > T - 1) guess = T - 1;

                // Adjust guess so that centers[left] <= p <= centers[right]
                int left = guess;
                if (p < centers[left] && left > 0) left--;
                while (left < T - 1 && p > centers[left + 1]) left++;

                int right = Math.Min(left + 1, T - 1);

                float cL = centers[left];
                float cR = centers[right];

                // Linear interpolation factor between neighboring centers
                // EPS prevents division by zero when centers coincide (single-tile case)
                float fx = (right == left) ? 0f : (p - cL) / Math.Max(EPS, (cR - cL));
                if (fx < 0f) fx = 0f; else if (fx > 1f) fx = 1f;

                i0[p] = left;
                i1[p] = right;
                w1[p] = fx;
                w0[p] = 1f - fx;
            }
            return (i0, i1, w0, w1);
        }
        /// <summary>
        /// Builds a CLAHE lookup table (LUT) for a single tile region:
        /// 1) computes luminance histogram (Rec.709) over the tile,
        /// 2) applies contrast limiting with excess redistribution,
        /// 3) converts clipped histogram to equalization LUT via CDF normalization.
        /// </summary>
        /// <param name="basePtr">Pointer to the first byte of the bitmap data.</param>
        /// <param name="stride">Stride in bytes of the bitmap row (may be &gt; width*4).</param>
        /// <param name="x0">Inclusive left bound of the tile.</param>
        /// <param name="x1">Exclusive right bound of the tile.</param>
        /// <param name="y0">Inclusive top bound of the tile.</param>
        /// <param name="y1">Exclusive bottom bound of the tile.</param>
        /// <param name="clipLimitFactor">
        /// Relative clip limit factor. Absolute limit per bin is
        /// <c>clip = clipLimitFactor * (tileArea / Bins)</c>. Use 0 to disable clipping.
        /// </param>
        /// <returns>Byte LUT of length 256 mapping input luminance → equalized luminance.</returns>
        private unsafe static byte[] BuildTileLUT(
            byte* basePtr, int stride,
            int x0, int x1, int y0, int y1,
            float clipLimitFactor)
        {
            int width = x1 - x0;
            int height = y1 - y0;
            int area = Math.Max(1, width * height);

            // 1) Luminance histogram (Rec.709 coefficients)
            int[] hist = new int[Bins];

            for (int y = y0; y < y1; y++)
            {
                byte* row = basePtr + y * stride;
                for (int x = x0; x < x1; x++)
                {
                    byte* px = row + (x << 2); // BGRA
                    float r = px[2], g = px[1], b = px[0];
                    float Y = 0.2126f * r + 0.7152f * g + 0.0722f * b;
                    int bin = (int)(Y + 0.5f);
                    if (bin < 0) bin = 0; else if (bin > 255) bin = 255;
                    hist[bin]++;
                }
            }

            // 2) Contrast limiting (clip) + redistribution of excess counts
            if (clipLimitFactor > 0f)
            {
                int clip = (int)(clipLimitFactor * (area / (float)Bins));
                if (clip < 1) clip = 1;

                int excess = 0;
                for (int i = 0; i < Bins; i++)
                {
                    if (hist[i] > clip)
                    {
                        excess += hist[i] - clip;
                        hist[i] = clip;
                    }
                }

                // Uniformly redistribute the excess into all bins (classic CLAHE)
                if (excess > 0)
                {
                    int increment = excess / Bins;
                    int remainder = excess - increment * Bins;

                    for (int i = 0; i < Bins; i++) hist[i] += increment;
                    // Distribute the leftover 1-counts into the first 'remainder' bins
                    for (int i = 0; i < remainder; i++) hist[i]++;
                }
            }

            // 3) CDF → LUT (normalize to [0..255], skipping leading empty mass via cdfMin)
            byte[] lut = new byte[Bins];
            int cdf = 0;
            int cdfMin = -1;

            for (int i = 0; i < Bins; i++)
            {
                cdf += hist[i];
                if (cdfMin < 0 && hist[i] != 0) cdfMin = cdf;
            }
            if (cdfMin < 0) cdfMin = 0; // flat tile: keep LUT at zeros

            int denom = Math.Max(1, area - cdfMin);
            int cum = 0;

            for (int i = 0; i < Bins; i++)
            {
                cum += hist[i];
                int val = (int)((cum - cdfMin) * 255L / denom);
                if (val < 0) val = 0; else if (val > 255) val = 255;
                lut[i] = (byte)val;
            }
            return lut;
        }
        #endregion
    }
}
