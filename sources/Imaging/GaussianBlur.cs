using System;
using System.Drawing;
using System.Drawing.Imaging;
using System.Threading.Tasks;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Gaussian blur filter.
    /// <remarks>
    /// This is a very fast pseudo-Gaussian blur via N separable box filters (default 3).
    /// </remarks>
    /// </summary>
    [Serializable]
    public class GaussianBlur : IBitmapFilter2, IBitmapFilter
    {
        #region Private data
        private int radiusX;     // ≈ 3*sigmaX
        private int radiusY;     // ≈ 3*sigmaY
        private float sigmaX;    // std for X (<=0 → ignore)
        private float sigmaY;    // std for Y (<=0 → ignore)
        private int boxes = 3;   // number of box passes (>=1)
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Gaussian blur filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="sigma">Sigma</param>
        /// <param name="boxes">Number of boxes</param>
        public GaussianBlur(SizeInt size, SizeFloat sigma, int boxes = 3) : this(size.Width, size.Height, sigma.Width, sigma.Height, boxes) { }
        /// <summary>
        /// Initializes the Gaussian blur filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="sigmaX">Standard deviation X (>0)</param>
        /// <param name="sigmaY">Standard deviation Y (>0)</param>
        /// <param name="boxes">Number of boxes</param>
        public GaussianBlur(int width, int height, float sigmaX = 0f, float sigmaY = 0f, int boxes = 3)
        {
            this.Size = new SizeInt(width, height);
            this.Sigma = new SizeFloat(sigmaX, sigmaY); // can be 0 → ignored
            this.Boxes = boxes;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(this.radiusX, this.radiusY);
            }
            set
            {
                this.radiusX = value.Width;
                this.radiusY = value.Height;
            }
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeFloat Sigma
        {
            get
            {
                return new SizeFloat(this.sigmaX, this.sigmaY);
            }
            set
            {
                this.sigmaX = value.Width;
                this.sigmaY = value.Height;
            }
        }
        /// <summary>
        /// Gets or sets a number of boxes.
        /// </summary>
        public int Boxes
        {
            get => this.boxes;
            set => this.boxes = Math.Max(1, value);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            using var Src = (Bitmap)Data.Clone();
            Apply(Data, Src);
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
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            var bmDst = BitmapFormat.Lock32bpp(Data);
            var bmSrc = BitmapFormat.Lock32bpp(Src);
            try
            {
                Apply(bmDst, bmSrc);
            }
            finally
            {
                BitmapFormat.Unlock(Data, bmDst);
                BitmapFormat.Unlock(Src, bmSrc);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width;
            int height = bmData.Height;
            int n = width * height;

            // Effective sigmas per axis
            float effSigmaX = (sigmaX > 0f && sigmaX < radiusX) ? sigmaX : Math.Max(0.3333f, radiusX / 3f);
            float effSigmaY = (sigmaY > 0f && sigmaY < radiusY) ? sigmaY : Math.Max(0.3333f, radiusY / 3f);

            // Compute box sizes → per-pass radii for each axis
            int[] sizesX = BoxesForGauss(effSigmaX, boxes); // odd sizes
            int[] sizesY = BoxesForGauss(effSigmaY, boxes);

            int[] radsX = new int[sizesX.Length];
            int[] radsY = new int[sizesY.Length];

            for (int i = 0; i < sizesX.Length; i++) radsX[i] = (sizesX[i] - 1) / 2;
            for (int i = 0; i < sizesY.Length; i++) radsY[i] = (sizesY[i] - 1) / 2;

            // Extract channels
            var B = new float[n];
            var G = new float[n];
            var R = new float[n];
            var A = new float[n];
            ExtractToFloats(bmSrc, B, G, R, A);

            // Temp buffer for in-place passes
            var tmp = new float[n];

            // Blur each channel with axis-specific radii per pass
            BoxBlurN(B, tmp, width, height, radsX, radsY);
            BoxBlurN(G, tmp, width, height, radsX, radsY);
            BoxBlurN(R, tmp, width, height, radsX, radsY);
            BoxBlurN(A, tmp, width, height, radsX, radsY);

            // Write back
            WriteBack(bmData, B, G, R, A);
        }
        #endregion

        #region Private voids

        /// <summary>
        /// Apply N separable box-filter passes in-place.
        /// For each pass i:
        ///   - run a horizontal box blur with radius radsX[i] into <paramref name="tmp"/>.
        ///   - run a vertical box blur with radius radsY[i] from <paramref name="tmp"/> back into <paramref name="src"/>.
        /// If a radius on a pass is 0, the source buffer is simply copied through for that direction.
        /// The final result remains in <paramref name="src"/>.
        /// </summary>
        /// <param name="src">Working buffer (and final destination) of size w*h.</param>
        /// <param name="tmp">Temporary buffer of size w*h (reused across passes).</param>
        /// <param name="w">Image width (pixels).</param>
        /// <param name="h">Image height (pixels).</param>
        /// <param name="radsX">Per-pass horizontal radii (length >= number of passes).</param>
        /// <param name="radsY">Per-pass vertical radii (length >= number of passes).</param>
        private static void BoxBlurN(float[] src, float[] tmp, int w, int h, int[] radsX, int[] radsY)
        {
            int passes = Math.Min(radsX.Length, radsY.Length);
            for (int i = 0; i < passes; i++)
            {
                int rH = radsX[i];
                int rV = radsY[i];

                // Horizontal pass: src -> tmp (or copy if radius=0)
                if (rH > 0) BoxBlurH(src, tmp, w, h, rH);
                else Array.Copy(src, tmp, src.Length);

                // Vertical pass: tmp -> src (or copy if radius=0)
                if (rV > 0) BoxBlurV(tmp, src, w, h, rV);
                else Array.Copy(tmp, src, src.Length);
            }
        }

        /// <summary>
        /// One horizontal box blur pass with edge clamping (replicate).
        /// Uses a sliding window (running sum) to achieve O(w) per row.
        /// </summary>
        /// <param name="src">Source buffer (w*h).</param>
        /// <param name="dst">Destination buffer (w*h).</param>
        /// <param name="w">Width.</param>
        /// <param name="h">Height.</param>
        /// <param name="r">Radius (window width = 2*r+1).</param>
        private static void BoxBlurH(float[] src, float[] dst, int w, int h, int r)
        {
            float div = 2 * r + 1;

            Parallel.For(0, h, y =>
            {
                int row = y * w;

                // Initialize window sum at x=0 with clamped indices.
                float sum = 0f;
                for (int ix = -r; ix <= r; ix++)
                {
                    int cx = ix; if (cx < 0) cx = 0; else if (cx >= w) cx = w - 1;
                    sum += src[row + cx];
                }
                dst[row + 0] = sum / div;

                // Slide window by 1 pixel to the right:
                //   add incoming sample at x+r, remove outgoing at x-r-1 (both clamped).
                for (int x = 1; x < w; x++)
                {
                    int addX = x + r; if (addX >= w) addX = w - 1;
                    int subX = x - r - 1; if (subX < 0) subX = 0;
                    sum += src[row + addX] - src[row + subX];
                    dst[row + x] = sum / div;
                }
            });
        }

        /// <summary>
        /// One vertical box blur pass with edge clamping (replicate).
        /// Uses a sliding window (running sum) to achieve O(h) per column.
        /// </summary>
        /// <param name="src">Source buffer (w*h).</param>
        /// <param name="dst">Destination buffer (w*h).</param>
        /// <param name="w">Width.</param>
        /// <param name="h">Height.</param>
        /// <param name="r">Radius (window height = 2*r+1).</param>
        private static void BoxBlurV(float[] src, float[] dst, int w, int h, int r)
        {
            float div = 2 * r + 1;

            Parallel.For(0, w, x =>
            {
                // Initialize window sum at y=0 with clamped indices.
                float sum = 0f;
                for (int iy = -r; iy <= r; iy++)
                {
                    int cy = iy; if (cy < 0) cy = 0; else if (cy >= h) cy = h - 1;
                    sum += src[cy * w + x];
                }
                dst[0 * w + x] = sum / div;

                // Slide window down:
                //   add incoming sample at y+r, remove outgoing at y-r-1 (both clamped).
                for (int y = 1; y < h; y++)
                {
                    int addY = y + r; if (addY >= h) addY = h - 1;
                    int subY = y - r - 1; if (subY < 0) subY = 0;
                    sum += src[addY * w + x] - src[subY * w + x];
                    dst[y * w + x] = sum / div;
                }
            });
        }

        /// <summary>
        /// Read 32bpp ARGB bitmap data into 4 float channel planes (BGRA order).
        /// Channels are stored as [0..255] floats, no normalization or gamma applied.
        /// </summary>
        /// <param name="bm">Locked 32bpp bitmap data.</param>
        /// <param name="B">Destination Blue plane (w*h).</param>
        /// <param name="G">Destination Green plane (w*h).</param>
        /// <param name="R">Destination Red plane (w*h).</param>
        /// <param name="A">Destination Alpha plane (w*h).</param>
        private static unsafe void ExtractToFloats(BitmapData bm, float[] B, float[] G, float[] R, float[] A)
        {
            int w = bm.Width, h = bm.Height, stride = bm.Stride;
            byte* p = (byte*)bm.Scan0.ToPointer();

            Parallel.For(0, h, y =>
            {
                int row = y * stride;
                int idx = y * w;
                for (int x = 0; x < w; x++, idx++)
                {
                    int k = row + (x << 2);     // 4 bytes per pixel
                    B[idx] = p[k + 0];          // BGRA layout
                    G[idx] = p[k + 1];
                    R[idx] = p[k + 2];
                    A[idx] = p[k + 3];
                }
            });
        }

        /// <summary>
        /// Write 4 float channel planes (BGRA order) back to a 32bpp ARGB bitmap.
        /// Values are clamped to [0,255] and rounded to nearest (0.5 up).
        /// No gamma correction is applied.
        /// </summary>
        /// <param name="bm">Locked 32bpp bitmap data.</param>
        /// <param name="B">Blue plane.</param>
        /// <param name="G">Green plane.</param>
        /// <param name="R">Red plane.</param>
        /// <param name="A">Alpha plane.</param>
        private static unsafe void WriteBack(BitmapData bm, float[] B, float[] G, float[] R, float[] A)
        {
            int w = bm.Width, h = bm.Height, stride = bm.Stride;
            byte* p = (byte*)bm.Scan0.ToPointer();

            static byte ToByte(float v) => v <= 0f ? (byte)0 : (v >= 255f ? (byte)255 : (byte)(v + 0.5f));

            Parallel.For(0, h, y =>
            {
                int row = y * stride;
                int idx = y * w;
                for (int x = 0; x < w; x++, idx++)
                {
                    int k = row + (x << 2);
                    p[k + 0] = ToByte(B[idx]);
                    p[k + 1] = ToByte(G[idx]);
                    p[k + 2] = ToByte(R[idx]);
                    p[k + 3] = ToByte(A[idx]);
                }
            });
        }

        /// <summary>
        /// Compute <paramref name="n"/> odd box sizes whose convolution approximates
        /// a Gaussian with standard deviation <paramref name="sigma"/>.
        /// Returns an array of sizes where the first m entries are the lower odd width (wl),
        /// and the rest are the higher odd width (wu = wl + 2). This is the classic
        /// "fast almost-Gaussian" 3-box generalization.
        /// </summary>
        /// <param name="sigma">Target standard deviation (pixels).</param>
        /// <param name="n">Number of boxes (>=1), typically 3.</param>
        /// <returns>Array of odd box widths of length <paramref name="n"/>.</returns>
        private static int[] BoxesForGauss(float sigma, int n)
        {
            n = Math.Max(1, n);

            // Ideal (possibly non-integer) width that yields the target variance when n boxes are convolved.
            double wIdeal = Math.Sqrt(12.0 * sigma * sigma / n + 1.0);

            // Nearest lower odd width and the next higher odd width.
            int wl = (int)Math.Floor(wIdeal);
            if ((wl & 1) == 0) wl--;      // force odd
            int wu = wl + 2;

            // Number of boxes that should use the lower width wl (the rest use wu).
            double mIdeal = (12.0 * sigma * sigma - n * wl * wl - 4.0 * n * wl - 3.0 * n)
                            / (-4.0 * wl - 4.0);
            int m = (int)Math.Round(mIdeal);
            if (m < 0) m = 0; if (m > n) m = n;

            var sizes = new int[n];
            for (int i = 0; i < n; i++) sizes[i] = (i < m) ? wl : wu;
            return sizes;
        }

        #endregion
    }
}
