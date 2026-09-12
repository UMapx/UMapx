using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the color transfer filter.
    /// </summary>
    [Serializable]
    public class ColorTransfer : IBitmapFilter2
    {
        #region Private data
        private float factor;
        private bool inverted;
        private Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color transfer filter.
        /// </summary>
        /// <param name="factor">Factor [0, 10].</param>
        /// <param name="inverted">Inverted or not.</param>
        /// <param name="space">Color space.</param>
        public ColorTransfer(float factor = 0.0f, bool inverted = false, Space space = Space.RGB)
        {
            Factor = factor;
            Space = space;
            Inverted = inverted;
        }
        /// <summary>
        /// Gets or sets the contrast factor in [0, 10].
        /// The direct gain is multiplied by 1 + Factor; the inverted gain is divided by it.
        /// </summary>
        public float Factor
        {
            get
            {
                return this.factor;
            }
            set
            {
                this.factor = value;
            }
        }
        /// <summary>
        /// Gets or sets whether to use the reciprocal contrast gain.
        /// When false, the gain is target deviation divided by source deviation.
        /// When true, it is source deviation divided by target deviation.
        /// Both modes shift the target channel mean to the source channel mean.
        /// </summary>
        public bool Inverted
        {
            get
            {
                return this.inverted;
            }
            set
            {
                this.inverted = value;
            }
        }
        /// <summary>
        /// Gets or sets the color space.
        /// </summary>
        public Space Space
        {
            get
            {
                return this.space;
            }
            set
            {
                this.space = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (bmData.PixelFormat != PixelFormat.Format32bppArgb || bmSrc.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            // filter
            switch (space)
            {
                case Space.HSB:
                    ApplyHSB(bmData, bmSrc);
                    break;
                case Space.HSL:
                    ApplyHSL(bmData, bmSrc);
                    break;
                case Space.YCbCr:
                    ApplyYCbCr(bmData, bmSrc);
                    break;
                case Space.RGB:
                    ApplyRGB(bmData, bmSrc);
                    break;
                default:
                    throw new NotSupportedException("Grayscale space is not supported for this filter");
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap.</param>
        /// <param name="Src">Bitmap.</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapFormat.Lock32bpp(Data);
            BitmapData bmSrc = BitmapFormat.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapFormat.Unlock(Data, bmData);
            BitmapFormat.Unlock(Src, bmSrc);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        private unsafe void ApplyRGB(BitmapData bmData, BitmapData bmSrc)
        {
            var target = BitmapMatrix.ToRGB(bmData);
            var source = BitmapMatrix.ToRGB(bmSrc);
            Reinhard(target, source, this.factor, this.inverted);
            BitmapMatrix.FromRGB(target, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        private unsafe void ApplyHSB(BitmapData bmData, BitmapData bmSrc)
        {
            var target = BitmapMatrix.ToHSB(bmData);
            var source = BitmapMatrix.ToHSB(bmSrc);
            Reinhard(target, source, this.factor, this.inverted);
            BitmapMatrix.FromHSB(target, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        private unsafe void ApplyHSL(BitmapData bmData, BitmapData bmSrc)
        {
            var target = BitmapMatrix.ToHSL(bmData);
            var source = BitmapMatrix.ToHSL(bmSrc);
            Reinhard(target, source, this.factor, this.inverted);
            BitmapMatrix.FromHSL(target, bmData);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data.</param>
        /// <param name="bmSrc">Bitmap data.</param>
        private unsafe void ApplyYCbCr(BitmapData bmData, BitmapData bmSrc)
        {
            var target = BitmapMatrix.ToYCbCr(bmData);
            var source = BitmapMatrix.ToYCbCr(bmSrc);
            Reinhard(target, source, this.factor, this.inverted);
            BitmapMatrix.FromYCbCr(target, bmData);
        }
        #endregion

        #region Specials
        /// <summary>
        /// Transfers channel means using the existing direct or reciprocal contrast gain.
        /// </summary>
        /// <param name="target">Destination color planes, modified in place.</param>
        /// <param name="source">Reference color planes, which may have a different size.</param>
        /// <param name="factor">Contrast factor in [0, 10].</param>
        /// <param name="inverted">Whether to use the reciprocal deviation ratio and factor.</param>
        /// <remarks>Statistics use every pixel with equal weight and population variance.
        /// If either channel is constant, the destination becomes the source mean:
        /// this defines the gain even when a deviation ratio is undefined. Alpha is not processed.</remarks>
        private static void Reinhard(float[][,] target, float[][,] source, float factor = 1.0f, bool inverted = false)
        {
            for (int i = 0; i < 3; i++)
            {
                ChannelStatistics(source[i], out double sourceMean, out double sourceStd);
                ChannelStatistics(target[i], out double targetMean, out double targetStd);

                // A constant channel must never enter a 0/0 or infinite-gain calculation.
                double gain = 0;
                if (sourceStd > 0 && targetStd > 0)
                    gain = inverted ? sourceStd / targetStd / (1.0 + factor)
                        : targetStd / sourceStd * (1.0 + factor);

                var channel = target[i];
                for (int y = 0; y < channel.GetLength(0); y++)
                    for (int x = 0; x < channel.GetLength(1); x++)
                        channel[y, x] = (float)(sourceMean + (channel[y, x] - targetMean) * gain);
            }
        }

        /// <summary>
        /// Computes the mean and population standard deviation over an entire color plane.
        /// </summary>
        /// <param name="channel">Finite channel samples, each representing one equally weighted pixel.</param>
        /// <param name="mean">The arithmetic mean, or zero for an empty plane.</param>
        /// <param name="standardDeviation">The nonnegative population deviation, or zero for an empty plane.</param>
        /// <remarks>Welford's centered update in double precision avoids cancellation
        /// in nearly constant channels. Taking a deviation of per-column deviations
        /// would instead measure the variation of column contrasts, not pixel contrast.</remarks>
        private static void ChannelStatistics(float[,] channel, out double mean, out double standardDeviation)
        {
            mean = 0;
            double squaredDeviations = 0;
            long count = 0;
            foreach (float value in channel)
            {
                double delta = value - mean;
                mean += delta / ++count;
                squaredDeviations += delta * (value - mean);
            }
            standardDeviation = count == 0 ? 0 : Math.Sqrt(Math.Max(0, squaredDeviations / count));
        }
        #endregion
    }
}
