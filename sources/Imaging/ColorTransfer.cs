using System;
using SkiaDrawing;
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
        /// <summary>
        /// Factor.
        /// </summary>
        protected float factor;
        /// <summary>
        /// Inverted or not.
        /// </summary>
        protected bool inverted;
        /// <summary>
        /// Color space.
        /// </summary>
        protected Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color transfer filter.
        /// </summary>
        /// <param name="factor">Factor [0, 10]</param>
        /// <param name="inverted">Inverted or not</param>
        /// <param name="space">Color space</param>
        public ColorTransfer(float factor = 0.0f, bool inverted = false, Space space = Space.RGB)
        {
            Factor = factor;
            Space = space;
            Inverted = inverted;
        }
        /// <summary>
        /// Gets or sets factor.
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
        /// Inverted or not.
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData, bmSrc);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData, bmSrc);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData, bmSrc);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData, bmSrc);
                    break;
                default:
                    throw new NotSupportedException("Grayscale space is not supported for this filter");
            }
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
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
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
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
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
        /// 
        /// </summary>
        /// <param name="target"></param>
        /// <param name="source"></param>
        /// <param name="factor"></param>
        /// <param name="inverted"></param>
        private static void Reinhard(float[][,] target, float[][,] source, float factor = 1.0f, bool inverted = false)
        {
            // do not use for alpha channel
            var channels = 3;

            // do job
            for (int i = 0; i < channels; i++)
            {
                // stats
                var sourceMean = source[i].Mean().Mean();
                var sourceStd = source[i].StnDev().StnDev();
                var targetMean = target[i].Mean().Mean();
                var targetStd = target[i].StnDev().StnDev();

                // process
                var temp = target[i].Sub(targetMean);
                temp = inverted ? temp.Mul(sourceStd / targetStd * 1.0f / (1 + factor)) 
                    : temp.Mul(targetStd / sourceStd * (1 + factor));
                temp = temp.Add(sourceMean);
                target[i] = temp;
            }
        }
        #endregion
    }
}
