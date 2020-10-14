using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the shadows and lights correction filter.
    /// <remarks>
    /// Shadow-Highlights correction is used to correct unevenly lit images. Unlike other local algorithms
    /// (for example, Single Scale Retinex, Homomorphic Enhancement, Flat-Field Correction) filter allows you to adjust the brightness values separately in dark and bright areas
    /// Images.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ShadowsHighlightsCorrection : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double shadows;
        private double lights;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shadows and lights correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="shadows">Shadows [0, 1]</param>
        /// <param name="highlights">Highlights [0, 1]</param>
        public ShadowsHighlightsCorrection(int radius, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(radius);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Initializes the shadows and lights correction filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="shadows">Shadows [0, 1]</param>
        /// <param name="highlights">Highlights [0, 1]</param>
        public ShadowsHighlightsCorrection(int width, int height, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(width, height);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Initializes the shadows and lights correction filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="shadows">Shadows [0, 1]</param>
        /// <param name="highlights">Highlights [0, 1]</param>
        public ShadowsHighlightsCorrection(SizeInt size, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(size);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Gets or sets the shadows value [0, 1].
        /// </summary>
        public double Shadows
        {
            get
            {
                return this.shadows;
            }
            set
            {
                this.shadows = Maths.Range(value, 0, 1);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the highlights value [0, 1].
        /// </summary>
        public double Highlights
        {
            get
            {
                return this.lights;
            }
            set
            {
                this.lights = Maths.Range(value, 0, 1);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            double s = (Intensity.log05 / Math.Log(0.5 - value2gamma(shadows)));
            double l = (Intensity.log05 / Math.Log(0.5 + value2gamma(lights)));
            this.values = Intensity.LogStretch(s, l, 256);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        private double value2gamma(double v)
        {
            return (v - Intensity.logEpsilon) / 2.0;
        }
        #endregion
    }
}
