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
        private float shadows;
        private float lights;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shadows and lights correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="shadows">Shadows [0, 1]</param>
        /// <param name="highlights">Highlights [0, 1]</param>
        public ShadowsHighlightsCorrection(int radius, Space space, float shadows = 0.4f, float highlights = 0.4f)
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
        public ShadowsHighlightsCorrection(int width, int height, Space space, float shadows = 0.4f, float highlights = 0.4f)
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
        public ShadowsHighlightsCorrection(SizeInt size, Space space, float shadows = 0.4f, float highlights = 0.4f)
        {
            gb = new BoxBlur(size);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Gets or sets the shadows value [0, 1].
        /// </summary>
        public float Shadows
        {
            get
            {
                return this.shadows;
            }
            set
            {
                this.shadows = MathsF.Range(value, 0, 1);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the highlights value [0, 1].
        /// </summary>
        public float Highlights
        {
            get
            {
                return this.lights;
            }
            set
            {
                this.lights = MathsF.Range(value, 0, 1);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            float s = (Intensity.log05 / (float)Math.Log(0.5 - value2gamma(shadows)));
            float l = (Intensity.log05 / (float)Math.Log(0.5 + value2gamma(lights)));
            this.values = Intensity.LogStretch(s, l, 256);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Converts an intensity value to its gamma-adjusted counterpart.
        /// </summary>
        /// <param name="v">Input intensity value</param>
        /// <returns>Gamma-corrected value</returns>
        private float value2gamma(float v)
        {
            return (v - Intensity.logEpsilon) / 2.0f;
        }
        #endregion
    }
}
