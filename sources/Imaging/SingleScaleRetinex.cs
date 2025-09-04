using System;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the Single Scale Retinex filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://dragon.larc.nasa.gov/background/pubabs/papers/gspx1.pdf
    /// </remarks>
    /// </summary>
    [Serializable]
    public class SingleScaleRetinex : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private float a;
        private float b;
        private float nbase;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Single Scale Retinex filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="nbase">Logarithm base</param>
        public SingleScaleRetinex(int radius, Space space, float a = 1, float b = 0, float nbase = Maths.Pi)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Initializes the Single Scale Retinex filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="nbase">Logarithm base</param>
        public SingleScaleRetinex(int width, int height, Space space, float a = 1, float b = 0, float nbase = Maths.Pi)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Initializes the Single Scale Retinex filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="nbase">Logarithm base</param>
        public SingleScaleRetinex(SizeInt size, Space space, float a = 1, float b = 0, float nbase = Maths.Pi)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Gets or sets the contrast [-1, 1].
        /// </summary>
        public float A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the offset value (0, 1].
        /// </summary>
        public float B
        {
            get
            {
                return this.b;
            }
            set
            {
                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the logarithm base.
        /// </summary>
        public float Base
        {
            get
            {
                return this.nbase;
            }
            set
            {
                if (value <= 0)
                    throw new ArgumentException("Logarithm base should be greater than 0");

                this.nbase = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.SingleScaleRetinex(this.nbase, this.a, this.b, 256);
        }
        #endregion
    }
}
