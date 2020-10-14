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
        private double a;
        private double b;
        private double nbase;
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
        public SingleScaleRetinex(int radius, Space space, double a = 1, double b = 0, double nbase = Math.PI)
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
        public SingleScaleRetinex(int width, int height, Space space, double a = 1, double b = 0, double nbase = Math.PI)
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
        public SingleScaleRetinex(SizeInt size, Space space, double a = 1, double b = 0, double nbase = Math.PI)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Gets or sets the contrast [-1, 1].
        /// </summary>
        public double A
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
        public double B
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
        public double Base
        {
            get
            {
                return this.nbase;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Logarithm base should be greater than 0");

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
