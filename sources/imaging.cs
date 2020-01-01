// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Collections;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Drawing.Imaging;
using System.IO;
using System.Threading.Tasks;
using UMapx.Colorspace;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Imaging
{
    // **************************************************************************
    //                          IMAGE PROCESSING TOOLBOX
    //                            UMAPX.NET FRAMEWORK
    // **************************************************************************
    // Image Processing Toolbox provides a comprehensive set of reference-
    // standard algorithms and workflow apps for image processing, analysis, 
    // visualization and algorithm development. You can perform image 
    // segmentation, image enhancement, noise reduction, geometric 
    // transformations and traditional image processing techniques.
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Interfaces
    /// <summary>
    /// Defines the interface for two images bitmap filter.
    /// </summary>
    public interface IBitmapFilter2
    {
        #region Filter components
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        void Apply(BitmapData bmData, BitmapData bmSrc);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        void Apply(Bitmap Data, Bitmap Src);
        #endregion
    }
    /// <summary>
    /// Defines the interface for bitmap filter.
    /// </summary>
    public interface IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        void Apply(BitmapData bmData);
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        void Apply(Bitmap Data);
        #endregion
    }
    /// <summary>
    /// Defines the interface of canvas.
    /// </summary>
    public interface ICanvas
    {
        #region Interface
        /// <summary>
        /// Gets or sets the width of the canvas.
        /// </summary>
        int Width { get; set; }
        /// <summary>
        /// Gets or sets the height of the canvas.
        /// </summary>
        int Height { get; set; }
        /// <summary>
        /// Creates canvas.
        /// </summary>
        /// <returns>Bitmap</returns>
        Bitmap Create();
        #endregion
    }
    #endregion

    #region Abstract classes
    /// <summary>
    /// Defines an abstract data rebuilding class.
    /// </summary>
    public abstract class Rebuilder
    {
        #region Protected components
        /// <summary>
        /// Use data rebuilding or not.
        /// </summary>
        protected bool rebuild = false;
        /// <summary>
        /// Implements the rebuilding of class data.
        /// </summary>
        protected abstract void Rebuild();
        #endregion
    }
    #endregion

    #region Rebuilder filters
    /// <summary>
    /// Defines the mask correction filter.
    /// </summary>
    public class Correction : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Values.
        /// </summary>
        protected double[] values;
        /// <summary>
        /// Color space.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the mask correction filter.
        /// </summary>
        /// <param name="values">Mask array</param>
        /// <param name="space">Color space</param>
        public Correction(double[] values, Space space)
        {
            Values = values;
            Space = space;
        }
        /// <summary>
        /// Initializes the mask correction filter.
        /// </summary>
        public Correction()
        {
            Values = new double[256];
            Space = Imaging.Space.RGB;
        }
        /// <summary>
        /// Gets or sets the mask array.
        /// </summary>
        public double[] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.Length != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
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
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData);
                    break;
                default:
                    ApplyGrayscale(bmData);
                    break;
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyRGB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    rgb = new RGB(p[k + 2], p[k + 1], p[k]);

                    rgb.Red = Maths.Byte(values[rgb.Red] * length);
                    rgb.Green = Maths.Byte(values[rgb.Green] * length);
                    rgb.Blue = Maths.Byte(values[rgb.Blue] * length);

                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSL(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                HSL hsl; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k]);
                    hsl.Lightness = values[(int)(hsl.Lightness * length)];
                    rgb = hsl.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                HSB hsb; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    hsb = HSB.FromRGB(p[k + 2], p[k + 1], p[k]);
                    hsb.Brightness = values[(int)(hsb.Brightness * length)];
                    rgb = hsb.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyYCbCr(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                YCbCr ycbcr; RGB rgb;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    ycbcr = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k]);
                    ycbcr.Y = values[(int)(ycbcr.Y * length)];
                    rgb = ycbcr.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyGrayscale(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.Length - 1;

            Parallel.For(0, height, y =>
            {
                int x, ystride, k;
                int luma;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    luma = RGB.Average(p[k + 2], p[k + 1], p[k]);
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(values[luma] * length);
                }
            }
            );

            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the RGBA mask correction filter.
    /// </summary>
    public class RGBACorrection : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Channel of RGBA.
        /// </summary>
        protected RGBA channel;
        /// <summary>
        /// Values.
        /// </summary>
        protected double[] values;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the RGBA mask correction filter.
        /// </summary>
        /// <param name="values">Mask array</param>
        /// <param name="channel">Channel</param>
        public RGBACorrection(double[] values, RGBA channel)
        {
            Values = values;
            Channel = channel;
        }
        /// <summary>
        /// Initializes the RGBA mask correction filter.
        /// </summary>
        public RGBACorrection()
        {
            Values = new double[256];
        }
        /// <summary>
        /// Gets or sets the mask array.
        /// </summary>
        public double[] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.Length != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
            }
        }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            // applying color filter:
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            int length = values.Length - 1;
            int c = (int)this.channel;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[c] = Maths.Byte(values[p[c]] * length);
                }
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the local mask correction filter.
    /// </summary>
    public class LocalCorrection : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Box blur filter.
        /// </summary>
        protected BoxBlur gb;
        /// <summary>
        /// Contrast.
        /// </summary>
        protected double[,] values;
        /// <summary>
        /// Color space.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        public LocalCorrection()
        {
            gb = new BoxBlur(3);
            Space = Space.RGB;
            Values = new double[256, 256];
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(int radius, double[,] values, Space space)
        {
            gb = new BoxBlur(radius);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(int width, int height, double[,] values, Space space)
        {
            gb = new BoxBlur(width, height);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Initializes the local mask correction filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="values">Matrix mask</param>
        public LocalCorrection(SizeInt size, double[,] values, Space space)
        {
            gb = new BoxBlur(size);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return gb.Size;
            }
            set
            {
                gb.Size = value;
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
        /// Gets or sets the matrix mask.
        /// </summary>
        public double[,] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.GetLength(0) != 256 || value.GetLength(1) != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            // box blur
            gb.Apply(bmSrc);

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
                    ApplyGrayscale(bmData, bmSrc);
                    break;
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
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
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    p[k2] = Maths.Byte(this.values[p[k2], pSrc[k2]] * length);
                    p[k1] = Maths.Byte(this.values[p[k1], pSrc[k1]] * length);
                    p[k] = Maths.Byte(this.values[p[k], pSrc[k]] * length);
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHSL(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                HSL lumI; HSL lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = HSL.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = HSL.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Lightness = this.values[(int)(lumI.Lightness * length), (int)(lumIx.Lightness * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHSB(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                HSB lumI; HSB lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = HSB.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = HSB.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Brightness = this.values[(int)(lumI.Brightness * length), (int)(lumIx.Brightness * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyYCbCr(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                YCbCr lumI; YCbCr lumIx; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    lumI = YCbCr.FromRGB(p[k2], p[k1], p[k]);
                    lumIx = YCbCr.FromRGB(pSrc[k2], pSrc[k1], pSrc[k]);
                    lumI.Y = this.values[(int)(lumI.Y * length), (int)(lumIx.Y * length)];
                    rgb = lumI.ToRGB;

                    p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyGrayscale(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;

            Parallel.For(0, height, j =>
            {
                int i, k, lumax, lumay, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    lumax = RGB.Average(p[k + 2], p[k + 1], p[k]);
                    lumay = RGB.Average(pSrc[k + 2], pSrc[k + 1], pSrc[k]);

                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(this.values[lumax, lumay] * length);
                }
            });

            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the local RGBA mask correction filter.
    /// </summary>
    public class RGBALocalCorrection : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Box blur filter.
        /// </summary>
        private BoxBlur gb;
        /// <summary>
        /// Contrast.
        /// </summary>
        protected double[,] values;
        /// <summary>
        /// Channel.
        /// </summary>
        protected RGBA channel;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local RGBA mask correction filter.
        /// </summary>
        public RGBALocalCorrection()
        {
            gb = new BoxBlur(3);
            Channel = RGBA.Red;
            Values = new double[256, 256];
        }
        /// <summary>
        /// Initializes the local RGBA mask correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="channel">Channel</param>
        /// <param name="values">Matrix mask</param>
        public RGBALocalCorrection(int radius, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(radius);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Initializes the local RGBA mask correction filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="channel">Channel</param>
        /// <param name="values">Matrix mask</param>
        public RGBALocalCorrection(int width, int height, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(width, height);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Initializes the local RGBA mask correction filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="channel">Channel</param>
        /// <param name="values">Matrix mask</param>
        public RGBALocalCorrection(SizeInt size, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(size);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return gb.Size;
            }
            set
            {
                gb.Size = value;
            }
        }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Gets or sets the matrix mask.
        /// </summary>
        public double[,] Values
        {
            get
            {
                return this.values;
            }
            set
            {
                if (value.GetLength(0) != 256 || value.GetLength(1) != 256)
                    throw new Exception("Mask size should be  256");

                this.values = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            // applying box blur:
            gb.Apply(bmSrc);

            // applying color filter:
            byte* p = (byte*)bmData.Scan0.ToPointer(), pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double length = values.GetLength(0) - 1;
            int c = (int)this.channel;

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    p[k + c] = Maths.Byte(this.values[p[k + c], pSrc[k + c]] * length);
                }
            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the point addition filter.
    /// </summary>
    public class PointAddition : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Point matrix.
        /// </summary>
        protected PointInt[,] points = new PointInt[0, 0];
        /// <summary>
        /// Image width.
        /// </summary>
        protected int width;
        /// <summary>
        /// Image height.
        /// </summary>
        protected int height;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the point addition filter.
        /// </summary>
        /// <param name="points">Array of ordered pairs of X and Y</param>
        public PointAddition(PointInt[,] points)
        {
            Points = points;
        }
        /// <summary>
        /// Initializes the point addition filter.
        /// </summary>
        public PointAddition() { }
        /// <summary>
        /// Gets or sets point matrix.
        /// </summary>
        public PointInt[,] Points
        {
            get
            {
                return this.points;
            }
            set
            {
                this.points = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get image sizes
            this.width = bmData.Width;
            this.height = bmData.Height;

            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            // do job
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                int x, y, i, c, x1, y1, jstride, k;
                jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    x = points[i, j].X;
                    y = points[i, j].Y;

                    x1 = i + x; y1 = j + y;

                    if (y1 >= 0 && y1 < height && x1 >= 0 && x1 < width)
                    {
                        c = (j + y) * stride + (i + x) * 4;

                        p[k] = pSrc[c]; // Blue
                        p[k + 1] = pSrc[c + 1]; // Green
                        p[k + 2] = pSrc[c + 2]; // Red
                        p[k + 3] = pSrc[c + 3]; // Alpha
                    }
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the point multiplication filter.
    /// </summary>
    public class PointMultiplication : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Point matrix.
        /// </summary>
        protected PointInt[,] points = new PointInt[0, 0];
        /// <summary>
        /// Image width.
        /// </summary>
        protected int width;
        /// <summary>
        /// Image height.
        /// </summary>
        protected int height;
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the point multiplication filter.
        /// </summary>
        /// <param name="points">Array of ordered pairs of X and Y</param>
        public PointMultiplication(PointInt[,] points)
        {
            Points = points;
        }
        /// <summary>
        /// Initializes the point multiplication filter.
        /// </summary>
        public PointMultiplication() { }
        /// <summary>
        /// Gets or sets point matrix.
        /// </summary>
        public PointInt[,] Points
        {
            get
            {
                return this.points;
            }
            set
            {
                this.points = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get image sizes
            this.width = bmData.Width;
            this.height = bmData.Height;

            // rebuild?
            if (rebuild == true)
                this.Rebuild(); this.rebuild = false;

            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                int x, y, i, c, jstride, k;
                jstride = j * stride;
                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    x = points[i, j].X;
                    y = points[i, j].Y;

                    if (y >= 0 && y < height && x >= 0 && x < width)
                    {
                        c = (y * stride) + (x * 4);

                        p[k] = pSrc[c]; // Blue
                        p[k + 1] = pSrc[c + 1]; // Green
                        p[k + 2] = pSrc[c + 2]; // Red
                        p[k + 3] = pSrc[c + 3]; // Alpha
                    }
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    #endregion

    #region Intensity adjustments
    /// <summary>
    /// Defines the channels inversion filter.
    /// </summary>
    public class InvertChannels : Correction, IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the channels inversion filter.
        /// </summary>
        public InvertChannels(Space space)
        {
            this.values = Intensity.Invert(256);
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild() { }
        #endregion
    }
    /// <summary>
    /// Defines the threshold filter.
    /// </summary>
    public class Threshold : Correction, IBitmapFilter
    {
        #region Private data
        private double threshold;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <param name="space">Color space</param>
        public Threshold(double threshold, Space space)
        {
            this.Value = threshold;
            this.Space = space;
        }
        /// <summary>
        /// Initializes the threshold filter.
        /// </summary>
        public Threshold()
        {
            Value = 0.5;
        }
        /// <summary>
        /// Gets or sets the threshold value [0, 1].
        /// </summary>
        public double Value
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bin(this.threshold, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the levels correction filter.
    /// <remarks>
    /// Filter usage example:
    /// https://digital-photography-school.com/using-levels-photoshop-image-correct-color-contrast/
    /// </remarks>
    /// </summary>
    public class LevelsCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble input;
        private RangeDouble output;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the levels correction filter.
        /// </summary>
        /// <param name="input">Input channel values</param>
        /// <param name="output">Output channel values</param>
        /// <param name="space">Color space</param>
        public LevelsCorrection(RangeDouble input, RangeDouble output, Space space)
        {
            Input = input; Output = output; this.Space = space;
        }
        /// <summary>
        /// Initializes the levels correction filter.
        /// </summary>
        public LevelsCorrection()
        {
            Input = new RangeDouble(0, 1);
            Output = new RangeDouble(0, 1);
        }
        /// <summary>
        /// Gets or sets input channel values.
        /// </summary>
        public RangeDouble Input
        {
            get
            {
                return this.input;
            }
            set
            {
                this.input = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets output channel values.
        /// </summary>
        public RangeDouble Output
        {
            get
            {
                return this.output;
            }
            set
            {
                this.output = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Levels(input, output, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the exposure correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Exposure_(photography)
    /// </remarks>
    /// </summary>
    public class ExposureCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double average;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the exposure correction filter.
        /// </summary>
        /// <param name="average">Average [0, 2500]</param>
        /// <param name="space">Color space</param>
        public ExposureCorrection(double average, Space space)
        {
            Average = average; this.Space = space;
        }
        /// <summary>
        /// Initializes the exposure correction filter.
        /// </summary>
        public ExposureCorrection()
        {
            Average = 128;
        }
        /// <summary>
        /// Gets or sets the average [0, 2500].
        /// </summary>
        public double Average
        {
            get
            {
                return this.average;
            }
            set
            {
                this.average = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Exposure(this.average, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the linear correction filter.
    /// </summary>
    public class LinearCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble range;
        private double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        /// <param name="range">Range values</param>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public LinearCorrection(RangeDouble range, double delta, Space space)
        {
            Range = range; Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        /// <param name="delta">Delta [-100, 100]</param>
        /// <param name="space">Color space</param>
        public LinearCorrection(double delta, Space space)
        {
            Range = new RangeDouble(0, 1); Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the linear correction filter.
        /// </summary>
        public LinearCorrection()
        {
            Range = new RangeDouble(0, 1); Delta = 0.5; this.Space = space;
        }
        /// <summary>
        /// Gets or sets range values.
        /// </summary>
        public RangeDouble Range
        {
            get
            {
                return this.range;
            }
            set
            {
                this.range = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public double Delta
        {
            get
            {
                return this.delta;
            }
            set
            {
                this.delta = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Linear(range, delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the sine correction filter.
    /// </summary>
    public class SinCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the sine correction filter.
        /// </summary>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public SinCorrection(double delta, Space space)
        {
            Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the sine correction filter.
        /// </summary>
        public SinCorrection()
        {
            Delta = 20;
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public double Delta
        {
            get
            {
                return this.delta;
            }
            set
            {
                this.delta = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Sin(delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the cosine correction filter.
    /// </summary>
    public class CosCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the cosine correction filter.
        /// </summary>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public CosCorrection(double delta, Space space)
        {
            Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the cosine correction filter.
        /// </summary>
        public CosCorrection()
        {
            Delta = 0.5;
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public double Delta
        {
            get
            {
                return this.delta;
            }
            set
            {
                this.delta = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Cos(delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the logarithmic correction filter.
    /// </summary>
    public class LogCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double nbase = 3.14f;
        private double delta = 0.5;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the logarithmic correction filter.
        /// </summary>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="delta">Delta [-1, 1]</param>
        /// <param name="space">Color space</param>
        public LogCorrection(double nbase, double delta, Space space)
        {
            Base = nbase; Delta = delta; this.Space = space;
        }
        /// <summary>
        /// Initializes the logarithmic correction filter.
        /// </summary>
        public LogCorrection() { }
        /// <summary>
        /// Gets or sets the base value of the logarithm.
        /// </summary>
        public double Base
        {
            get
            {
                return this.nbase;
            }
            set
            {
                this.nbase = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the delta value [-1, 1].
        /// </summary>
        public double Delta
        {
            get
            {
                return this.delta;
            }
            set
            {
                this.delta = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Log(nbase, delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the brightness correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    /// </summary>
    public class BrightnessCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double brightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the brightness correction filter.
        /// </summary>
        /// <param name="brightness">Brightness [-1, 1]</param>
        /// <param name="space">Color space</param>
        public BrightnessCorrection(double brightness, Space space)
        {
            Brightness = brightness; this.Space = space;
        }
        /// <summary>
        /// Initializes the brightness correction filter.
        /// </summary>
        public BrightnessCorrection()
        {
            Brightness = 0.5;
        }
        /// <summary>
        /// Gets or sets the brightness value [-1, 1].
        /// </summary>
        public double Brightness
        {
            get
            {
                return this.brightness;
            }
            set
            {
                this.brightness = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.brightness / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the contrast correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    /// </summary>
    public class ContrastCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the contrast correction filter.
        /// </summary>
        /// <param name="value">Contrast [-1, 1]</param>
        /// <param name="space">Color space</param>
        public ContrastCorrection(double value, Space space)
        {
            Contrast = value; this.Space = space;
        }
        /// <summary>
        /// Initializes the contrast correction filter.
        /// </summary>
        public ContrastCorrection()
        {
            Contrast = 0.5;
        }
        /// <summary>
        /// Gets or sets the contrast value.
        /// </summary>
        public double Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Contrast(this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the gamma correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Gamma_correction
    /// </remarks>
    /// </summary>
    public class GammaCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double g;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the gamma correction filter.
        /// </summary>
        /// <param name="g">Gamma [0, 20]</param>
        /// <param name="space">Color space</param>
        public GammaCorrection(double g, Space space)
        {
            Gamma = g; Space = space;
        }
        /// <summary>
        /// Initializes the gamma correction filter.
        /// </summary>
        public GammaCorrection()
        {
            Gamma = 2.2;
        }
        /// <summary>
        /// Gets or sets the gamma value [0, 20].
        /// </summary>
        public double Gamma
        {
            get
            {
                return g;
            }
            set
            {
                g = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Gamma(this.g, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the shift correction filter.
    /// </summary>
    public class ShiftCorrection : Correction, IBitmapFilter
    {
        #region Private data
        private double offset;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shift correction filter.
        /// </summary>
        /// <param name="offset">Offset (-0.5, 0.5)</param>
        /// <param name="space">Color space</param>
        public ShiftCorrection(double offset, Space space)
        {
            Offset = offset; Space = space;
        }
        /// <summary>
        /// Gets or sets the offset value (-0.5, 0.5).
        /// </summary>
        public double Offset
        {
            get
            {
                return this.offset;
            }
            set
            {
                this.offset = Maths.Range(value, -0.49998, 0.49998);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Shift(this.offset, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the global contrast enhancement filter.
    /// </summary>
    public class ContrastEnhancement : Correction, IBitmapFilter
    {
        #region Private data
        private double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the global contrast enhancement filter.
        /// </summary>
        /// <param name="contrast">Contrast [-1, 1]</param>
        /// <param name="space">Color space</param>
        public ContrastEnhancement(double contrast, Space space)
        {
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the contrast coefficent value [-1, 1].
        /// </summary>
        public double Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LogContrast(1 + this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the quantization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://en.wikipedia.org/wiki/Posterization
    /// </remarks>
    /// </summary>
    public class Quantization : Correction, IBitmapFilter
    {
        #region Private data
        private int levels;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the quantization filter.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="space">Color space</param>
        public Quantization(int levels, Space space)
        {
            Levels = levels; Space = space;
        }
        /// <summary>
        /// Initializes the quantization filter.
        /// </summary>
        public Quantization()
        {
            Levels = 4;
        }
        /// <summary>
        /// Gets or sets number of levels.
        /// </summary>
        public int Levels
        {
            get
            {
                return levels;
            }
            set
            {
                this.levels = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Quantize(this.levels, 256);
        }
        #endregion
    }
    #endregion

    #region Local intensity adjustments
    /// <summary>
    /// Defines the Bradley local threshold filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf
    /// </remarks>
    /// </summary>
    public class LocalThreshold : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double difference;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(int radius, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(radius);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(int width, int height, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(width, height);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Initializes the Bradley local threshold filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="difference">Difference [0, 1]</param>
        public LocalThreshold(SizeInt size, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(size);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Gets or sets the difference [0, 1].
        /// </summary>
        public double Difference
        {
            get
            {
                return difference;
            }
            set
            {
                difference = Maths.Double(value);
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bradley(this.difference, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the local contrast enhancement filter.
    /// <remarks>
    /// This filter is also known as "Unsharp Masking."
    /// More information can be found on the website:
    /// http://www.cambridgeincolour.com/tutorials/local-contrast-enhancement.htm
    /// Filter usage example:
    /// http://www.knowhowtransfer.com/photoshop-professional-plugins/alce-local-contrast-enhancer/
    /// </remarks>
    /// </summary>
    public class LocalContrastEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(int radius, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(radius);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(int width, int height, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(width, height);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local contrast enhancement filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [-1, 1]</param>
        public LocalContrastEnhancement(SizeInt size, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(size);
            Space = space;Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the contrast value [-1, 1].
        /// </summary>
        public double Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastEnhancement(this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the local contrast inversion filter.
    /// <remarks>
    /// This filter is used to equalize the illumination of images by averaging the brightness.
    /// </remarks>
    /// </summary>
    public class LocalContrastInversion : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(int radius, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(radius);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(int width, int height, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(width, height);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Initializes the local contrast inversion filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public LocalContrastInversion(SizeInt size, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(size);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Gets or sets the contrast value (0, 1].
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
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastInversion(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the contrast enhancement filter.
    /// </summary>
    public class KsiContrastEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the contrast enhancement filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        public KsiContrastEnhancement(int radius, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the contrast enhancement filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        public KsiContrastEnhancement(int width, int height, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the contrast enhancement filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        public KsiContrastEnhancement(SizeInt size, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the contrast value [-1, 1].
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
        /// Gets or sets the offset value [-1, 1].
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
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.KsiContrastEnchancement(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the filter for homomorphic processing.
    /// <remarks>
    /// A homomorphic filter is most often used to equalize the illumination of images.
    /// It simultaneously normalizes the brightness of the image and increases the contrast.
    /// 
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Homomorphic_filtering
    /// </remarks>
    /// </summary>
    public class HomomorphicEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(int radius, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(int width, int height, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Initializes the filter for homomorphic processing.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        public HomomorphicEnhancement(SizeInt size, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Gets or sets the contrast value [-1, 1].
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
                if (value <= 0)
                    throw new Exception("Invalid argument value");

                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.HomomorphicEnhancement(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Single Scale Retinex filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://dragon.larc.nasa.gov/background/pubabs/papers/gspx1.pdf
    /// </remarks>
    /// </summary>
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
    /// <summary>
    /// Defines the shadows and lights correction filter.
    /// <remarks>
    /// Shadow-Highlights correction is used to correct unevenly lit images. Unlike other local algorithms
    /// (for example, Single Scale Retinex, Homomorphic Enhancement, Flat-Field Correction) filter allows you to adjust the brightness values separately in dark and bright areas
    /// Images.
    /// </remarks>
    /// </summary>
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
    /// <summary>
    /// Defines the flat-field correction filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://imagej.net/Image_Intensity_Processing
    /// </remarks>
    /// </summary>
    public class FlatFieldCorrection : IBitmapFilter2
    {
        #region Private data
        private BoxBlur gb;     // box blur filter,
        private double mR;      // mean of red channel,
        private double mG;      // mean of green channel,
        private double mB;      // mean of blue channel.
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the flat-field correction filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public FlatFieldCorrection(int radius = 15)
        {
            gb = new BoxBlur(radius);
        }
        /// <summary>
        /// Initializes the flat-field correction filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public FlatFieldCorrection(int width, int height)
        {
            gb = new BoxBlur(width, height);
        }
        /// <summary>
        /// Initializes the flat-field correction filter.
        /// </summary>
        /// <param name="size">Radius</param>
        public FlatFieldCorrection(SizeInt size)
        {
            gb = new BoxBlur(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return gb.Size;
            }
            set
            {
                gb.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            gb.Apply(bmSrc);
            flatfield(bmData, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void flatfield(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            this.globalmeans(bmSrc); // calculating medians.

            // the formula of transforms: I(x,y) * mean(J) / J(x, y)
            for (x = 0; x < width; x++)
            {
                for (y = 0; y < height; y++, p += 4, pSrc += 4)
                {
                    if (pSrc[2] != 0)
                    {
                        p[2] = Maths.Byte(p[2] * mR / pSrc[2]);
                    }
                    if (pSrc[1] != 0)
                    {
                        p[1] = Maths.Byte(p[1] * mG / pSrc[1]);
                    }
                    if (pSrc[0] != 0)
                    {
                        p[0] = Maths.Byte(p[0] * mB / pSrc[0]);
                    }
                }
            }

            return;
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Array</returns>
        private unsafe void globalmeans(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            double total = width * height;
            double r = 0, g = 0, b = 0;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    r += p[2];
                    g += p[1];
                    b += p[0];
                }
            }

            this.mR = r / total;
            this.mG = g / total;
            this.mB = b / total;
            return;
        }
        #endregion
    }
    #endregion

    #region Grayscale filters
    /// <summary>
    /// Defines the grayscale filter based on the HSB structure.
    /// <remarks>
    /// The filter discolors the specified part of the image.
    /// </remarks>
    /// </summary>
    public class HSBGrayscale : IBitmapFilter
    {
        #region Private data
        private int min;
        private int max;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grayscale filter based on the HSB structure.
        /// </summary>
        /// <param name="hue">Hue range [0, 359]</param>
        public HSBGrayscale(RangeInt hue)
        {
            Hue = hue;
        }
        /// <summary>
        /// Initializes the grayscale filter based on the HSB structure.
        /// </summary>
        /// <param name="min">Lower bound [0, 359]</param>
        /// <param name="max">Upper bound [0, 359]</param>
        public HSBGrayscale(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the hue range [0, 359].
        /// </summary>
        public RangeInt Hue
        {
            get
            {
                return new RangeInt(min, max);
            }
            set
            {
                this.min = value.Min;
                this.max = value.Max;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSB hsb; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // This function modifies a given image in order to keep a specific hue
                    // (given too) and to desaturate the rest of the image. This procedure 
                    // originates a image with black and white colormap, excluding the parts
                    // colored with that hue.
                    // Victor Martnez Cagigal, 23/02/2015
                    //
                    // Designed for UMapx.NET by Asiryan Valeriy, 2018.

                    k = jstride + i * 4;

                    // Convert to hsb:
                    hsb = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);

                    // Getting hue and saturation parameters:
                    double hue = hsb.Hue, saturation = hsb.Saturation;

                    // Applying filter:
                    if (min < max)
                    {
                        hsb.Saturation = (hue > min && hue < max) ? saturation : 0;
                    }
                    else
                    {
                        hsb.Saturation = ((hue > min && hue <= 360) || (hue < max && hue >= 0)) ? saturation : 0;
                    }

                    // Convert to rgb:
                    rgb = hsb.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the grayscale filter based on the HSL structure.
    /// <remarks>
    /// The filter discolors the specified part of the image.
    /// </remarks>
    /// </summary>
    public class HSLGrayscale : IBitmapFilter
    {
        #region Private data
        private int min;
        private int max;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grayscale filter based on the HSL structure.
        /// </summary>
        /// <param name="hue">Hue range [0, 359]</param>
        public HSLGrayscale(RangeInt hue)
        {
            Hue = hue;
        }
        /// <summary>
        /// Initializes the grayscale filter based on the HSL structure.
        /// </summary>
        /// <param name="min">Lower bound [0, 359]</param>
        /// <param name="max">Upper bound [0, 359]</param>
        public HSLGrayscale(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Gets or sets the hue range [0, 359].
        /// </summary>
        public RangeInt Hue
        {
            get
            {
                return new RangeInt(min, max);
            }
            set
            {
                this.min = value.Min;
                this.max = value.Max;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSL hsl; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // This function modifies a given image in order to keep a specific hue
                    // (given too) and to desaturate the rest of the image. This procedure 
                    // originates a image with black and white colormap, excluding the parts
                    // colored with that hue.
                    // Victor Martnez Cagigal, 23/02/2015
                    //
                    // Designed for UMapx.NET by Asiryan Valeriy, 2018.

                    k = jstride + i * 4;

                    // Convert to hsl:
                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);

                    // Getting hue and saturation parameters:
                    double hue = hsl.Hue, saturation = hsl.Saturation;

                    // Applying filter:
                    if (min < max)
                    {
                        hsl.Saturation = (hue > min && hue < max) ? saturation : 0;
                    }
                    else
                    {
                        hsl.Saturation = ((hue > min && hue <= 360) || (hue < max && hue >= 0)) ? saturation : 0;
                    }

                    // Convert to rgb:
                    rgb = hsl.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the grayscale filter.
    /// </summary>
    public class Grayscale : IBitmapFilter
    {
        #region Private data
        private double cr;
        private double cg;
        private double cb;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grayscale filter.
        /// </summary>
        /// <param name="cr">Red</param>
        /// <param name="cg">Green</param>
        /// <param name="cb">Blue</param>
        public Grayscale(double cr, double cg, double cb)
        {
            Cr = cr;
            Cg = cg;
            Cb = cb;
        }
        /// <summary>
        /// Initializes the grayscale filter.
        /// </summary>
        public Grayscale()
        {
            Cr = 0.333f;
            Cg = 0.333f;
            Cb = 0.333f;
        }
        /// <summary>
        /// Gets or sets the red channel value.
        /// </summary>
        public double Cr
        {
            get
            {
                return this.cr;
            }
            set
            {
                this.cr = value;
            }
        }
        /// <summary>
        /// Gets or sets the green channel value.
        /// </summary>
        public double Cg
        {
            get
            {
                return this.cg;
            }
            set
            {
                this.cg = value;
            }
        }
        /// <summary>
        /// Gets or sets the blue channel value.
        /// </summary>
        public double Cb
        {
            get
            {
                return this.cb;
            }
            set
            {
                this.cb = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[0] = p[1] = p[2] = Maths.Byte(cr * p[2] + cg * p[1] + cb * p[0]);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Initializes the grayscale filter (BT709).
        /// </summary>
        public static Grayscale BT709
        {
            get
            {
                return new Grayscale(0.212f, 0.715f, 0.072f);
            }
        }
        /// <summary>
        /// Initializes the grayscale filter (R-Y).
        /// </summary>
        public static Grayscale RY
        {
            get
            {
                return new Grayscale(0.5f, 0.419f, 0.081f);
            }
        }
        /// <summary>
        /// Initializes the grayscale filter Y.
        /// </summary>
        public static Grayscale Y
        {
            get
            {
                return new Grayscale(0.299f, 0.587f, 0.114f);
            }
        }
        #endregion

        #region IsGrayscale
        /// <summary>
        /// Checks if Bitmap is a grayscale image.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Boolean</returns>
        public static bool IsGrayscale(Bitmap b)
        {
            bool isgrayscale = false;
            if (b.PixelFormat == PixelFormat.Format8bppIndexed)
            {
                isgrayscale = true;
                ColorPalette palette = b.Palette;
                Color colour;
                for (int i = 0; i < 256; i++)
                {
                    colour = palette.Entries[i];
                    if ((colour.R != i) || (colour.G != i) || (colour.B != i))
                    {
                        isgrayscale = false;
                        break;
                    }
                }
            }
            return isgrayscale;
        }
        #endregion
    }
    #endregion

    #region Color adjustments
    /// <summary>
    /// Defines the color photo filter.
    /// </summary>
    public class PhotoFilter : IBitmapFilter
    {
        #region Private data
        private double s;
        private Color color;
        private IDoubleMesh blendf;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        /// <param name="blendf">Blend function</param>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public PhotoFilter(IDoubleMesh blendf, Color color, double strength = 0.5)
        {
            BlendFunction = blendf;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public PhotoFilter(Color color, double strength = 0.5)
        {
            BlendFunction = BlendMode.Pegtop;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Initializes the color photo filter.
        /// </summary>
        public PhotoFilter()
        {
            BlendFunction = BlendMode.Pegtop;
            Color = Color.White;
            Strength = 0.5;
        }
        /// <summary>
        /// Gets or sets the blend function.
        /// </summary>
        public IDoubleMesh BlendFunction
        {
            get
            {
                return this.blendf;
            }
            set
            {
                this.blendf = value;
            }
        }
        /// <summary>
        /// gets or sets the filter color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets filter strenght [0, 1].
        /// </summary>
        public double Strength
        {
            get
            {
                return this.s;
            }
            set
            {
                this.s = Maths.Double(value);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;
            int stride = bmData.Stride;

            // filter color
            double rf = color.R / 255.0;
            double gf = color.G / 255.0;
            double bf = color.B / 255.0;

            // do job
            Parallel.For(0, height, y =>
            {
                // bitmap color
                double lm;
                double rb;
                double gb;
                double bb;

                int x, ystride, k;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    // luma
                    lm = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0;

                    // blending
                    rb = this.blendf(lm, rf) * 255.0;
                    gb = this.blendf(lm, gf) * 255.0;
                    bb = this.blendf(lm, bf) * 255.0;

                    // recording
                    p[k + 2] = Maths.Byte(rb * s + p[k + 2] * (1.0 - s));
                    p[k + 1] = Maths.Byte(gb * s + p[k + 1] * (1.0 - s));
                    p[k    ] = Maths.Byte(bb * s + p[k    ] * (1.0 - s));
                }
            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Initializes the cold filter (82).
        /// </summary>
        public static PhotoFilter Cold82
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 181, 255));
            }
        }
        /// <summary>
        /// Initializes the cold filter LBB.
        /// </summary>
        public static PhotoFilter ColdLBB
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 93, 255));
            }
        }
        /// <summary>
        /// Initializes the hot filter (81).
        /// </summary>
        public static PhotoFilter Warm81
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(235, 177, 19));
            }
        }
        /// <summary>
        /// Initializes the hot filter LBA.
        /// </summary>
        public static PhotoFilter WarmLBA
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(250, 150, 0));
            }
        }
        /// <summary>
        /// Initializes the sepia filter.
        /// </summary>
        public static PhotoFilter Sepia
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Initializes the red filter.
        /// </summary>
        public static PhotoFilter Red
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(234, 26, 26));
            }
        }
        /// <summary>
        /// Initializes the blue filter.
        /// </summary>
        public static PhotoFilter Blue
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(29, 53, 234));
            }
        }
        /// <summary>
        /// Initializes the green filter.
        /// </summary>
        public static PhotoFilter Green
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(25, 201, 25));
            }
        }
        /// <summary>
        /// Initializes the "underwater" filter.
        /// </summary>
        public static PhotoFilter Underwater
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 194, 177));
            }
        }
        /// <summary>
        /// Initializes the purple filter.
        /// </summary>
        public static PhotoFilter Purple
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(227, 24, 227));
            }
        }
        /// <summary>
        /// Initializes the orange filter.
        /// </summary>
        public static PhotoFilter Orange
        {
            get
            {
                return new PhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Initializes the yellow filter.
        /// </summary>
        public static PhotoFilter Yellow
        {
            get
            {
                return new PhotoFilter(Color.Yellow);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the temperature correction filter.
    /// <remarks>
    /// The filter uses an approximation of the Planck curve.
    /// </remarks>
    /// </summary>
    public class TemperatureCorrection : PhotoFilter, IBitmapFilter
    {
        #region Private data
        private double temperature;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the temperature correction filter.
        /// </summary>
        /// <param name="temperature">Temperature [1E3K, 1E4K]</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public TemperatureCorrection(double temperature, double strength = 0.5)
        {
            Temperature = temperature; Strength = strength;
        }
        /// <summary>
        /// Initializes the temperature correction filter.
        /// </summary>
        public TemperatureCorrection()
        {
            Temperature = 1000; Strength = 0.5;
        }
        /// <summary>
        /// Gets or sets the temperature [1E3K, 1E4K].
        /// </summary>
        public double Temperature
        {
            get
            {
                return this.temperature;
            }
            set
            {
                this.temperature = value;
                this.Color = RGB.Temp2RGB(this.temperature);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the color filter based on the YUV structure.
    /// </summary>
    public class YUVPhotoFilter : IBitmapFilter
    {
        #region Private data
        private double s;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color filter based on the YUV structure.
        /// </summary>
        /// <param name="color">Color</param>
        /// <param name="strength">Strenght [0, 1]</param>
        public YUVPhotoFilter(Color color, double strength = 0.5)
        {
            Color = color; Strength = strength;
        }
        /// <summary>
        /// Initializes the color filter based on the YUV structure.
        /// </summary>
        public YUVPhotoFilter()
        {
            Color = Color.White; Strength = 0.5;
        }
        /// <summary>
        /// Gets or sets the filter color.
        /// </summary>
        public Color Color
        {
            get
            {
                return color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets the filter strength [0, 1].
        /// </summary>
        public double Strength
        {
            get
            {
                return this.s;
            }
            set
            {
                this.s = Maths.Double(value);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double z = 1 - s;

            Parallel.For(0, height, y =>
            {
                YUV nYUV, iYUV; RGB rgb;
                int nR, nG, nB, iR, iG, iB;
                int x, ystride, k, luma;

                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;

                    iR = p[k + 2]; iG = p[k + 1]; iB = p[k];

                    luma = RGB.HDTV(iR, iG, iB);
                    nYUV = YUV.FromRGB(luma, luma, luma);
                    iYUV = YUV.FromRGB(color.R, color.G, color.B);
                    rgb  = AddColor(nYUV, iYUV).ToRGB;
                    nR   = rgb.Red; nG = rgb.Green; nB = rgb.Blue;

                    p[k + 2] = Maths.Byte(nR * s + iR * z);
                    p[k + 1] = Maths.Byte(nG * s + iG * z);
                    p[k    ] = Maths.Byte(nB * s + iB * z);
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        /// <summary>
        /// Checks if the color is a shade of gray.
        /// </summary>
        /// <param name="color">Color</param>
        /// <returns>Boolean</returns>
        public static bool IsGrayColor(Color color)
        {
            if (color.R == color.G && color.G == color.B)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Blend two colors in YUV space.
        /// </summary>
        /// <param name="yuv1">First color</param>
        /// <param name="yuv2">Second color</param>
        /// <returns>YUV</returns>
        public static YUV AddColor(YUV yuv1, YUV yuv2)
        {
            return new YUV(yuv1.Y, yuv1.U + yuv2.U, yuv1.V + yuv2.V);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Initializes the sepia filter.
        /// </summary>
        public static YUVPhotoFilter Sepia
        {
            get
            {
                return new YUVPhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Initializes the orange filter.
        /// </summary>
        public static YUVPhotoFilter Orange
        {
            get
            {
                return new YUVPhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Initializes the yellow filter.
        /// </summary>
        public static YUVPhotoFilter Yellow
        {
            get
            {
                return new YUVPhotoFilter(Color.Yellow);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the saturation correction filter.
    /// </summary>
    public class SaturationCorrection : IBitmapFilter
    {
        #region Private data
        private double saturation;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the saturation correction filter.
        /// </summary>
        /// <param name="saturation">Saturation [-100, 100]</param>
        public SaturationCorrection(double saturation)
        {
            Saturation = saturation;
        }
        /// <summary>
        /// Initializes the saturation correction filter.
        /// </summary>
        public SaturationCorrection()
        {
            Saturation = 20;
        }
        /// <summary>
        /// Gets or sets the saturation value [-100, 100].
        /// </summary>
        public double Saturation
        {
            get
            {
                return this.saturation;
            }
            set
            {
                this.saturation = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            double s = this.saturation / 255.0;

            Parallel.For(0, height, y =>
            {
                RGB rgb;
                int x, ystride, k;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    k = ystride + x * 4;
                    rgb = RGB.Saturation(p[k + 2], p[k + 1], p[k], s);

                    p[k    ] = rgb.Red;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the color replacement filter.
    /// </summary>
    public class ColorReplace : IBitmapFilter
    {
        #region Private data
        private Color input = Color.Transparent;
        private Color output = Color.Transparent;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the color replacement filter.
        /// </summary>
        /// <param name="input">Input color</param>
        /// <param name="output">Output color</param>
        public ColorReplace(Color input, Color output)
        {
            Input = input;
            Output = output;
        }
        /// <summary>
        /// Initializes the color replacement filter.
        /// </summary>
        public ColorReplace() { }
        /// <summary>
        /// Input color.
        /// </summary>
        public Color Input
        {
            get
            {
                return this.input;
            }
            set
            {
                this.input = value;
            }
        }
        /// <summary>
        /// Output color.
        /// </summary>
        public Color Output
        {
            get
            {
                return this.output;
            }
            set
            {
                this.output = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte r1 = input.R, g1 = input.G, b1 = input.B, a1 = input.A;
            byte r2 = output.R, g2 = output.G, b2 = output.B, a2 = output.A;
            int y, x, width = bmData.Width, height = bmData.Height;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    if ((p[3] == a1) && (p[2] == r1) && (p[1] == g1) && (p[0] == b1))
                    {
                        p[3] = a2; p[2] = r2; p[1] = g2; p[0] = b2;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    #endregion

    #region Color reductions
    /// <summary>
    /// Defines the tone diffusion dithering filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    public class ToneDiffusionDithering : IBitmapFilter
    {
        #region Private data
        private double[,] matrix;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the tone diffusion dithering filter.
        /// </summary>
        public ToneDiffusionDithering()
        {
            this.matrix = new double[4, 4] {
			                            {  15, 143,  47, 175 },
			                            { 207,  79, 239, 111 },
			                            {  63, 191,  31, 159 },
			                            { 255, 127, 223,  95 }};
        }
        /// <summary>
        /// Initializes the tone diffusion dithering filter.
        /// </summary>
        /// <param name="matrix">Matrix</param>
        public ToneDiffusionDithering(double[,] matrix)
        {
            Matrix = matrix;
        }
        /// <summary>
        /// Gets or sets the tone diffusion dithering matrix. 
        /// </summary>
        public double[,] Matrix
        {
            get
            {
                return this.matrix;
            }
            set
            {
                this.matrix = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int rows = matrix.GetLength(0), cols = matrix.GetLength(1);
            int y, x, width = bmData.Width, height = bmData.Height;
            byte n;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    n = (byte)(matrix[(y % rows), (x % cols)]);
                    p[2] = (byte)((p[2] <= n) ? 0 : 255);
                    p[1] = (byte)((p[1] <= n) ? 0 : 255);
                    p[0] = (byte)((p[0] <= n) ? 0 : 255);
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        private static System.Random rand = new System.Random();
        /// <summary>
        /// Initializes the order dithering filter.
        /// <remarks>
        /// More information can be found on the website:
        /// http://en.wikipedia.org/wiki/Ordered_dithering
        /// Filter usage example:
        /// https://en.wikipedia.org/wiki/Dither
        /// </remarks>
        /// </summary>
        /// <param name="radius">Radius [0, 255]</param>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Order(int radius)
        {
            byte c = (byte)(256 / radius / radius + 1), d = 0;
            double[,] table = new double[radius, radius];
            int i, j;

            for (i = 0; i < radius; i++)
            {
                for (j = 0; j < radius; j++, d += c)
                {
                    table[i, j] = d;
                }
            }

            return new ToneDiffusionDithering(table);
        }
        /// <summary>
        /// Initializes the random dithering filter.
        /// </summary>
        /// <param name="radius">Radius [0, 255]</param>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Random(int radius)
        {
            double[,] table = new double[radius, radius];
            int i, j;

            for (i = 0; i < radius; i++)
            {
                for (j = 0; j < radius; j++)
                {
                    table[i, j] = (byte)rand.Next(0, 255);
                }
            }

            return new ToneDiffusionDithering(table);
        }
        /// <summary>
        /// Initializes the classic dithering filter.
        /// </summary>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Basic()
        {
            return new ToneDiffusionDithering(new double[4, 4] {
			                            {  15, 143,  47, 175 },
			                            { 207,  79, 239, 111 },
			                            {  63, 191,  31, 159 },
			                            { 255, 127, 223,  95 }});
        }
        /// <summary>
        /// Initializes the Bayer dithering filter.
        /// </summary>
        /// <returns>Tone diffusion dithering filter</returns>
        public static ToneDiffusionDithering Bayer()
        {
            return new ToneDiffusionDithering(new double[,] {
                                      {   0, 192,  48, 240 },
                                      { 128,  64, 176, 112 },
                                      {  32, 224,  16, 208 },
                                      { 160,  96, 144,  80 } });
        }
        #endregion
    }
    /// <summary>
    /// Defines the error diffusion dithering filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    public class ErrorDiffusionDithering : Rebuilder, IBitmapFilter
    {
        #region Private data
        private int x;
        private int y;
        private int width;
        private int height;
        private int stride;
        private double[][] matrix;
        private double summary;
        private int levels;
        private double[] table;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the error diffusion dithering filter.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="matrix">Matrix</param>
        public ErrorDiffusionDithering(int levels, double[][] matrix)
        {
            this.Levels = levels;
            this.Matrix = matrix;
        }
        /// <summary>
        /// Gets or sets the number of levels.
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
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets the matrix.
        /// </summary>
        public double[][] Matrix
        {
            get
            {
                return this.matrix;
            }
            set
            {
                this.matrix = value;
                this.summary = 0;
                int n = matrix.Length;
                double[] row;
                int i, j, k;

                for (i = 0; i < n; i++)
                {
                    row = matrix[i];
                    k = row.Length;

                    for (j = 0; j < k; j++)
                    {
                        summary += row[j];
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Bitmap</returns>
        public unsafe void Apply(BitmapData bmData)
        {
            // rebuild?
            if (rebuild == true)
            {
                this.Rebuild(); this.rebuild = false;
            }

            // params
            this.width = bmData.Width;
            this.height = bmData.Height;
            this.stride = bmData.Stride;
            int length = table.Length;
            byte* ptr = (byte*)bmData.Scan0.ToPointer();
            int r, g, b;
            Color color;

            // for each line
            for (y = 0; y < height; y++)
            {
                // for each pixels
                for (x = 0; x < width; x++, ptr += 4)
                {
                    // current
                    r = ptr[2];
                    g = ptr[1];
                    b = ptr[0];

                    // get color from palette, which is the closest to current pixel's value
                    color = GetColor(r, g, b, table);
                    ptr[2] = color.R;
                    ptr[1] = color.G;
                    ptr[0] = color.B;

                    // do error diffusion
                    Diffuse(r - color.R, g - color.G, b - color.B, ptr);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Bitmap</returns>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.table = Intensity.Quantize(this.levels, 256).Mul(255);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="rError"></param>
        /// <param name="gError"></param>
        /// <param name="bError"></param>
        /// <param name="ptr"></param>
        protected unsafe void Diffuse(int rError, int gError, int bError, byte* ptr)
        {
            double edR;	// error diffusion
            double edG;	// error diffusion
            double edB;	// error diffusion

            // do error diffusion to right-standing neighbors
            double[] row = matrix[0];
            int jI, jP, i, k, jC, n;
            int length = matrix.Length;


            for (jI = 1, jP = 4, jC = 0, k = row.Length; jC < k; jI++, jC++, jP += 4)
            {
                if (x + jI >= width)
                    break;

                edR = ptr[jP + 2] + (rError * row[jC]) / summary;
                ptr[jP + 2] = Maths.Byte(edR);

                edG = ptr[jP + 1] + (gError * row[jC]) / summary;
                ptr[jP + 1] = Maths.Byte(edG);

                edB = ptr[jP + 0] + (bError * row[jC]) / summary;
                ptr[jP + 0] = Maths.Byte(edB);
            }

            // do error diffusion to bottom neigbors
            for (i = 1, n = length; i < n; i++)
            {
                if (y + i >= height)
                    break;

                // move pointer to next image line
                ptr += stride;

                // get coefficients of the row
                row = matrix[i];

                // process the row
                for (jC = 0, k = row.Length, jI = -(k >> 1), jP = -(k >> 1) * 4; jC < k; jI++, jC++, jP += 4)
                {
                    if (x + jI >= width)
                        break;
                    if (x + jI < 0)
                        continue;

                    edR = ptr[jP + 2] + (rError * row[jC]) / summary;
                    ptr[jP + 2] = Maths.Byte(edR);

                    edG = ptr[jP + 1] + (gError * row[jC]) / summary;
                    ptr[jP + 1] = Maths.Byte(edG);

                    edB = ptr[jP + 0] + (bError * row[jC]) / summary;
                    ptr[jP + 0] = Maths.Byte(edB);

                }
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="red"></param>
        /// <param name="green"></param>
        /// <param name="blue"></param>
        /// <param name="table"></param>
        /// <returns></returns>
        private Color GetColor(int red, int green, int blue, double[] table)
        {
            byte r = Maths.Byte(table[red]);
            byte g = Maths.Byte(table[green]);
            byte b = Maths.Byte(table[blue]);
            return Color.FromArgb(r, g, b);
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Initializes the Atkinson dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Atkinson
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[3][] {
                new double[2]   {          1, 1 },
                new double[5]   { 0, 1, 1, 1, 0 },
                new double[5]   { 0, 0, 1, 0, 0 } });
            }
        }
        /// <summary>
        /// Initializes the Burkes dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Burkes
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                    new double[2] { 8, 4 },
                    new double[5] { 2, 4, 8, 4, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Fan dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Fan
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                    new double[1] { 8 },
                    new double[4] { 1, 1, 2, 4 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra lite dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering SierraLite
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                new double[1] { 2 },
                new double[2] { 1, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Sierra
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[3][] {
                new double[2] { 5, 3 },
                new double[5] { 2, 4, 5, 4, 2 },
                new double[3] { 2, 3, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Sierra lite dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering SierraTowsRows
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                new double[2] { 4, 3 },
                new double[5] { 1, 2, 3, 2, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Flyd-Steinberg dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering FloydSteinberg
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                new double[1] {       7 },
                new double[3] { 3, 5, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Jarvis-Judice-Ninke dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering JarvisJudiceNinke
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[3][] {
                new double[2] {          7, 5 },
                new double[5] { 3, 5, 7, 5, 3 },
                new double[5] { 1, 3, 5, 3, 1 } });
            }
        }
        /// <summary>
        /// Initializes the Stevenson dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Stevenson
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[3][] {
                new double[4] { 12, 26, 30, 16 },
                new double[3] { 12, 26, 12   },
                new double[4] { 5, 12, 12, 5 } });
            }
        }
        /// <summary>
        /// Initializes the Shiau dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Shiau
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[2][] {
                new double[1] { 4 },
                new double[3] { 1, 1, 2 } });
            }
        }
        /// <summary>
        /// Initializes the Stucki dithering filter.
        /// </summary>
        /// <returns>Error diffusion dithering filter</returns>
        public static ErrorDiffusionDithering Stucki
        {
            get
            {
                return new ErrorDiffusionDithering(2,
                    new double[3][] {
                new double[2] { 8, 4 },
                new double[5] { 2, 4, 8, 4, 2 },
                new double[5] { 1, 2, 4, 2, 1 } });
            }
        }
        #endregion
    }
    #endregion

    #region Noise filters
    /// <summary>
    /// Defines the additive noise filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Gaussian_noise
    /// </remarks>
    /// </summary>
    public class AdditiveNoise : IBitmapFilter
    {
        #region Private data
        private int amount = 10;
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the additive noise filter.
        /// </summary>
        public AdditiveNoise() { }
        /// <summary>
        /// Initializes the additive noise filter.
        /// </summary>
        /// <param name="amount">Amount [0, 100]</param>
        public AdditiveNoise(int amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Gets or sets the amout value [0, 100].
        /// </summary>
        public int Amount
        {
            get
            {
                return this.amount;
            }
            set
            {
                this.amount = Math.Max(0, Math.Min(100, value));
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            int stride = bmData.Stride;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(p[2] + generator.Next(-amount, amount));
                    p[1] = Maths.Byte(p[1] + generator.Next(-amount, amount));
                    p[0] = Maths.Byte(p[0] + generator.Next(-amount, amount));
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the salt and pepper noise filter.
    /// <remarks>
    /// Filter usage example:
    /// https://en.wikipedia.org/wiki/Salt-and-pepper_noise
    /// </remarks>
    /// </summary>
    public class SaltAndPepper : IBitmapFilter
    {
        #region Private data
        private double amount = 10;
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the salt and pepper noise filter.
        /// </summary>
        public SaltAndPepper() { }
        /// <summary>
        /// Initializes the salt and pepper noise filter.
        /// </summary>
        /// <param name="amount">Amount [0, 100].</param>
        public SaltAndPepper(double amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Gets or sets the amout value [0, 100].
        /// </summary>
        public double Amount
        {
            get
            {
                return amount;
            }
            set
            {
                amount = Math.Max(0, Math.Min(100, value));
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            int noisyPixels = (int)((width * height * amount) / 100);
            byte[] values = new byte[2] { 0, 255 };
            int stride = bmData.Stride;
            int i, colorPlane;

            for (i = 0; i < noisyPixels; i++)
            {
                x = generator.Next(width);
                y = generator.Next(height);
                colorPlane = generator.Next(3);
                p[y * stride + x * 4 + colorPlane] = values[generator.Next(2)];
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    #endregion

    #region Channels
    /// <summary>
    /// Defines the transparency correction filter.
    /// </summary>
    public class TransparencyCorrection : RGBACorrection, IBitmapFilter
    {
        #region Private data
        private double transparency;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the transparency correction filter.
        /// </summary>
        /// <param name="transparency">Transparency [-1, 1]</param>
        public TransparencyCorrection(double transparency)
        {
            this.Channel = RGBA.Alpha;
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the transparency correction filter.
        /// </summary>
        public TransparencyCorrection()
        {
            this.Channel = RGBA.Alpha;
            Transparency = 0;
        }
        /// <summary>
        /// Gets or sets the transparency value [-1, 1].
        /// </summary>
        public double Transparency
        {
            get
            {
                return this.transparency;
            }
            set
            {
                this.transparency = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.transparency, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the channel level correction filter.
    /// <remarks>
    /// Filter usage example:
    /// https://digital-photography-school.com/using-levels-photoshop-image-correct-color-contrast/
    /// </remarks>
    /// </summary>
    public class LevelsChannelCorrection : RGBACorrection, IBitmapFilter
    {
        #region Private data
        private RangeDouble input;
        private RangeDouble output;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel level correction filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        /// <param name="input">Input channel values</param>
        /// <param name="output">Output channel values</param>
        public LevelsChannelCorrection(RGBA channel, RangeDouble input, RangeDouble output)
        {
            this.Channel = channel; Input = input; Output = output;
        }
        /// <summary>
        /// Initializes the channel level correction filter.
        /// </summary>
        public LevelsChannelCorrection()
        {
            Input = new RangeDouble(0, 255);
            Output = new RangeDouble(0, 255);
        }
        /// <summary>
        /// Gets or sets input channel values.
        /// </summary>
        public RangeDouble Input
        {
            get
            {
                return this.input;
            }
            set
            {
                this.input = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Gets or sets output channel values.
        /// </summary>
        public RangeDouble Output
        {
            get
            {
                return this.output;
            }
            set
            {
                this.output = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Levels(input, output, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the channel rotation filter.
    /// </summary>
    public class RotateChannel : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the channel rotation filter.
        /// </summary>
        public RotateChannel() { }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte red, green, blue;
            int y, x, height = bmData.Height, width = bmData.Width;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    blue = p[0];
                    green = p[1];
                    red = p[2];

                    p[0] = red;
                    p[1] = blue;
                    p[2] = green;
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the channel equalization filter.
    /// </summary>
    public class EqualizeChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel equalization filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public EqualizeChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel equalization filter.
        /// </summary>
        public EqualizeChannel() { }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int c1 = (int)this.channel;

            for (int y = 0; y < bmData.Height; y++)
            {
                for (int x = 0; x < bmData.Width; x++, p += 4)
                {
                    p[2] = p[c1];
                    p[1] = p[c1];
                    p[0] = p[c1];
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the channel hide filter.
    /// </summary>
    public class HideChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel hide filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public HideChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel hide filter.
        /// </summary>
        public HideChannel() { }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int c1 = (int)this.channel;

            for (int y = 0; y < bmData.Height; y++)
            {
                for (int x = 0; x < bmData.Width; x++, p += 4)
                {
                    p[c1] = 0;
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the channel show filter.
    /// </summary>
    public class ShowChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the channel show filter.
        /// </summary>
        /// <param name="channel">Channel of RGBA model</param>
        public ShowChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Initializes the channel show filter.
        /// </summary>
        public ShowChannel() { }
        /// <summary>
        /// Gets or sets the channel of the RGBA model.
        /// </summary>
        public RGBA Channel
        {
            get
            {
                return this.channel;
            }
            set
            {
                this.channel = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            for (int y = 0; y < bmData.Height; y++)
            {
                for (int x = 0; x < bmData.Width; x++, p += 4)
                {
                    switch (channel)
                    {
                        case RGBA.Blue:
                            p[1] = 0; p[2] = 0;
                            break;

                        case RGBA.Green:
                            p[2] = 0; p[0] = 0;
                            break;

                        case RGBA.Red:
                            p[1] = 0; p[0] = 0;
                            break;

                        case RGBA.Alpha:
                            p[0] = p[1] = p[2] = 255;
                            break;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    #endregion

    #region Canvas
    /// <summary>
    /// Defines the color canvas class.
    /// </summary>
    public class CanvasColor : ICanvas
    {
        #region Private data
        private Color color = Color.White;
        private Bitmap bitmap = new Bitmap(256, 256);
        private int width = 256, height = 256;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the color canvas class.
        /// </summary>
        public CanvasColor() { }
        /// <summary>
        /// Initializes the color canvas class.
        /// </summary>
        /// <param name="width">Canvas width</param>
        /// <param name="height">Canvas height</param>
        /// <param name="color">Color</param>
        public CanvasColor(int width, int height, Color color)
        {
            Width = width; Height = height; Color = color;
        }
        /// <summary>
        /// Gets or sets canvas color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Gets or sets the width of the canvas.
        /// </summary>
        public int Width
        {
            get
            {
                return this.width;
            }
            set
            {
                this.width = value;
            }
        }
        /// <summary>
        /// Gets or sets the height of the canvas.
        /// </summary>
        public int Height
        {
            get
            {
                return this.height;
            }
            set
            {
                this.height = value;
            }
        }
        /// <summary>
        /// Creates canvas.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap Create()
        {
            bitmap = new Bitmap(width, height);
            Graphics graphics = Graphics.FromImage(bitmap);
            SolidBrush brush = new SolidBrush(color);
            graphics.FillRectangle(brush, 0, 0, width, height);
            graphics.Dispose();
            brush.Dispose();
            return bitmap;
        }
        #endregion
    }
    /// <summary>
    /// Defines the gradient canvas class.
    /// </summary>
    public class CanvasGradient : ICanvas
    {
        #region Private data
        private Color color1 = Color.White;
        private Color color2 = Color.Black;
        private Bitmap bitmap = new Bitmap(256, 256);
        private int width = 256, height = 256;
        private double angle = 0;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the gradient canvas class.
        /// </summary>
        public CanvasGradient() { }
        /// <summary>
        /// Initializes the gradient canvas class.
        /// </summary>
        /// <param name="width">Canvas width</param>
        /// <param name="height">Canvas height</param>
        /// <param name="angle">Angle</param>
        /// <param name="color1">First color</param>
        /// <param name="color2">Second color</param>
        public CanvasGradient(int width, int height, double angle, Color color1, Color color2)
        {
            Width = width; Height = height; Angle = angle; Color1 = color1; Color2 = color2;
        }
        /// <summary>
        /// Gets or sets the first color.
        /// </summary>
        public Color Color1
        {
            get
            {
                return this.color1;
            }
            set
            {
                this.color1 = value;
            }
        }
        /// <summary>
        /// Gets or sets the second color.
        /// </summary>
        public Color Color2
        {
            get
            {
                return this.color2;
            }
            set
            {
                this.color2 = value;
            }
        }
        /// <summary>
        /// Gets or sets the angle value.
        /// </summary>
        public double Angle
        {
            get
            {
                return this.angle;
            }
            set
            {
                this.angle = value;
            }
        }
        /// <summary>
        /// Gets or sets the width of the canvas.
        /// </summary>
        public int Width
        {
            get
            {
                return this.width;
            }
            set
            {
                this.width = value;
            }
        }
        /// <summary>
        /// Gets or sets the height of the canvas.
        /// </summary>
        public int Height
        {
            get
            {
                return this.height;
            }
            set
            {
                this.height = value;
            }
        }
        /// <summary>
        /// Creates canvas.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap Create()
        {
            bitmap = new Bitmap(width, height);
            Graphics graphics = Graphics.FromImage(bitmap);
            Rectangle rectangle = new Rectangle(0, 0, width, height);
            LinearGradientBrush brush = new LinearGradientBrush(rectangle, color1, color2, (float)angle);
            graphics.FillRectangle(brush, rectangle);
            graphics.Dispose();
            brush.Dispose();
            return bitmap;
        }
        #endregion
    }
    #endregion

    #region Points
    /// <summary>
    /// Defines the pixelation filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Pixelate : PointAddition, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the pixelation filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Pixelate(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the pixelation filter.
        /// </summary>
        public Pixelate() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Pixelate(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Defines the grid filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Grid : PointAddition, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the grid filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Grid(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the grid filter.
        /// </summary>
        public Grid() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Grid(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Defines the jitter filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Jitter : PointAddition, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the jitter filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Jitter(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the jitter filter.
        /// </summary>
        public Jitter() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Noise(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Defines the water filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Water : PointMultiplication, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the water filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Water(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the water filter.
        /// </summary>
        public Water() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Water(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Defines the noise filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Noise : PointMultiplication, IBitmapFilter2
    {
        #region Private data
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the noise filter.
        /// </summary>
        /// <param name="value">Value [0, 100]</param>
        public Noise(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Initializes the noise filter.
        /// </summary>
        public Noise() { }
        /// <summary>
        /// Gets or sets the value.
        /// </summary>
        public int Value
        {
            get
            {
                return this.value;
            }
            set
            {
                this.value = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Noise(this.width, this.height, this.value);
        }
        #endregion
    }
    #endregion

    #region Morphs
    /// <summary>
    /// Defines the shift filter.
    /// </summary>
    public class Shift : IBitmapFilter2
    {
        #region Private data
        private int x;
        private int y;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        /// <param name="x">Offset value of axis X</param>
        /// <param name="y">Offset value of axis Y</param>
        public Shift(int x, int y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        /// <param name="point">A pair of integers representing an ordered pair of X and Y coordinates</param>
        public Shift(PointInt point)
        {
            X = point.X;
            Y = point.Y;
        }
        /// <summary>
        /// Initializes the shift filter.
        /// </summary>
        public Shift() { }
        /// <summary>
        /// Gets or sets the offset value of axis X.
        /// </summary>
        public int X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Gets or sets the offset value of axis Y.
        /// </summary>
        public int Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // image properties:
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int stride = bmSrc.Stride;

            // exception!
            if (bmData.Width  != width ||
                bmData.Height != height)
                throw new Exception("Input image sizes must be the same");

            // applying only for X:
            if (x != 0 && y == 0)
                ShiftX(bmData, bmSrc, width, height, stride);

            // applying only for Y:
            else if (y != 0 && x == 0)
                ShiftY(bmData, bmSrc, width, height, stride);

            // applying for two sides:
            else
            {
                ShiftX(bmSrc, bmData, width, height, stride);
                ShiftY(bmData, bmSrc, width, height, stride);
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="stride">Stride</param>
        private unsafe void ShiftY(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1;
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j + y; // Offset of axis Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i; // Offset of axis X

                    if (sY < height && sY >= 0)
                    {
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front    ];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sY < 0) ? -1 : 1;

                        reverse = Math.Abs(sX * 4 + (sY - f * height) * stride);
                        p[0] = pSrc[reverse    ];
                        p[1] = pSrc[reverse + 1];
                        p[2] = pSrc[reverse + 2];
                        p[3] = pSrc[reverse + 3];
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="stride">Stride</param>
        private unsafe void ShiftX(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1;
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j; // Offset of axis Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i + x; // Offset of axis X

                    if (sX < width && sX >= 0)
                    {
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sX < 0) ? -1 : 1;
                        
                        reverse = Math.Abs((sX - f * width) * 4 + sY * stride);
                        p[0] = pSrc[reverse];
                        p[1] = pSrc[reverse + 1];
                        p[2] = pSrc[reverse + 2];
                        p[3] = pSrc[reverse + 3];
                    }
                }
            }
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the flip filter.
    /// </summary>
    public class Flip : IBitmapFilter
    {
        #region Private data
        private bool x = true;
        private bool y = true;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the flip filter.
        /// </summary>
        public Flip() { }
        /// <summary>
        /// Initializes the flip filter.
        /// </summary>
        /// <param name="x">Flip X</param>
        /// <param name="y">Flip Y</param>
        public Flip(bool x, bool y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Gets or sets flip X.
        /// </summary>
        public bool X
        {
            get
            {
                return this.x;
            }
            set
            {
                this.x = value;
            }
        }
        /// <summary>
        /// Gets or sets flip Y.
        /// </summary>
        public bool Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            #region Data
            int pixel = 4, start = 0;
            int width = bmData.Width, height = bmData.Height;
            int offset = width * pixel, stride = bmData.Stride;
            #endregion

            #region FlipY
            // Flip X:
            if (this.x == true)
            {
                byte* p = (byte*)bmData.Scan0.ToPointer();
                byte* pSrc = (byte*)bmData.Scan0.ToPointer();
                int s0 = stride - (width >> 1) * pixel;
                int s1 = stride + (width >> 1) * pixel;
                pSrc += (width - 1) * pixel;
                int l, w2;
                byte b;

                for (int k = 0; k < height; k++)
                {
                    l = 0; w2 = (width >> 1);

                    while (l < w2)
                    {
                        b = p[2];
                        p[2] = pSrc[2];
                        pSrc[2] = b;

                        b = p[1];
                        p[1] = pSrc[1];
                        pSrc[1] = b;

                        b = *p;
                        *p = *pSrc;
                        *pSrc = b;

                        l++;
                        p += 4;
                        pSrc -= 4;
                    }
                    p += s0;
                    pSrc += s1;
                }
            }
            #endregion

            #region FlipX
            // Flip Y:
            if (this.y == true)
            {
                byte* p = (byte*)bmData.Scan0.ToPointer();
                byte* pSrc = (byte*)bmData.Scan0.ToPointer(); 
                pSrc += (height - 1) * stride;
                int offset2 = stride - width * pixel;
                int m = 0, n = start;
                int h2 = (height >> 1);
                byte b2;

                while (m < h2)
                {
                    n = start;
                    while (n < offset)
                    {
                        b2 = *p;
                        *p = *pSrc;
                        *pSrc = b2;
                        n++;
                        p++;
                        pSrc++;
                    }
                    p += offset2;
                    pSrc += offset2 - stride - stride;
                    m++;
                }
            }
            #endregion
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the rotation filter.
    /// </summary>
    public class Rotate : IBitmapFilter2
    {
        #region Private data
        private double angle;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the rotation filter.
        /// </summary>
        /// <param name="angle">Angle</param>
        /// <param name="color">Background color</param>
        public Rotate(double angle, Color color)
        {
            Angle = angle; Color = color;
        }
        /// <summary>
        /// Gets or sets angle value.
        /// </summary>
        public double Angle
        {
            get
            {
                return this.angle;
            }
            set
            {
                this.angle = value;
            }
        }
        /// <summary>
        /// Gets or sets background color.
        /// </summary>
        public Color Color
        {
            get
            {
                return this.color;
            }
            set
            {
                this.color = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // exception!
            if (bmSrc.Width != bmData.Width || bmSrc.Height != bmData.Height)
                throw new Exception("Input image sizes must be the same");

            // get source image size
            int width = bmSrc.Width, height = bmSrc.Height, stride = bmSrc.Stride;
            double xradius = (width - 1) / 2.0, yradius = (height - 1) / 2.0;

            // angle's sine and cosine
            double angleRad = -angle * Math.PI / 180;
            double angleCos = Math.Cos(angleRad);
            double angleSin = Math.Sin(angleRad);

            // destination pixel's coordinate relative to image center
            double cx, cy = -yradius;

            // fill values
            byte fillA = color.A;
            byte fillR = color.R;
            byte fillG = color.G;
            byte fillB = color.B;

            // do the job
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* p;

            // source pixel's coordinates
            int ox, oy, y, x;

            for (y = 0; y < height; y++)
            {
                cx = -xradius;
                for (x = 0; x < width; x++, dst += 4)
                {
                    // coordinate of the nearest point
                    ox = (int)( angleCos * cx + angleSin * cy + xradius);
                    oy = (int)(-angleSin * cx + angleCos * cy + yradius);

                    // validate source pixel's coordinates
                    if ((ox < 0) || (oy < 0) || (ox >= width) || (oy >= height))
                    {
                        // fill destination image with filler
                        dst[3] = fillA;
                        dst[2] = fillR;
                        dst[1] = fillG;
                        dst[0] = fillB;
                    }
                    else
                    {
                        // fill destination image with pixel from source image
                        p = src + oy * stride + ox * 4;

                        dst[3] = p[3];
                        dst[2] = p[2];
                        dst[1] = p[1];
                        dst[0] = p[0];
                    }
                    cx++;
                }
                cy++;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the crop filter.
    /// </summary>
    public class Crop : IBitmapFilter2
    {
        #region Private data
        Rectangle rectangle;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the crop filter.
        /// </summary>
        /// <param name="rectangle">Rectangle</param>
        public Crop(Rectangle rectangle)
        {
            Rectangle = rectangle;
        }
        /// <summary>
        /// Initializes the crop filter.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public Crop(int x, int y, int width, int height)
        {
            Rectangle = new Rectangle(x, y, width, height);
        }
        /// <summary>
        /// Gets or sets rectangle.
        /// </summary>
        public Rectangle Rectangle
        {
            get
            {
                return rectangle;
            }
            set
            {
                rectangle = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image:
            int width = bmSrc.Width;
            int height = bmSrc.Height;

            // images strides:
            int srcStride = bmSrc.Stride;
            int dstStride = bmData.Stride;

            // do the job:
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();
            byte* dst, src, p;

            // source pixel's coordinates from rectangle:
            int startX = Maths.Range(rectangle.X, 0, width);
            int startY = Maths.Range(rectangle.Y, 0, height);
            int endX   = Maths.Range(rectangle.Width - startX, 0, width);
            int endY   = Maths.Range(rectangle.Height - startY, 0, height);

            // pixel offsets:
            int x, y, i;

            // for each line:
            for (y = 0; y < endY; y++)
            {
                dst = pDst + dstStride * y;
                src = pSrc + srcStride * (y + startX);

                // for each pixel:
                for (x = 0; x < endX; x++)
                {
                    p = src + 4 * (x + startX);

                    for (i = 0; i < 4; i++, dst++, p++)
                    {
                        *dst = *p;
                    }
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        #endregion
    }
    /// <summary>
    /// Defines the resize filter.
    /// </summary>
    public class Resize : IBitmapFilter2
    {
        #region Private data
        int newWidth;
        int newHeight;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the resize filter.
        /// </summary>
        /// <param name="width">Width</param>
        /// <param name="height">Height</param>
        public Resize(int width = 512, int height = 512)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the resize filter.
        /// </summary>
        /// <param name="size">Size</param>
        public Resize(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Gets or sets image size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(newWidth, newHeight);
            }
            set
            {
                this.newWidth = value.Width;
                this.newHeight = value.Height;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image:
            int width = Maths.Range(newWidth, 0, bmData.Width);
            int height = Maths.Range(newHeight, 0, bmData.Height);
            int srcStride = bmSrc.Stride, dstStride = bmData.Stride;
            double xFactor = (double)bmSrc.Width / newWidth;
            double yFactor = (double)bmSrc.Height / newHeight;

            // do the job:
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            byte* pDst = (byte*)bmData.Scan0.ToPointer();
            byte* dst, src, p;

            // source pixel's coordinates
            int x, y, i;

            // for each line
            for (y = 0; y < height; y++)
            {
                dst = pDst + dstStride * y;
                src = pSrc + srcStride * ((int)(y * yFactor));

                // for each pixel
                for (x = 0; x < width; x++)
                {
                    p = src + 4 * ((int)(x * xFactor));

                    for (i = 0; i < 4; i++, dst++, p++)
                    {
                        *dst = *p;
                    }
                }
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        #endregion
    }
    /// <summary>
    /// Defines the merge filter.
    /// </summary>
    public class Merge : IBitmapFilter2
    {
        #region Private data
        private int transparency = 255;
        private PointInt point;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        /// <param name="point">A pair of integers representing an ordered pair of X and Y coordinates</param>
        /// <param name="transparency">Transparency [0, 255]</param>
        public Merge(PointInt point, int transparency = 128)
        {
            Point = point;
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="transparency">Transparency [0, 255]</param>
        public Merge(int x, int y, int transparency)
        {
            Point = new PointInt(x, y);
            Transparency = transparency;
        }
        /// <summary>
        /// Initializes the merge filter.
        /// </summary>
        public Merge() { }
        /// <summary>
        /// Gets or sets the transparency value [0, 255].
        /// </summary>
        public int Transparency
        {
            get
            {
                return this.transparency;
            }
            set
            {
                this.transparency = value;
            }
        }
        /// <summary>
        /// Gets or sets a pair of integers representing an ordered pair of X and Y coordinates.
        /// </summary>
        public PointInt Point
        {
            get
            {
                return point;
            }
            set
            {
                point = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image:
            int srcStride = bmSrc.Stride, dstStride = bmData.Stride;

            // source pixel's coordinates from rectangle:
            int startX = (int)Maths.Max(point.X, 0);
            int startY = (int)Maths.Max(point.Y, 0);
            int endX   = (int)Maths.Min(bmSrc.Width,  bmData.Width - startX);
            int endY   = (int)Maths.Min(bmSrc.Height, bmData.Height - startY);

            // do the job:
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            // transparency:
            double t = transparency / 255.0;
            int a0, a1;

            // pixel offsets:
            int x, y, k, l;

            // for each line:
            for (y = 0; y < endY; y++)
            {
                // for each pixel:
                for (x = 0; x < endX; x++)
                {
                    // local data:
                    k = dstStride * (y + startY) + 4 * (x + startX);
                    l = srcStride * (y         ) + 4 * (x         );

                    // calculating transparency:
                    a1 = (int)(src[3] * t); a0 = 255 - a1;

                    // applying filter:
                    dst[k + 0] = merge(dst[k + 0], src[l + 0], a0, a1);
                    dst[k + 1] = merge(dst[k + 1], src[l + 1], a0, a1);
                    dst[k + 2] = merge(dst[k + 2], src[l + 2], a0, a1);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        #endregion

        #region Merging function components
        /// <summary>
        /// Merge function.
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="a0"></param>
        /// <param name="a1"></param>
        /// <returns></returns>
        private static byte merge(byte x, byte y, int a0, int a1)
        {
            return Maths.Byte((x * a0 + y * a1) / 255);
        }
        #endregion
    }
    /// <summary>
    /// Defines the texturing filter.
    /// </summary>
    public class Texturer : IBitmapFilter
    {
        #region Private data
        private double[,] texture = null;
        private double depth = 0.5;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the texturing filter.
        /// </summary>
        /// <param name="texture">Matrix</param>
        public Texturer(double[,] texture)
        {
            Texture = texture;
        }
        /// <summary>
        /// Initializes the texturing filter.
        /// </summary>
        /// <param name="texture">Matrix</param>
        /// <param name="depth">Depth [0, 1]</param>
        public Texturer(double[,] texture, double depth = 1.0)
        {
            Texture = texture; Depth = depth;
        }
        /// <summary>
        /// Gets or sets the texture matrix.
        /// </summary>
        public double[,] Texture
        {
            get { return texture; }
            set { texture = value; }
        }
        /// <summary>
        /// Gets or sets the depth value [0, 1].
        /// </summary>
        public double Depth
        {
            get { return this.depth; }
            set { this.depth = Maths.Double(value); }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int widthToProcess = Math.Min(width, texture.GetLength(1));
            int heightToProcess = Math.Min(height, texture.GetLength(0));
            double z = 1.0 - this.depth;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, heightToProcess, y =>
            {
                int x, ystride, k;
                byte red, green, blue;
                double t;

                ystride = y * stride;

                for (x = 0; x < widthToProcess; x++)
                {
                    k = ystride + x * 4;
                    red = p[k + 2]; green = p[k + 1]; blue = p[k];
                    t = texture[y, x];

                    p[k + 2] = Maths.Byte((z * red) + (this.depth * red) * t);
                    p[k + 1] = Maths.Byte((z * green) + (this.depth * green) * t);
                    p[k] = Maths.Byte((z * blue) + (this.depth * blue) * t);
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        private static Random rand = new Random();
        private static int r;

        /// <summary>
        /// Implements the construction of a wood texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="rings">Rings</param>
        /// <returns>Matrix</returns>
        public static Texturer Wood(int m, int l, double rings = 12)
        {
            PerlinNoise noise = new PerlinNoise(8, 0.5, 1.0 / 32, 0.05); r = rand.Next(5000);
            double[,] texture = new double[m, l];
            int w2 = l / 2, h2 = m / 2;

            Parallel.For(0, m, y =>
            {
                double xv, yv;
                int x;

                for (x = 0; x < l; x++)
                {
                    xv = (double)(x - w2) / l;
                    yv = (double)(y - h2) / m;

                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(Math.Sin((Math.Sqrt(xv * xv + yv * yv) + noise.Function2D(x + r, y + r)) * Math.PI * 2 * rings)));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a textile texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Textile(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(3, 0.65, 1.0 / 8, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;
                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Max(0.0f, Math.Min(1.0f, (
                                (double)Math.Sin(x + noise.Function2D(x + r, y + r)) +
                                (double)Math.Sin(y + noise.Function2D(x + r, y + r))) * 0.25f + 0.5f));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a marble texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="yPeriod">Y-period</param>
        /// <param name="xPeriod">X-period</param>
        /// <returns>Matrix</returns>
        public static Texturer Marble(int m, int l, double yPeriod = 10.0, double xPeriod = 5.0)
        {
            PerlinNoise noise = new PerlinNoise(2, 0.65, 1.0 / 32, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];
            double xFact = xPeriod / l;
            double yFact = yPeriod / m;

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(Math.Sin((x * xFact + y * yFact + noise.Function2D(x + r, y + r)) * Math.PI)));
                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a labyrinth texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Labyrinth(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(1, 0.65, 1.0 / 16, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Min(1.0f, (float)Math.Abs(noise.Function2D(x + r, y + r)));

                }
            }
            );

            return new Texturer(texture);
        }
        /// <summary>
        /// Implements the construction of a clouds texture.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Texturer Clouds(int m, int l)
        {
            PerlinNoise noise = new PerlinNoise(8, 0.5, 1.0 / 32, 1.0); r = rand.Next(5000);
            double[,] texture = new double[m, l];

            Parallel.For(0, m, y =>
            {
                int x;

                for (x = 0; x < l; x++)
                {
                    texture[y, x] = Math.Max(0.0f, Math.Min(1.0f, (double)noise.Function2D(x + r, y + r) * 0.5f + 0.5f));

                }
            }
            );

            return new Texturer(texture);
        }
        #endregion
    }
    #endregion

    #region Operations
    /// <summary>
    /// Defines the linear operation filter.
    /// <remarks>
    /// This filter works according to the following algorithm: C (x, y) = a * A (x, y) + b * B (x, y), where A, B are the original images,
    /// a, b are the coefficients.
    /// </remarks>
    /// </summary>
    public class Operation : IBitmapFilter2
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the linear operation filter.
        /// </summary>
        /// <param name="a">First image coefficient</param>
        /// <param name="b">Second image coefficient</param>
        public Operation(double a, double b)
        {
            A = a; B = b; 
        }
        /// <summary>
        /// Gets or sets the first image coefficient.
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
            }
        }
        /// <summary>
        /// Gets or sets the second image coefficient.
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
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (x = 0; x < width; x++)
            {
                for (y = 0; y < height; y++, p += 4, pSrc += 4)
                {
                    p[2] = Maths.Byte(a * p[2] + b * pSrc[2]);
                    p[1] = Maths.Byte(a * p[1] + b * pSrc[1]);
                    p[0] = Maths.Byte(a * p[0] + b * pSrc[0]);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        #endregion

        #region Public static methods
        /// <summary>
        /// Addition filter.
        /// </summary>
        public static Operation Addition
        {
            get
            {
                return new Operation(1, 1);
            }
        }
        /// <summary>
        /// Subtraction filter.
        /// </summary>
        public static Operation Subtraction
        {
            get
            {
                return new Operation(1, -1);
            }
        }
        /// <summary>
        /// Averaging filter.
        /// </summary>
        public static Operation Averaging
        {
            get
            {
                return new Operation(0.5, 0.5);
            }
        }
        #endregion
    }
    /// <summary>
    /// Defines the stereo effect filter for a pair of images.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
    /// </remarks>
    /// </summary>
    public class StereoAnaglyph : IBitmapFilter2
    {
        #region Private data
        private Anaglyph algorithm = Anaglyph.Gray;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the stereo effect filter for a pair of images.
        /// </summary>
        /// <param name="algorithm">Algorithm</param>
        public StereoAnaglyph(Anaglyph algorithm)
        {
            this.algorithm = algorithm;
        }
        /// <summary>
        /// Gets or sets the algorithm.
        /// </summary>
        public Anaglyph Algorithm
        {
            get { return algorithm; }
            set { algorithm = value; }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            int width = bmData.Width, height = bmData.Height;
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)bmSrc.Scan0.ToPointer();
            int x, y;

            switch (algorithm)
            {
                case Anaglyph.True:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = 0;
                            p[0] = (byte)(pSrc[2] * 0.299 + pSrc[1] * 0.587 + pSrc[0] * 0.114);
                        }
                    }
                    break;

                case Anaglyph.Gray:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = (byte)(pSrc[2] * 0.299 + pSrc[1] * 0.587 + pSrc[0] * 0.114);
                            p[0] = p[1];
                        }
                    }
                    break;

                case Anaglyph.Color:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            // keep Red as it is and take only Green and Blue from the second image
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;

                case Anaglyph.HalfColor:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114);
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;

                case Anaglyph.Optimized:
                    // for each line
                    for (y = 0; y < height; y++)
                    {
                        // for each pixel
                        for (x = 0; x < width; x++, p += 4, pSrc += 4)
                        {
                            p[2] = (byte)(p[1] * 0.7 + p[0] * 0.3);
                            p[1] = pSrc[1];
                            p[0] = pSrc[0];
                        }
                    }
                    break;
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        #endregion

        #region Enums
        /// <summary>
        /// Defines the stereo effect creation algorithm.
        /// </summary>
        /// <remarks>
        /// More information can be found on the website:
        /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
        /// </remarks>
        public enum Anaglyph
        {
            /// <summary>
            /// Creates a stereo effect for a pair of images according to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            True,
            /// <summary>
            /// Creates a stereo effect for a pair of images according to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Gray,
            /// <summary>
            /// Creates a stereo effect for a pair of images according to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=R<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Color,
            /// <summary>
            /// Creates a stereo effect for a pair of images according to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            HalfColor,
            /// <summary>
            /// Creates a stereo effect for a pair of images according to the following calculations:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.7*G<sub>l</sub>+0.3*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Optimized
        }
        #endregion
    }
    /// <summary>
    /// Defines the oil filter.
    /// <remarks>
    /// More information can be found on the website:
    /// https://www.codeproject.com/articles/471994/oilpainteffect
    /// </remarks>
    /// </summary>
    public class OilPainting : IBitmapFilter2
    {
        #region Private data
        private double depth;
        private int radius0;
        private int radius1;
        private int l0;
        private int l1;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(int radius = 3, double depth = 1.0)
        {
            Size = new SizeInt(radius, radius);
            Depth = depth;
        }
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="height">Filter height</param>
        /// <param name="width">Filter width</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(int width, int height, double depth = 1.0)
        {
            Size = new SizeInt(width, height);
            Depth = depth;
        }
        /// <summary>
        /// Initializes the oil filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="depth">Value [0, 1]</param>
        public OilPainting(SizeInt size, double depth = 1.0)
        {
            Size = size;
            Depth = depth;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(this.l0, this.l1);
            }
            set
            {
                this.l0 = value.Width;
                this.l1 = value.Height;
                this.radius0 = l0 >> 1;
                this.radius1 = l1 >> 1;
            }
        }
        /// <summary>
        /// Gets or sets the depth value [0, 1].
        /// </summary>
        public double Depth
        {
            get
            {
                return this.depth;
            }
            set
            {
                this.depth = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            double strenght = this.depth * 255.0;
            double strenghtGlobal = this.depth / 3.0;
            #endregion

            Parallel.For(0, height, y =>
            {
                int[] Red = new int[256];
                int[] Green = new int[256];
                int[] Blue = new int[256];
                int[] Intensity = new int[256];

                int red, green, blue, intensity, max = 0, index = 0;
                int x, i, j, ir, jr, yr, xr, n;
                int ystride, v;
                byte* p;

                yr = y - radius0;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - radius1;
                    v = ystride + x * 4;
                    red = green = blue = 0;

                    #region Convolution filtering
                    for (i = 0; i < l0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        for (j = 0; j < l1; j++)
                        {
                            jr = xr + j;

                            if (jr < 0) continue; if (jr >= width) break;

                            p = &src[ir * stride + jr * 4];
                            red = p[2]; green = p[1]; blue = p[0];
                            intensity = (int)((red + green + blue) * strenghtGlobal);

                            Red[intensity] += red;
                            Green[intensity] += green;
                            Blue[intensity] += blue;
                            Intensity[intensity]++;
                        }
                    }
                    #endregion

                    #region Frequent intensity
                    for (n = 0; n < strenght; n++)
                    {
                        if (Intensity[n] > max)
                        {
                            max = Intensity[n];
                            index = n;
                        }
                    }
                    #endregion

                    #region Result pixel calculation
                    if (max != 0)
                    {
                        blue = Blue[index] / max;
                        green = Green[index] / max;
                        red = Red[index] / max;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = (byte)red;
                    dst[v + 1] = (byte)green;
                    dst[v] = (byte)blue;
                    #endregion

                    #region Clear data
                    Array.Clear(Red, 0, 256);
                    Array.Clear(Green, 0, 256);
                    Array.Clear(Blue, 0, 256);
                    Array.Clear(Intensity, 0, 256);
                    max = 0; index = 0;
                    #endregion
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the convolution filter.
    /// </summary>
    public class Convolution : IBitmapFilter2
    {
        #region Private data
        private double[][] kernel;
        private double offset;
        private int radius0;
        private int radius1;
        private int l0;
        private int l1;
        private bool bilateral;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the convolution filter.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="offset">Offset</param>
        /// <param name="bilateral">Bilateral processing or not</param>
        public Convolution(double[,] m, double offset = 0, bool bilateral = false)
        {
            Matrix = m; Offset = offset; Bilateral = bilateral;
        }
        /// <summary>
        /// Initializes the convolution filter.
        /// </summary>
        public Convolution()
        {
            Matrix = Matrice.One(3, 3); Offset = 0; Bilateral = false;
        }
        /// <summary>
        /// Gets or sets the convolution matrix.
        /// </summary>
        public double[,] Matrix
        {
            get
            {
                return Jagged.FromJagged(this.kernel);
            }
            set
            {
                Data(value);
            }
        }
        /// <summary>
        /// Gets or sets the offset value.
        /// </summary>
        public double Offset
        {
            get
            {
                return this.offset;
            }
            set
            {
                this.offset = value;
            }
        }
        /// <summary>
        /// Bilateral processing or not.
        /// </summary>
        public bool Bilateral
        {
            get
            {
                return this.bilateral;
            }
            set
            {
                this.bilateral = value;
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="m">Matrix</param>
        private void Data(double[,] m)
        {
            this.l0 = m.GetLength(0);
            this.l1 = m.GetLength(1);
            this.radius0 = l0 >> 1;
            this.radius1 = l1 >> 1;
            this.kernel = Jagged.ToJagged(m);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            if (this.l0 != this.l1 && this.bilateral == true)
                throw new Exception("Matrix must be squared");

            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            #endregion

            if (!this.bilateral)
            {
                Parallel.For(0, height, y =>
                {
                    double red, green, blue, div, k;
                    int x, i, j, ir, jr, yr, xr;
                    int ystride, v;
                    byte* p;

                    yr = y - radius0;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        xr = x - radius1;
                        v = ystride + x * 4;

                        red = green = blue = div = 0;

                        #region Convolution filtering
                        for (i = 0; i < l0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < l1; j++)
                            {
                                jr = xr + j;

                                if (jr < 0) continue; if (jr >= width) break;

                                k = kernel[i][j];

                                if (k != 0)
                                {
                                    p = &src[ir * stride + jr * 4];
                                    div += k;
                                    red += k * p[2]; green += k * p[1]; blue += k * p[0];
                                }
                            }
                        }
                        #endregion

                        #region Divider and offset
                        if (div != 0)
                        {
                            red /= div;
                            green /= div;
                            blue /= div;
                        }
                        if (offset != 0)
                        {
                            red += offset;
                            green += offset;
                            blue += offset;
                        }
                        #endregion

                        #region Recording pixel
                        dst[v + 2] = Maths.Byte(red);
                        dst[v + 1] = Maths.Byte(green);
                        dst[v] = Maths.Byte(blue);
                        #endregion
                    }
                });
            }
            else
            {
                Parallel.For(0, height, y =>
                {
                    double red1, green1, blue1, div1, k1;
                    double red2, green2, blue2, div2, k2;
                    int x, i, j, ir, jr, yr, xr;
                    int ystride, v;
                    byte* p;

                    yr = y - radius0;
                    ystride = y * stride;

                    for (x = 0; x < width; x++)
                    {
                        xr = x - radius1;
                        v = ystride + x * 4;

                        red1 = green1 = blue1 = div1 = 0;
                        red2 = green2 = blue2 = div2 = 0;

                        #region Convolution filtering
                        for (i = 0; i < l0; i++)
                        {
                            ir = yr + i;
                            if (ir < 0) continue; if (ir >= height) break;

                            for (j = 0; j < l1; j++)
                            {
                                jr = xr + j;

                                if (jr < 0) continue; if (jr >= width) break;

                                k1 = kernel[i][j]; k2 = kernel[j][i];

                                p = &src[ir * stride + jr * 4];

                                div1 += k1; div2 += k2;
                                red1 += k1 * p[2]; green1 += k1 * p[1]; blue1 += k1 * p[0];
                                red2 += k2 * p[2]; green2 += k2 * p[1]; blue2 += k2 * p[0];
                            }
                        }
                        #endregion

                        #region Divider and offset
                        if (div1 != 0)
                        {
                            red1 /= div1;
                            green1 /= div1;
                            blue1 /= div1;
                        }
                        if (div2 != 0)
                        {
                            red2 /= div2;
                            green2 /= div2;
                            blue2 /= div2;
                        }
                        if (offset != 0)
                        {
                            red1 += offset;
                            green1 += offset;
                            blue1 += offset;

                            red2 += offset;
                            green2 += offset;
                            blue2 += offset;
                        }
                        #endregion

                        #region Recording pixel
                        dst[v + 2] = Maths.Byte(G(red1, red2));
                        dst[v + 1] = Maths.Byte(G(green1, green2));
                        dst[v] = Maths.Byte(G(blue1, blue2));
                        #endregion
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
        }
        #endregion

        #region Sobel's gradient components
        /// <summary>
        /// Gets the value of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Double precision floating point number</returns>
        public static double G(double Gx, double Gy)
        {
            return Math.Sqrt(Gx * Gx + Gy * Gy);
        }
        /// <summary>
        /// Gets the angle of the gradient operator.
        /// </summary>
        /// <param name="Gx">Gradient X</param>
        /// <param name="Gy">Gradient Y</param>
        /// <returns>Double precision floating point number</returns>
        public static double Tetta(double Gx, double Gy)
        {
            return Math.Atan(Gx / Gy);
        }
        #endregion

        #region Public static methods
        #region Radius matrix
        /// <summary>
        /// mplements the construction of the inverted Gausssian filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution LoGaussian(int m, int l, double sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            double[,] H = new double[m, l];
            double sigma2 = sigma * sigma;
            double f0 = -1.0 / (Math.PI * sigma2 * sigma2);
            double f1 = 2.0 * sigma2;
            double kernel;
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    kernel = (x * x + y * y) / f1;
                    H[i, j] = f0 * (1.0 - kernel) * Math.Exp(-kernel);
                }
            }
            return new Convolution(H);
        }
        /// <summary>
        /// mplements the construction of the Gaussian blur filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Gaussian(int m, int l, double sigma)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            double[,] H = new double[m, l];
            int i, j, x, y;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    H[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }
            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the "unsharp masking" filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="sigma">Standard deviation (!=0)</param>
        /// <returns>Matrix</returns>
        public static Convolution Unsharp(int m, int l, double sigma)
        {
            double[,] G = new double[m, l];
            int i, j, x, y;
            int r1 = m / 2;
            int r2 = l / 2;

            for (y = -r1, i = 0; i < m; y++, i++)
            {
                for (x = -r2, j = 0; j < l; x++, j++)
                {
                    G[i, j] = Kernel.Gaussian(x, sigma) * Kernel.Gaussian(y, sigma);
                }
            }

            double[,] invG = new double[m, l];
            double max = G[0, 0];
            double summary = 0;
            double v, iv;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    v = G[i, j] / max;
                    iv = v > byte.MaxValue ? byte.MaxValue : (int)v;
                    invG[i, j] = -iv;
                    summary += iv;
                }
            }

            invG[m / 2, l / 2] = 2 * summary - invG[m / 2, l / 2];
            return new Convolution(invG);
        }
        /// <summary>
        /// Implements the construction of the high-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <param name="boost">Boost</param>
        /// <returns>Matrix</returns>
        public static Convolution HighPass(int m, int l, double boost)
        {
            int r1 = m / 2;
            int r2 = l / 2;
            double[,] H = new double[m, l];
            int i, j;

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < l; j++)
                {
                    H[i, j] = -1;
                }
            }
            H[r1, r2] = boost;
            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the low-pass filter.
        /// </summary>
        /// <param name="m">Height</param>
        /// <param name="l">Width</param>
        /// <returns>Matrix</returns>
        public static Convolution LowPass(int m, int l)
        {
            return Convolution.HighPass(m, l, 1);
        }
        /// <summary>
        /// Implements the construction of the emboss filter.
        /// </summary>
        /// <param name="n">Size</param>
        /// <returns>Matrix</returns>
        public static Convolution Emboss(int n)
        {
            double[,] H = new double[n, n];
            int r = n - 1, r2 = r / 2;

            H[0, 0] = -2; H[0, r2] = -1; //         0;
            H[r2, 0] = -1; H[r2, r2] = 1; H[r2, r] = 1;
            H[r, r2] = 1; H[r, r] = 2; //         0;

            return new Convolution(H);
        }
        #endregion

        #region Fixed radius matrix
        /// <summary>
        /// Implements the construction of the Roberts operator [2 x 2].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Roberts()
        {
            return new Convolution(new double[2, 2] { { 1, 0 }, { 0, -1 } });
        }
        /// <summary>
        /// Implements the construction of the Prewitt operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Prewitt()
        {
            return new Convolution(new double[3, 3] { { -1, -1, -1 }, { 0, 0, 0 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Implements the construction of the Sobel operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Sobel()
        {
            return new Convolution(new double[3, 3] { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } });
        }
        /// <summary>
        /// Implements the construction of the Scharr operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Scharr()
        {
            return new Convolution(new double[3, 3] { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } });
        }
        /// <summary>
        /// Implements the construction of the Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution Laplacian()
        {
            return new Convolution(new double[3, 3] { { 0, 1, 0 }, { 1, -4, 1 }, { 0, 1, 0 } });
        }
        /// <summary>
        /// Implements the construction of the diagonal Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianDiagonal()
        {
            return new Convolution(new double[3, 3] { { 1, 1, 1 }, { 1, -8, 1 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Implements the construction of the inverted Laplacian operator [3 x 3].
        /// </summary>
        /// <returns>Matrix</returns>
        public static Convolution LaplacianInvert()
        {
            return new Convolution(new double[3, 3] { { -1, 0, -1 }, { 0, 4, 0 }, { -1, 0, -1 } });
        }
        #endregion

        #region Fixed radius compass matrix
        /// <summary>
        /// Implements the construction of the Kirsh operator [3 x 3].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Kirsh(Gradient direction)
        {
            double[,] H = new double[3, 3];

            if (direction == Gradient.North)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.NorthWest)
            {
                H[0, 0] = 5; H[0, 1] = 5; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 5; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = -3; H[2, 2] = -3;
            }
            else if (direction == Gradient.SouthWest)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = 5; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = -3;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = -3;
                H[2, 0] = 5; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.SouthEast)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = -3;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = 5; H[2, 2] = 5;
            }
            else if (direction == Gradient.East)
            {
                H[0, 0] = -3; H[0, 1] = -3; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = 5;
            }
            else if (direction == Gradient.NorthEast)
            {
                H[0, 0] = -3; H[0, 1] = 5; H[0, 2] = 5;
                H[1, 0] = -3; H[1, 1] = 0; H[1, 2] = 5;
                H[2, 0] = -3; H[2, 1] = -3; H[2, 2] = -3;
            }

            return new Convolution(H);
        }
        /// <summary>
        /// Implements the construction of the Roberts operator [3 x 3]. [2 x 2].
        /// </summary>
        /// <param name="direction">Gradient direction</param>
        /// <returns>Matrix</returns>
        public static Convolution Roberts(Gradient direction)
        {
            double[,] H = new double[2, 2];

            if (direction == Gradient.North)
            {
                H[0, 0] = 0; H[0, 1] = 1;
                H[1, 0] = -1; H[1, 1] = 0;
            }
            else if (direction == Gradient.West)
            {
                H[0, 0] = 1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = -1;
            }
            else if (direction == Gradient.South)
            {
                H[0, 0] = 0; H[0, 1] = -1;
                H[1, 0] = 1; H[1, 1] = 0;
            }
            else // Properties.Gradient.East
            {
                H[0, 0] = -1; H[0, 1] = 0;
                H[1, 0] = 0; H[1, 1] = 1;
            }

            return new Convolution(H);
        }
        #endregion
        #endregion
    }
    /// <summary>
    /// Defines the morphology filter.
    /// </summary>
    public class Morphology : IBitmapFilter2
    {
        #region Private data
        private int threshold;
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(int radius = 3, int threshold = 0)
        {
            Size = new SizeInt(radius, radius);
            Threshold = threshold;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(int width, int height, int threshold)
        {
            Size = new SizeInt(width, height);
            Threshold = threshold;
        }
        /// <summary>
        /// Initializes the morphology filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        /// <param name="threshold">Threshold</param>
        public Morphology(SizeInt size, int threshold = 0)
        {
            Size = size;
            Threshold = threshold;
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
        /// Gets or sets the threshold.
        /// </summary>
        public int Threshold
        {
            get
            {
                return this.threshold;
            }
            set
            {
                this.threshold = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyHorizontal(bmSrc, bmData);
                ApplyVertical(bmData, bmSrc);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData, bmSrc);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData, bmSrc);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rh >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte[] red, green, blue;
                int x, i, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y - r;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x;
                    v = ystride + x * 4;

                    red = new byte[rh];
                    green = new byte[rh];
                    blue = new byte[rh];

                    #region Convolution filtering
                    for (i = 0; i < rh; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = xr + 0;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red[i] = p[2];
                        green[i] = p[1];
                        blue[i] = p[0];
                    }
                    #endregion

                    #region Morphology filtering
                    Array.Sort(red);
                    Array.Sort(green);
                    Array.Sort(blue);
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red[threshold];
                    dst[v + 1] = green[threshold];
                    dst[v] = blue[threshold];
                    #endregion
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rw >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte[] red, green, blue;
                int x, j, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - r;
                    v = ystride + x * 4;

                    red = new byte[rw];
                    green = new byte[rw];
                    blue = new byte[rw];

                    #region Convolution filtering
                    ir = yr + 0;
                    if (ir < 0) continue; if (ir >= height) break;

                    for (j = 0; j < rw; j++)
                    {
                        jr = xr + j;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red[j] = p[2];
                        green[j] = p[1];
                        blue[j] = p[0];
                    }
                    #endregion

                    #region Morphology filtering
                    Array.Sort(red);
                    Array.Sort(green);
                    Array.Sort(blue);
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red[threshold];
                    dst[v + 1] = green[threshold];
                    dst[v] = blue[threshold];
                    #endregion
                }
            }
            );

            return;
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Initializes the median filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Median(int radius)
        {
            int threshold = radius / 2;
            return new Morphology(radius, radius, threshold);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Erosion(int radius)
        {
            return new Morphology(radius, radius, 0);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public static Morphology Dilatation(int radius)
        {
            return new Morphology(radius, radius, radius - 1);
        }
        #endregion
    }
    /// <summary>
    /// Defines the dilatation filter.
    /// </summary>
    public class Dilatation : IBitmapFilter2
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Dilatation(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Dilatation(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the dilatation filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Dilatation(SizeInt size)
        {
            Size = size;
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
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyHorizontal(bmSrc, bmData);
                ApplyVertical(bmData, bmSrc);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData, bmSrc);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData, bmSrc);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rh >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, i, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y - r;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x;
                    v = ystride + x * 4;

                    red = green = blue = 0;

                    #region Convolution filtering
                    for (i = 0; i < rh; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = xr + 0;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] > red ? p[2] : red;
                        green = p[1] > green ? p[1] : green;
                        blue = p[0] > blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rw >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, j, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - r;
                    v = ystride + x * 4;

                    red = green = blue = 0;

                    #region Convolution filtering
                    ir = yr + 0;
                    if (ir < 0) continue; if (ir >= height) break;

                    for (j = 0; j < rw; j++)
                    {
                        jr = xr + j;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] > red ? p[2] : red;
                        green = p[1] > green ? p[1] : green;
                        blue = p[0] > blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the erosion filter.
    /// </summary>
    public class Erosion : IBitmapFilter2
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Erosion(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Erosion(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the erosion filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Erosion(SizeInt size)
        {
            Size = size;
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
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyHorizontal(bmSrc, bmData);
                ApplyVertical(bmData, bmSrc);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData, bmSrc);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData, bmSrc);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rh >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, i, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y - r;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x;
                    v = ystride + x * 4;

                    red = green = blue = 255;

                    #region Convolution filtering
                    for (i = 0; i < rh; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;

                        jr = xr + 0;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] < red ? p[2] : red;
                        green = p[1] < green ? p[1] : green;
                        blue = p[0] < blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int r = rw >> 1;
            #endregion

            Parallel.For(0, height, y =>
            {
                byte red, green, blue;
                int x, j, ir, jr, yr, xr;
                int ystride, v;
                byte* p;

                yr = y;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - r;
                    v = ystride + x * 4;

                    red = green = blue = 255;

                    #region Convolution filtering
                    ir = yr + 0;
                    if (ir < 0) continue; if (ir >= height) break;

                    for (j = 0; j < rw; j++)
                    {
                        jr = xr + j;

                        if (jr < 0) continue; if (jr >= width) break;

                        p = &src[ir * stride + jr * 4];
                        red = p[2] < red ? p[2] : red;
                        green = p[1] < green ? p[1] : green;
                        blue = p[0] < blue ? p[0] : blue;
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = red;
                    dst[v + 1] = green;
                    dst[v] = blue;
                    #endregion
                }
            }
            );

            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the top-hat filter.
    /// </summary>
    public class TopHat : IBitmapFilter2
    {
        #region Private data
        private Opening opening = new Opening();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the top-hat filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public TopHat(int radius = 3)
        {
            opening = new Opening(radius);
        }
        /// <summary>
        /// Initializes the top-hat filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public TopHat(int width, int height)
        {
            opening = new Opening(width, height);
        }
        /// <summary>
        /// Initializes the top-hat filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public TopHat(SizeInt size)
        {
            opening = new Opening(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return opening.Size;
            }
            set
            {
                opening.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            Bitmap Src0 = (Bitmap)BitmapConverter.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapConverter.Lock32bpp(Src0);

            // Filter applying:
            opening.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapConverter.Unlock(Src0, bmSrc0);
            Src0.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the bottom-hat filter.
    /// </summary>
    public class BottomHat : IBitmapFilter2
    {
        #region Private data
        private Closing closing = new Closing();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public BottomHat(int radius = 3)
        {
            closing = new Closing(radius);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public BottomHat(int width, int height)
        {
            closing = new Closing(width, height);
        }
        /// <summary>
        /// Initializes the bottom-hat filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public BottomHat(SizeInt size)
        {
            closing = new Closing(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return closing.Size;
            }
            set
            {
                closing.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            Bitmap Src0 = (Bitmap)BitmapConverter.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapConverter.Lock32bpp(Src0);

            // Filter applying:
            closing.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapConverter.Unlock(Src0, bmSrc0);
            Src0.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the closing filter.
    /// </summary>
    public class Closing : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Dilatation dilatation = new Dilatation();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Closing(int radius = 3)
        {
            erosion = new Erosion(radius);
            dilatation = new Dilatation(radius);
        }
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Closing(int width, int height)
        {
            erosion = new Erosion(width, height);
            dilatation = new Dilatation(width, height);
        }
        /// <summary>
        /// Initializes the closing filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Closing(SizeInt size)
        {
            erosion = new Erosion(size);
            dilatation = new Dilatation(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return erosion.Size;
            }
            set
            {
                erosion.Size = value;
                dilatation.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            dilatation.Apply(bmSrc, bmData);
            erosion.Apply(bmData, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the opening filter.
    /// </summary>
    public class Opening : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Dilatation dilatation = new Dilatation();
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the opening filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public Opening(int radius = 3)
        {
            erosion = new Erosion(radius);
            dilatation = new Dilatation(radius);
        }
        /// <summary>
        /// Initializes the opening filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public Opening(int width, int height)
        {
            erosion = new Erosion(width, height);
            dilatation = new Dilatation(width, height);
        }
        /// <summary>
        /// Initializes the opening filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public Opening(SizeInt size)
        {
            erosion = new Erosion(size);
            dilatation = new Dilatation(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return erosion.Size;
            }
            set
            {
                erosion.Size = value;
                dilatation.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            erosion.Apply(bmSrc, bmData);
            dilatation.Apply(bmData, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the edge glow filter.
    /// </summary>
    public class EdgeGlow : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public EdgeGlow(int radius = 3)
        {
            erosion = new Erosion(radius);
        }
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public EdgeGlow(int width, int height)
        {
            erosion = new Erosion(width, height);
        }
        /// <summary>
        /// Initializes the edge glow filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public EdgeGlow(SizeInt size)
        {
            erosion = new Erosion(size);
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return this.erosion.Size;
            }
            set
            {
                this.erosion.Size = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // Creating resources:
            Bitmap Src0 = (Bitmap)BitmapConverter.Bitmap(bmSrc).Clone();
            BitmapData bmSrc0 = BitmapConverter.Lock32bpp(Src0);

            // Filter applying:
            erosion.Apply(bmSrc, bmSrc0);
            subtraction.Apply(bmData, bmSrc);

            // Delete resources:
            BitmapConverter.Unlock(Src0, bmSrc0);
            Src0.Dispose();
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
        }
        #endregion
    }
    /// <summary>
    /// Defines the box blur filter.
    /// </summary>
    public class BoxBlur : IBitmapFilter
    {
        #region Private data
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public BoxBlur(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public BoxBlur(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the box blur filter.
        /// </summary>
        /// <param name="size">Filter size</param>
        public BoxBlur(SizeInt size)
        {
            Size = size;
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
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            if (rw >= 2 && rh >= 2)
            {
                ApplyHorizontal(bmData);
                ApplyVertical(bmData);
            }
            else if (rw >= 2 && rh < 2)
            {
                ApplyHorizontal(bmData);
            }
            else if (rw < 2 && rh >= 2)
            {
                ApplyVertical(bmData);
            }
            else return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyVertical(BitmapData bmData)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rh >= height ? height - 1 : rh;
            int v = h >> 1;
            int dl = height - v;
            #endregion

            Parallel.For(0, width, x =>
            {
                double r = 0;
                double g = 0;
                double b = 0;
                int p, w, q, y;
                int xx = x * 4;

                for (p = xx, y = 0; y < h; y++, p += stride)
                {
                    r += dst[p + 2];
                    g += dst[p + 1];
                    b += dst[p + 0];
                }

                for (p = xx, y = 0; y < v; y++, p += stride)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    y = v,
                    p = xx,
                    q = xx + (y + 0) * stride,
                    w = xx + (y + v) * stride;

                    y < dl;

                    y++,
                    p += stride,
                    q += stride,
                    w += stride)
                {
                    r = r - dst[p + 2] + dst[w + 2];
                    g = g - dst[p + 1] + dst[w + 1];
                    b = b - dst[p + 0] + dst[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    y = dl,
                    w = xx + (y - v) * stride,
                    p = xx + (y + 0) * stride;

                    y < height;

                    y++,
                    w += stride,
                    p += stride)
                {
                    r = r - dst[w + 2] + dst[p + 2];
                    g = g - dst[w + 1] + dst[p + 1];
                    b = b - dst[w + 0] + dst[p + 0];

                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

            });

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHorizontal(BitmapData bmData)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            int h = rw >= width ? width - 1 : rw;
            int v = h >> 1;
            int dl = width - v;
            #endregion

            Parallel.For(0, height, y =>
            {
                double r = 0;
                double g = 0;
                double b = 0;
                int p, q, w, x;
                int yy = y * stride;

                for (p = yy, x = 0; x < h; x++, p += 4)
                {
                    r += dst[p + 2];
                    g += dst[p + 1];
                    b += dst[p + 0];
                }

                for (p = yy, x = 0; x < v; x++, p += 4)
                {
                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

                for (
                    x = v,
                    p = yy,
                    q = yy + (x + 0) * 4,
                    w = yy + (x + v) * 4;

                    x < dl;

                    x++,
                    p += 4,
                    q += 4,
                    w += 4)
                {
                    r = r - dst[p + 2] + dst[w + 2];
                    g = g - dst[p + 1] + dst[w + 1];
                    b = b - dst[p + 0] + dst[w + 0];

                    dst[q + 2] = Maths.Byte(r / h);
                    dst[q + 1] = Maths.Byte(g / h);
                    dst[q + 0] = Maths.Byte(b / h);
                }

                for (
                    x = dl,
                    w = (x - v) * 4 + yy,
                    p = (x + 0) * 4 + yy;

                    x < width;

                    x++,
                    p += 4,
                    w += 4)
                {
                    r = r - dst[w + 2] + dst[p + 2];
                    g = g - dst[w + 1] + dst[p + 1];
                    b = b - dst[w + 0] + dst[p + 0];

                    dst[p + 2] = Maths.Byte(r / h);
                    dst[p + 1] = Maths.Byte(g / h);
                    dst[p + 0] = Maths.Byte(b / h);
                }

            });
        }
        #endregion
    }
    #endregion

    #region Colorspace filters
    /// <summary>
    /// Defines the HSB filter.
    /// </summary>
    public class HSBFilter : IBitmapFilter
    {
        #region Private data
        private int hue;
        private double saturation;
        private double brightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the HSB filter.
        /// </summary>
        /// <param name="hue">Hue [0, 359]</param>
        /// <param name="saturation">Saturation [-1, 1]</param>
        /// <param name="brightness">Brightness [-1, 1]</param>
        public HSBFilter(int hue, double saturation, double brightness)
        {
            Hue = hue;
            Saturation = saturation;
            Brightness = brightness;
        }
        /// <summary>
        /// Initializes the HSB filter.
        /// </summary>
        public HSBFilter()
        {
            new HSBFilter(0, 0, 0);
        }
        /// <summary>
        /// Hue [0, 359].
        /// </summary>
        public int Hue
        {
            get
            {
                return this.hue;
            }
            set
            {
                this.hue = value;
            }
        }
        /// <summary>
        /// Saturation [-1, 1].
        /// </summary>
        public double Saturation
        {
            get
            {
                return this.saturation;
            }
            set
            {
                this.saturation = value;
            }
        }
        /// <summary>
        /// Brightness [-1, 1].
        /// </summary>
        public double Brightness
        {
            get
            {
                return this.brightness;
            }
            set
            {
                this.brightness = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSB hsb; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    hsb = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    hsb.Hue = hue;
                    hsb.Saturation += saturation;
                    hsb.Brightness += brightness;
                    rgb = hsb.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the HSL filter.
    /// </summary>
    public class HSLFilter : IBitmapFilter
    {
        #region Private data
        private int hue;
        private double saturation;
        private double lightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the HSL filter.
        /// </summary>
        /// <param name="hue">Hue [0, 359]</param>
        /// <param name="saturation">Saturation [-1, 1]</param>
        /// <param name="lightness">Lightness [-1, 1]</param>
        public HSLFilter(int hue, double saturation, double lightness)
        {
            Hue = hue;
            Saturation = saturation;
            Lightness = lightness;
        }
        /// <summary>
        /// Initializes the HSL filter.
        /// </summary>
        public HSLFilter()
        {
            new HSLFilter(0, 0, 0);
        }
        /// <summary>
        /// Hue [0, 359].
        /// </summary>
        public int Hue
        {
            get
            {
                return this.hue;
            }
            set
            {
                this.hue = value;
            }
        }
        /// <summary>
        /// Saturation [-1, 1].
        /// </summary>
        public double Saturation
        {
            get
            {
                return this.saturation;
            }
            set
            {
                this.saturation = value;
            }
        }
        /// <summary>
        /// Lightness [-1, 1].
        /// </summary>
        public double Lightness
        {
            get
            {
                return this.lightness;
            }
            set
            {
                this.lightness = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                HSL hsl; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    hsl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    hsl.Hue = hue;
                    hsl.Saturation += saturation;
                    hsl.Lightness += lightness;
                    rgb = hsl.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Defines the YCbCr filter.
    /// </summary>
    public class YCbCrFilter : IBitmapFilter
    {
        #region Private data
        private double y;
        private double cb;
        private double cr;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the YCbCr filter.
        /// </summary>
        /// <param name="y">Y [-1, 1]</param>
        /// <param name="cb">Cb [-1, 1]</param>
        /// <param name="cr">Cr [-1, 1]</param>
        public YCbCrFilter(double y, double cb, double cr)
        {
            Y = y;
            Cb = cb;
            Cr = cr;
        }
        /// <summary>
        /// Initializes the YCbCr filter.
        /// </summary>
        public YCbCrFilter()
        {
            new YCbCrFilter(0, 0, 0);
        }
        /// <summary>
        /// Y [-1, 1].
        /// </summary>
        public double Y
        {
            get
            {
                return this.y;
            }
            set
            {
                this.y = value;
            }
        }
        /// <summary>
        /// Cb [-1, 1].
        /// </summary>
        public double Cb
        {
            get
            {
                return this.cb;
            }
            set
            {
                this.cb = value;
            }
        }
        /// <summary>
        /// Cr [-1, 1].
        /// </summary>
        public double Cr
        {
            get
            {
                return this.cr;
            }
            set
            {
                this.cr = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                YCbCr ycbcr; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    ycbcr = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    ycbcr.Y += y;
                    ycbcr.Cb += cb;
                    ycbcr.Cr += cr;
                    rgb = ycbcr.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the CMYK filter.
    /// </summary>
    public class CMYKFilter : IBitmapFilter
    {
        #region Private data
        private double cyan;
        private double magenta;
        private double yellow;
        private double keycolor;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the CMYK filter.
        /// </summary>
        /// <param name="cyan">Cyan [-1, 1]</param>
        /// <param name="magenta">Magenta [-1, 1]</param>
        /// <param name="yellow">Yellow [-1, 1]</param>
        /// <param name="keycolor">Keycolor [-1, 1]</param>
        public CMYKFilter(double cyan, double magenta, double yellow, double keycolor)
        {
            Cyan = cyan;
            Magenta = magenta;
            Yellow = yellow;
            Keycolor = keycolor;
        }
        /// <summary>
        /// Initializes the CMYK filter.
        /// </summary>
        public CMYKFilter()
        {
            new CMYKFilter(0, 0, 0, 0);
        }
        /// <summary>
        /// Cyan [-1, 1].
        /// </summary>
        public double Cyan
        {
            get
            {
                return this.cyan;
            }
            set
            {
                this.cyan = value;
            }
        }
        /// <summary>
        /// Magenta [-1, 1].
        /// </summary>
        public double Magenta
        {
            get
            {
                return this.magenta;
            }
            set
            {
                this.magenta = value;
            }
        }
        /// <summary>
        /// Yellow [-1, 1].
        /// </summary>
        public double Yellow
        {
            get
            {
                return this.yellow;
            }
            set
            {
                this.yellow = value;
            }
        }
        /// <summary>
        /// Keycolor [-1, 1].
        /// </summary>
        public double Keycolor
        {
            get
            {
                return this.keycolor;
            }
            set
            {
                this.keycolor = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            Parallel.For(0, height, j =>
            {
                CMYK cmyk; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    cmyk = CMYK.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    cmyk.Cyan += cyan;
                    cmyk.Magenta += magenta;
                    cmyk.Yellow += yellow;
                    cmyk.Keycolor += keycolor;
                    rgb = cmyk.ToRGB;

                    p[k + 0] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the RGB filter.
    /// </summary>
    public class RGBFilter : IBitmapFilter
    {
        #region Private data
        private int red;
        private int green;
        private int blue;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the RGB filter.
        /// </summary>
        /// <param name="red">Red [-255, 255]</param>
        /// <param name="green">Green [-255, 255]</param>
        /// <param name="blue">Blue [-255, 255]</param>
        public RGBFilter(int red, int green, int blue)
        {
            Red = red;
            Green = green;
            Blue = blue;
        }
        /// <summary>
        /// Initializes the RGB filter.
        /// </summary>
        public RGBFilter()
        {
            new RGBFilter(0, 0, 0);
        }
        /// <summary>
        /// Red [-255, 255].
        /// </summary>
        public int Red
        {
            get
            {
                return this.red;
            }
            set
            {
                this.red = value;
            }
        }
        /// <summary>
        /// Green [-255, 255].
        /// </summary>
        public int Green
        {
            get
            {
                return this.green;
            }
            set
            {
                this.green = value;
            }
        }
        /// <summary>
        /// Blue [-255, 255].
        /// </summary>
        public int Blue
        {
            get
            {
                return this.blue;
            }
            set
            {
                this.blue = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap</param>
        public unsafe void Apply(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(p[2] + red);
                    p[1] = Maths.Byte(p[1] + green);
                    p[0] = Maths.Byte(p[0] + blue);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    #endregion

    #region Histogram equalizations
    /// <summary>
    /// Defines the global histogram equalization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.cromwell-intl.com/3d/histogram/
    /// </remarks>
    /// </summary>
    public class HistogramEqualization : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Initializes the global histogram equalization filter.
        /// </summary>
        public HistogramEqualization() { }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData)
        {
            int[] H = Statistics.Histogram(bmData);
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            double[] table = Statistics.Equalize(H);
            double length = table.Length - 1;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    p[2] = Maths.Byte(table[p[2]] * length);
                    p[1] = Maths.Byte(table[p[1]] * length);
                    p[0] = Maths.Byte(table[p[0]] * length);
                }
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the local histogram equalization filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://angeljohnsy.blogspot.com/2011/06/local-histogram-equalization.html
    /// </remarks>
    /// </summary>
    public class LocalHistogramEqualization : IBitmapFilter2
    {
        #region Private data
        private int l0;
        private int l1;
        private int rw;
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        public LocalHistogramEqualization(int radius = 10)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        public LocalHistogramEqualization(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Initializes the local histogram equalization filter.
        /// </summary>
        /// <param name="size">Radius</param>
        public LocalHistogramEqualization(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return new SizeInt(l0, l1);
            }
            set
            {
                Data(value);
            }
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="size">Radius</param>
        private void Data(SizeInt size)
        {
            this.l0 = size.Width;
            this.l1 = size.Height;
            this.rw = l0 >> 1;
            this.rh = l1 >> 1;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            #endregion

            Parallel.For(0, height, y =>
            {
                int[] H = new int[256];
                int[] cdf = new int[256];
                int n = l0 * l1;
                int brightness;
                double dn = 255.0 / n;

                int x, i, j, ir, jr, yr, xr, irstride;
                int ystride, v;
                byte* p;

                yr = y - rw;
                ystride = y * stride;

                for (x = 0; x < width; x++)
                {
                    xr = x - rh;
                    v = ystride + x * 4;

                    #region Convolution filtering
                    for (i = 0; i < l0; i++)
                    {
                        ir = yr + i;
                        if (ir < 0) continue; if (ir >= height) break;
                        irstride = ir * stride;

                        for (j = 0; j < l1; j++)
                        {
                            jr = xr + j;

                            if (jr < 0) continue; if (jr >= width) break;

                            p = &src[irstride + jr * 4];
                            brightness = (p[2] + p[1] + p[0]) / 3;
                            H[brightness]++;
                        }
                    }
                    #endregion

                    #region Density function
                    cdf[0] = H[0];

                    for (i = 1; i < 256; i++)
                    {
                        cdf[i] = H[i] + cdf[i - 1];
                    }
                    #endregion

                    #region Recording pixel
                    dst[v + 2] = Maths.Byte(cdf[src[v + 2]] * dn);
                    dst[v + 1] = Maths.Byte(cdf[src[v + 1]] * dn);
                    dst[v    ] = Maths.Byte(cdf[src[v]] * dn);
                    #endregion

                    #region Clear data
                    Array.Clear(H, 0, 256);
                    Array.Clear(cdf, 0, 256);
                    #endregion
                }
            }
            );

            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="Src">Bitmap</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Src = (Bitmap)Data.Clone();
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
            Src.Dispose();
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the global histogram stretch filter.
    /// </summary>
    public class HistogramStretch : Correction, IBitmapFilter
    {
        #region Private data
        private RangeDouble range;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the global histogram stretch filter.
        /// </summary>
        /// <param name="min">Minimum intensity [0, 1]</param>
        /// <param name="max">Maximum intensity [0, 1]</param>
        /// <param name="space">Color space</param>
        public HistogramStretch(double min, double max, Space space)
        {
            Range = new RangeDouble(min, max);
            Space = space;
        }
        /// <summary>
        /// Initializes the global histogram stretch filter.
        /// </summary>
        /// <param name="range">Intensity range</param>
        /// <param name="space">Color space</param>
        public HistogramStretch(RangeDouble range, Space space)
        {
            Range = range;
            Space = space;
        }
        /// <summary>
        /// Gets or sets the intensity range.
        /// </summary>
        public RangeDouble Range
        {
            get
            {
                return this.range;
            }
            set
            {
                this.range = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Implements filter rebuilding.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Equalize(range.Min, range.Max, 256);
        }
        #endregion
    }
    /// <summary>
    /// Defines the local histogram stretch filter.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.academia.edu/7629047/Image_enhancement_by_local_histogram_stretching
    /// </remarks>
    /// </summary>
    public class LocalHistogramStretch : IBitmapFilter
    {
        #region Private data
        private BoxBlur gb = new BoxBlur();
        private Erosion er = new Erosion();
        private Dilatation di = new Dilatation();
        private Space space;
        private double contrast;
        private bool smoothing;
        #endregion

        #region Filter componetns
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="radius">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(int radius, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = new SizeInt(radius, radius);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="width">Filter width</param>
        /// <param name="height">Filter height</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(int width, int height, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = new SizeInt(width, height);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Initializes the local histogram stretch filter.
        /// </summary>
        /// <param name="size">Radius</param>
        /// <param name="space">Color space</param>
        /// <param name="contrast">Contrast [0, 1]</param>
        /// <param name="smoothing">Smoothing</param>
        public LocalHistogramStretch(SizeInt size, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = size;
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Gets or sets the filter size.
        /// </summary>
        public SizeInt Size
        {
            get
            {
                return er.Size;
            }
            set
            {
                er.Size = value; // local min  filter -> r,
                di.Size = value; // local max  filter -> r,
                gb.Size =        // local mean filter -> 2r.
                    new SizeInt( // we got only even values for radius.
                        value.Width * 2,
                        value.Height * 2);
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
        /// Gets or sets the contrast value [0, 1].
        /// </summary>
        public double Contrast
        {
            get
            {
                return this.contrast;
            }
            set
            {
                this.contrast = Maths.Double(value);
            }
        }
        /// <summary>
        /// Smoothing.
        /// </summary>
        public bool Smoothing
        {
            get
            {
                return this.smoothing;
            }
            set
            {
                this.smoothing = value;
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            Bitmap Max = (Bitmap)BitmapConverter.Bitmap(bmData).Clone();
            Bitmap Min = (Bitmap)Max.Clone();

            di.Apply(Max); er.Apply(Min);

            BitmapData bmMax = BitmapConverter.Lock32bpp(Max);
            BitmapData bmMin = BitmapConverter.Lock32bpp(Min);

            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            Apply(bmData, bmMax, bmMin);

            BitmapConverter.Unlock(Max, bmMax);
            BitmapConverter.Unlock(Min, bmMin);

            Max.Dispose(); Min.Dispose();
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            Bitmap Max = (Bitmap)Data.Clone();
            Bitmap Min = (Bitmap)Data.Clone();

            di.Apply(Max); er.Apply(Min);

            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmMax = BitmapConverter.Lock32bpp(Max);
            BitmapData bmMin = BitmapConverter.Lock32bpp(Min);

            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            Apply(bmData, bmMax, bmMin);

            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Max, bmMax);
            BitmapConverter.Unlock(Min, bmMin);

            Max.Dispose(); Min.Dispose();
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void Apply(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData, bmMax, bmMin);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData, bmMax, bmMin);
                    break;
                default:
                    ApplyGrayscale(bmData, bmMax, bmMin);
                    break;
            }
            return;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyRGB(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            double required = 1.0 - this.contrast;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, q, v, jstride = j * stride;
                double mag, max, min;
                double num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    // Local function:
                    for (q = 0; q < 3; q++)
                    {
                        v = k + q;

                        mag = p[v] / 255.0;
                        max = pMax[v] / 255.0;
                        min = pMin[v] / 255.0;

                        num1 = max - min;

                        if (num1 < required)
                        {
                            num2 = min + (required - num1) * min / (num1 - 1f);
                            min = Maths.Double(num2);
                            max = Maths.Double(num2 + required);
                        }

                        num1 = max - min;
                        num3 = mag - min;

                        if (num1 > 0)
                        {
                            p[v] = Maths.Byte(255 * num3 / num1);
                        }
                    }
                    // end local function.
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyHSB(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            double required = 1.0 - this.contrast;

            Parallel.For(0, height, j =>
            {
                HSB imag; HSB imax; HSB imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                double mag, max, min;
                double num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = HSB.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSB.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSB.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Brightness;
                    max = imax.Brightness;
                    min = imin.Brightness;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Brightness = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyHSL(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            double required = 1.0 - this.contrast;

            Parallel.For(0, height, j =>
            {
                HSL imag; HSL imax; HSL imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                double mag, max, min;
                double num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = HSL.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSL.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSL.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Lightness;
                    max = imax.Lightness;
                    min = imin.Lightness;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Lightness = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyYCbCr(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            double required = 1.0 - this.contrast;

            Parallel.For(0, height, j =>
            {
                YCbCr imag; YCbCr imax; YCbCr imin; RGB rgb;
                int i, k, k1, k2, jstride = j * stride;
                double mag, max, min;
                double num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    imag = YCbCr.FromRGB(p[k2], p[k1], p[k]);
                    imax = YCbCr.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = YCbCr.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    mag = imag.Y;
                    max = imax.Y;
                    min = imin.Y;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        imag.Y = num3 / num1;
                        rgb = imag.ToRGB;
                        p[k2] = rgb.Red; p[k1] = rgb.Green; p[k] = rgb.Blue;
                    }
                }
            }
            );
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmMax">Bitmap data</param>
        /// <param name="bmMin">Bitmap data</param> 
        private unsafe void ApplyGrayscale(BitmapData bmData, BitmapData bmMax, BitmapData bmMin)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pMax = (byte*)bmMax.Scan0.ToPointer();
            byte* pMin = (byte*)bmMin.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;

            double required = 1.0 - this.contrast;

            Parallel.For(0, height, j =>
            {
                int i, k, k1, k2, v, jstride = j * stride;
                double mag, max, min;
                double num1, num2, num3;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4; k1 = k + 1; k2 = k + 2;

                    // Local function:
                    v = k;

                    mag = RGB.Average(p[v], p[v + 1], p[v + 2]) / 255.0;
                    max = RGB.Average(pMax[v], pMax[v + 1], pMax[v + 2]) / 255.0;
                    min = RGB.Average(pMin[v], pMin[v + 1], pMin[v + 2]) / 255.0;

                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    num1 = max - min;
                    num3 = mag - min;

                    if (num1 > 0)
                    {
                        p[v] = p[v + 1] = p[v + 2] = Maths.Byte(255 * num3 / num1);
                    }
                    // end local function.
                }
            }
            );
        }
        #endregion
    }
    #endregion

    #region Special filters
    /// <summary>
    /// Defines the bitmap filter.
    /// </summary>
    public class BitmapFilter : IBitmapFilter
    {
        #region Private data
        private IFilter filter;
        private Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bitmap filter.
        /// </summary>
        /// <param name="filter">Filter</param>
        /// <param name="space">Color space</param>
        public BitmapFilter(IFilter filter, Space space = Space.RGB)
        {
            this.filter = filter;
            this.space = space;
        }
        /// <summary>
        /// Gets or sets the filter.
        /// </summary>
        public IFilter Filter
        {
            get
            {
                return this.filter;
            }
            set
            {
                this.filter = value;
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
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    ApplyHSB(bmData);
                    break;
                case Imaging.Space.HSL:
                    ApplyHSL(bmData);
                    break;
                case Imaging.Space.YCbCr:
                    ApplyYCbCr(bmData);
                    break;
                case Imaging.Space.RGB:
                    ApplyRGB(bmData);
                    break;
                default:
                    ApplyGrayscale(bmData);
                    break;
            }
            return;
        }
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyRGB(BitmapData bmData)
        {
            double[][,] rgb = BitmapConverter.ToRGB(bmData, true);

            this.filter.Apply(rgb[0]);
            this.filter.Apply(rgb[1]);
            this.filter.Apply(rgb[2]);

            BitmapConverter.FromRGB(rgb, bmData);
            return;
        }
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSB(BitmapData bmData)
        {
            double[][,] hsb = BitmapConverter.ToHSB(bmData, true);

            this.filter.Apply(hsb[0]);
            this.filter.Apply(hsb[1]);
            this.filter.Apply(hsb[2]);

            BitmapConverter.FromHSB(hsb, bmData);
            return;
        }
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyHSL(BitmapData bmData)
        {
            double[][,] hsl = BitmapConverter.ToHSL(bmData, true);

            this.filter.Apply(hsl[0]);
            this.filter.Apply(hsl[1]);
            this.filter.Apply(hsl[2]);

            BitmapConverter.FromHSL(hsl, bmData);
            return;
        }
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyYCbCr(BitmapData bmData)
        {
            double[][,] ycbcr = BitmapConverter.ToYCbCr(bmData, true);

            this.filter.Apply(ycbcr[0]);
            this.filter.Apply(ycbcr[1]);
            this.filter.Apply(ycbcr[2]);

            BitmapConverter.FromHSL(ycbcr, bmData);
            return;
        }
        /// <summary>
        /// Appy filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        private unsafe void ApplyGrayscale(BitmapData bmData)
        {
            double[,] y = BitmapConverter.ToGrayscale(bmData);
            this.filter.Apply(y);
            BitmapConverter.FromGrayscale(y, bmData);
            return;
        }
        #endregion
    }
    /// <summary>
    /// Defines the bitmap blending filter.
    /// <remarks>
    /// Filter can be used to implement HDR effects.
    /// More information can be found on the website:
    /// https://web.stanford.edu/class/cs231m/lectures/lecture-5-stitching-blending.pdf
    /// (стр. 65-75)
    /// </remarks>
    /// </summary>
    public class BitmapBlender
    {
        #region Private data
        private IBlendFilter filter;
        private Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Initializes the bitmap blending filter.
        /// </summary>
        /// <param name="filter">Blending filter</param>
        /// <param name="space">Color space</param>
        public BitmapBlender(IBlendFilter filter, Space space = Space.RGB)
        {
            this.filter = filter;
            this.space = space;
        }
        /// <summary>
        /// Gets or sets the blending filter.
        /// </summary>
        public IBlendFilter Filter
        {
            get
            {
                return this.filter;
            }
            set
            {
                this.filter = value;
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
        /// <param name="images">Bitmap array</param>
        /// <returns>Bitmap</returns>
        public Bitmap Apply(params Bitmap[] images)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    return ApplyHSB(images);

                case Imaging.Space.HSL:
                    return ApplyHSL(images);

                case Imaging.Space.YCbCr:
                    return ApplyYCbCr(images);

                case Imaging.Space.RGB:
                    return ApplyRGB(images);

                default:
                    return ApplyGrayscale(images);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="images">Bitmap array</param>
        /// <returns>Bitmap</returns>
        public Bitmap Apply(params BitmapData[] images)
        {
            // filter
            switch (space)
            {
                case Imaging.Space.HSB:
                    return ApplyHSB(images);

                case Imaging.Space.HSL:
                    return ApplyHSL(images);

                case Imaging.Space.YCbCr:
                    return ApplyYCbCr(images);

                case Imaging.Space.RGB:
                    return ApplyRGB(images);

                default:
                    return ApplyGrayscale(images);
            }
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyRGB(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] r = new double[length][,];
            double[][,] g = new double[length][,];
            double[][,] b = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToRGB(Data[i], false);

                r[i] = t[2];
                g[i] = t[1];
                b[i] = t[0];
            }

            return BitmapConverter.FromRGB(new double[][,] { this.filter.Apply(b), this.filter.Apply(g), this.filter.Apply(r) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyRGB(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] r = new double[length][,];
            double[][,] g = new double[length][,];
            double[][,] b = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToRGB(bmData[i], false);

                r[i] = t[2];
                g[i] = t[1];
                b[i] = t[0];
            }

            return BitmapConverter.FromRGB(new double[][,] { this.filter.Apply(b), this.filter.Apply(g), this.filter.Apply(r) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyGrayscale(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] r = new double[length][,];

            for (int i = 0; i < length; i++)
            {
                r[i] = BitmapConverter.ToGrayscale(Data[i]);
            }

            return BitmapConverter.FromGrayscale(this.filter.Apply(r));
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyGrayscale(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] r = new double[length][,];

            for (int i = 0; i < length; i++)
            {
                r[i] = BitmapConverter.ToGrayscale(bmData[i]);
            }

            return BitmapConverter.FromGrayscale(this.filter.Apply(r));
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyYCbCr(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] yy = new double[length][,];
            double[][,] cb = new double[length][,];
            double[][,] cr = new double[length][,];
            double[][,] tt;

            for (int i = 0; i < length; i++)
            {
                tt = BitmapConverter.ToYCbCr(Data[i], false);

                yy[i] = tt[0];
                cb[i] = tt[1];
                cr[i] = tt[2];
            }

            return BitmapConverter.FromYCbCr(new double[][,] { this.filter.Apply(yy), BoxFilterOptions.boxf(cb), BoxFilterOptions.boxf(cr) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyYCbCr(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] yy = new double[length][,];
            double[][,] cb = new double[length][,];
            double[][,] cr = new double[length][,];
            double[][,] tt;

            for (int i = 0; i < length; i++)
            {
                tt = BitmapConverter.ToYCbCr(bmData[i], false);

                yy[i] = tt[0];
                cb[i] = tt[1];
                cr[i] = tt[2];
            }

            return BitmapConverter.FromYCbCr(new double[][,] { this.filter.Apply(yy), BoxFilterOptions.boxf(cb), BoxFilterOptions.boxf(cr) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyHSB(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] h = new double[length][,];
            double[][,] s = new double[length][,];
            double[][,] b = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToHSB(Data[i], false);

                h[i] = t[0];
                s[i] = t[1];
                b[i] = t[2];
            }

            return BitmapConverter.FromHSB(new double[][,] { BoxFilterOptions.boxf(h), BoxFilterOptions.boxf(s), this.filter.Apply(b) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyHSB(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] h = new double[length][,];
            double[][,] s = new double[length][,];
            double[][,] b = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToHSB(bmData[i], false);

                h[i] = t[0];
                s[i] = t[1];
                b[i] = t[2];
            }

            return BitmapConverter.FromHSB(new double[][,] { BoxFilterOptions.boxf(h), BoxFilterOptions.boxf(s), this.filter.Apply(b) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyHSL(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] h = new double[length][,];
            double[][,] s = new double[length][,];
            double[][,] l = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToHSL(Data[i], false);

                h[i] = t[0];
                s[i] = t[1];
                l[i] = t[2];
            }

            return BitmapConverter.FromHSL(new double[][,] { BoxFilterOptions.boxf(h), BoxFilterOptions.boxf(s), this.filter.Apply(l) });
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap array</param>
        /// <returns>Bitmap</returns>
        private Bitmap ApplyHSL(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] h = new double[length][,];
            double[][,] s = new double[length][,];
            double[][,] l = new double[length][,];
            double[][,] t;

            for (int i = 0; i < length; i++)
            {
                t = BitmapConverter.ToHSL(bmData[i], false);

                h[i] = t[0];
                s[i] = t[1];
                l[i] = t[2];
            }

            return BitmapConverter.FromHSL(new double[][,] { BoxFilterOptions.boxf(h), BoxFilterOptions.boxf(s), this.filter.Apply(l) });
        }
        #endregion
    }
    #endregion

    #region Bitmap conversions
    /// <summary>
    /// Uses to work with bitmaps.
    /// </summary>
    public static class BitmapConverter
    {
        #region Bitmap convert components
        /// <summary>
        /// Converts Bitmap to icon file.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="size">Size</param>
        /// <returns>Icon</returns>
        public static Icon ToIco(this Bitmap b, int size)
        {
            Bitmap bmp = new Bitmap(b, new Size(size, size));
            MemoryStream pngstream = new MemoryStream();
            MemoryStream icostream = new MemoryStream();

            byte[] pngicon = new byte[] { 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 24, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            byte[] png;

            bmp.Save(pngstream, System.Drawing.Imaging.ImageFormat.Png);
            pngstream.Position = 0;
            png = pngstream.ToArray();

            if (size >= 256) size = 0;
            pngicon[6] = (byte)size;
            pngicon[7] = (byte)size;
            pngicon[14] = (byte)(png.Length & 255);
            pngicon[15] = (byte)(png.Length / 256);
            pngicon[18] = (byte)(pngicon.Length);

            icostream.Write(pngicon, 0, pngicon.Length);
            icostream.Write(png, 0, png.Length);
            icostream.Position = 0;
            bmp.Dispose();
            pngstream.Dispose();

            return new Icon(icostream);
        }
        /// <summary>
        /// Converts Bitmap to JPEG format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToJpeg(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Jpeg);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to BMP format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToBmp(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Bmp);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to GIF format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToGif(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Gif);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to PNG format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToPng(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Png);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Converts Bitmap to TIFF format.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap</returns>
        public static Bitmap ToTiff(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Tiff);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Gets the Bitmap from the BitmapData.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Bitmap(this BitmapData bmData)
        {
            return new Bitmap(bmData.Width, bmData.Height, bmData.Stride, bmData.PixelFormat, bmData.Scan0);
        }
        /// <summary>
        /// Converts Bitmap to a specific format
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="pixelformat">Pixel format</param>
        /// <returns>Bitmap</returns>
        public static Bitmap Bitmap(this Bitmap b, PixelFormat pixelformat)
        {
            return b.Clone(new Rectangle(0, 0, b.Width, b.Height), pixelformat);
        }
        #endregion

        #region BitmapData voids
        /// <summary>
        /// Blocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap data</returns>
        public static BitmapData Lock32bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format32bppArgb);
        }
        /// <summary>
        /// Blocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <returns>Bitmap data</returns>
        public static BitmapData Lock8bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format8bppIndexed);
        }
        /// <summary>
        /// Unblocks Bitmap in system memory.
        /// </summary>
        /// <param name="b">Bitmap</param>
        /// <param name="bmData">Bitmap data</param>
        public static void Unlock(this Bitmap b, BitmapData bmData)
        {
            b.UnlockBits(bmData);
            return;
        }
        #endregion

        #region RGB
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public static double[][,] ToRGB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToRGB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an RGB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>RGBA structure array</returns>
        public unsafe static double[][,] ToRGB(BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // color
                        x[j, i] = p[k + 0] / 255.0;
                        y[j, i] = p[k + 1] / 255.0;
                        z[j, i] = p[k + 2] / 255.0;
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    // color
                    x[j, i] = p[k + 0] / 255.0;
                    y[j, i] = p[k + 1] / 255.0;
                    z[j, i] = p[k + 2] / 255.0;
                }
            });

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromRGB(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = (byte)255;
                    }
                });
            }


            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromRGB(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        // recording model:
                        p[k + 0] = Maths.Byte(x[j, i] * 255.0);
                        p[k + 1] = Maths.Byte(y[j, i] * 255.0);
                        p[k + 2] = Maths.Byte(z[j, i] * 255.0);
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an RGB structure to a color image.
        /// </summary>
        /// <param name="array">RGBA structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromRGB(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromRGB(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region HSB
        /// <summary>
        /// Converts a Bitmap to an HSB structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSB structure array</returns>
        public static double[][,] ToHSB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSB structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSB structure array</returns>
        public unsafe static double[][,] ToHSB(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

                Parallel.For(0, height, j =>
                {
                    HSB mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Hue;
                        y[j, i] = mdl.Saturation;
                        z[j, i] = mdl.Brightness;
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                HSB mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = HSB.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Hue;
                    y[j, i] = mdl.Saturation;
                    z[j, i] = mdl.Brightness;
                }
            });

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSB(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSB(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSB(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an HSB structure to a color image.
        /// </summary>
        /// <param name="array">HSB structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromHSB(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromHSB(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region HSL
        /// <summary>
        /// Converts a Bitmap to an HSL structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSL structure array</returns>
        public static double[][,] ToHSL(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSL(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an HSL structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>HSL structure array</returns>
        public unsafe static double[][,] ToHSL(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

                Parallel.For(0, height, j =>
                {
                    HSL mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Hue;
                        y[j, i] = mdl.Saturation;
                        z[j, i] = mdl.Lightness;
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                HSL mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = HSL.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Hue;
                    y[j, i] = mdl.Saturation;
                    z[j, i] = mdl.Lightness;
                }
            });

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromHSL(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromHSL(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new HSL(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an HSL structure to a color image.
        /// </summary>
        /// <param name="array">HSL structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromHSL(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromHSL(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region YCbCr
        /// <summary>
        /// Converts a Bitmap to an YCbCr structure with or without alpha-channel.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>YCbCr structure array</returns>
        public static double[][,] ToYCbCr(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToYCbCr(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts a Bitmap to an YCbCr structure with or without alpha-channel.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="alpha">Alpha-channel</param>
        /// <returns>YCbCr structure array</returns>
        public unsafe static double[][,] ToYCbCr(this BitmapData bmData, bool alpha = false)
        {
            // params
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // matrices
            double[,] x = new double[height, width];
            double[,] y = new double[height, width];
            double[,] z = new double[height, width];

            // with alpha channel
            if (alpha)
            {
                double[,] a = new double[height, width];

                Parallel.For(0, height, j =>
                {
                    YCbCr mdl;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        mdl = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                        x[j, i] = mdl.Y;
                        y[j, i] = mdl.Cb;
                        z[j, i] = mdl.Cr;
                        a[j, i] = p[k + 3] / 255.0;
                    }
                });

                return new double[][,] { x, y, z, a };
            }

            // without alpha channel
            Parallel.For(0, height, j =>
            {
                YCbCr mdl;

                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * 4;

                    mdl = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k + 0]);
                    x[j, i] = mdl.Y;
                    y[j, i] = mdl.Cb;
                    z[j, i] = mdl.Cr;
                }
            });

            return new double[][,] { x, y, z };
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromYCbCr(this double[][,] array)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromYCbCr(this double[][,] array, BitmapData bmData)
        {
            // matrices
            double[,] x = array[0];
            double[,] y = array[1];
            double[,] z = array[2];

            // params
            int width = x.GetLength(1), height = x.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            // alpha
            bool alpha = array.Length == 4;

            if (alpha)
            {
                double[,] a = array[3];

                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = Maths.Byte(a[j, i] * 255.0);
                    }
                });
            }
            else
            {
                Parallel.For(0, height, j =>
                {
                    RGB rgb;

                    int i, k, jstride = j * stride;

                    for (i = 0; i < width; i++)
                    {
                        // shift:
                        k = jstride + i * 4;

                        rgb = new YCbCr(x[j, i], y[j, i], z[j, i]).ToRGB;

                        // recording model:
                        p[k + 0] = rgb.Blue;
                        p[k + 1] = rgb.Green;
                        p[k + 2] = rgb.Red;
                        p[k + 3] = (byte)255;
                    }
                });
            }

            return;
        }
        /// <summary>
        /// Converts an YCbCr structure to a color image.
        /// </summary>
        /// <param name="array">YCbCr structure array</param>
        /// <param name="Data">Bitmap</param>
        public static void FromYCbCr(this double[][,] array, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromYCbCr(array, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion

        #region Grayscale
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Matrix</returns>
        public static double[,] ToGrayscale(this Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[,] rgb = ToGrayscale(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Converts Bitmap to averaged channel value matrix.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Matrix</returns>
        public unsafe static double[,] ToGrayscale(this BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] rgb = new double[height, width];
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    rgb[j, i] = RGB.Average(p[k + 2], p[k + 1], p[k]) / 255.0;
                }
            });

            return rgb;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <returns>Bitmap</returns>
        public unsafe static Bitmap FromGrayscale(this double[,] m)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0);
                    p[k + 3] = 255;
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="bmData">Bitmap data</param>
        public unsafe static void FromGrayscale(this double[,] m, BitmapData bmData)
        {
            int width = m.GetLength(1), height = m.GetLength(0);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = p[k + 1] = p[k] = Maths.Byte(m[j, i] * 255.0);
                    p[k + 3] = 255;
                }
            });

            return;
        }
        /// <summary>
        /// Converts a matrix of channel values to a monochrome Bitmap.
        /// </summary>
        /// <param name="m">Matrix</param>
        /// <param name="Data">Bitmap</param>
        public static void FromGrayscale(this double[,] m, Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            FromGrayscale(m, bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        #endregion
    }
    #endregion

    #region Hough transforms
    /// <summary>
    /// Defines the Hough line.
    /// </summary>
    public struct HoughLine : IComparable
    {
        #region Struct components
        /// <summary>
        /// Slope of the line [0, 180).
        /// <remarks>
        /// It is the angle between the polar axis and the radius of the line.
        /// </remarks>
        /// </summary>
        public readonly double Theta;
        /// <summary>
        /// Line distance from the center of the image (-inf, +inf).
        /// <remarks>
        /// A negative radius line means that the line is at the bottom of the polar coordinate system. Therefore
        /// the angle θ should be increased by 180 degrees, and Radius should be positive.
        /// </remarks>
        /// </summary>
        public readonly short Radius;
        /// <summary>
        /// Absolute line intensity (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Relative line intensity (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Initializes the Hough line.
        /// </summary>
        /// <param name="theta">Slope of the line [0, 180)</param>
        /// <param name="radius">Radius (-inf, +inf)</param>
        /// <param name="intensity">Absolute line intensity (0, +inf)</param>
        /// <param name="relativeIntensity">Relative line intensity (0, 1]</param>
        public HoughLine(double theta, short radius, short intensity, double relativeIntensity)
        {
            Theta = theta;
            Radius = radius;
            Intensity = intensity;
            RelativeIntensity = relativeIntensity;
        }
        /// <summary>
        /// Compares object to another instance of this class.
        /// </summary>
        /// <param name="value">Object</param>
        /// <returns>Integer number</returns>
        public int CompareTo(object value)
        {
            return (-Intensity.CompareTo(((HoughLine)value).Intensity));
        }
        #endregion
    }
    /// <summary>
    /// Defines the Hough circle.
    /// </summary>
    public struct HoughCircle : IComparable
    {
        #region Struct components
        /// <summary>
        /// Coordinate X.
        /// </summary>
        public readonly int X;
        /// <summary>
        /// Coordinate Y.
        /// </summary>
        public readonly int Y;
        /// <summary>
        /// Radius.
        /// </summary>
        public readonly int Radius;
        /// <summary>
        /// Absolute line intensity (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Relative line intensity (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Initializes the Hough circle.
        /// </summary>
        /// <param name="x">Coordinate X</param>
        /// <param name="y">Coordinate Y</param>
        /// <param name="radius">Radius</param>
        /// <param name="intensity">Absolute line intensity (0, +inf)</param>
        /// <param name="relativeIntensity">Relative line intensity (0, 1]</param>
        public HoughCircle(int x, int y, int radius, short intensity, double relativeIntensity)
        {
            X = x;
            Y = y;
            Radius = radius;
            Intensity = intensity;
            RelativeIntensity = relativeIntensity;
        }
        /// <summary>
        /// Compares object to another instance of this class.
        /// </summary>
        /// <param name="value">Object</param>
        /// <returns>Integer number</returns>
        public int CompareTo(object value)
        {
            return (-Intensity.CompareTo(((HoughCircle)value).Intensity));
        }
        #endregion
    }
    /// <summary>
    /// Defines the Hough line transform filter.
    /// </summary>
    public class HoughLineTransform : IBitmapFilter
    {
        #region Private data
        private int stepsPerDegree;
        private int houghHeight;
        private double thetaStep;

        private double[] sinMap;
        private double[] cosMap;
        
        private short[,] houghMap;
        private short maxMapIntensity = 0;

        private int localPeakRadius = 4;
        private short minLineIntensity = 10;
        private ArrayList lines = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hough line transform filter.
        /// </summary>
        public HoughLineTransform()
        {
            StepsPerDegree = 1;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Gets or sets steps per degree.
        /// </summary>
        public int StepsPerDegree
        {
            get { return stepsPerDegree; }
            set
            {
                stepsPerDegree = Math.Max(1, Math.Min(10, value));
                houghHeight = 180 * stepsPerDegree;
                thetaStep = Math.PI / houghHeight;

                // precalculate Sine and Cosine values
                sinMap = new double[houghHeight];
                cosMap = new double[houghHeight];

                for (int i = 0; i < houghHeight; i++)
                {
                    sinMap[i] = Math.Sin(i * thetaStep);
                    cosMap[i] = Math.Cos(i * thetaStep);
                }
            }
        }
        /// <summary>
        /// Gets the minimum line intensity.
        /// </summary>
        public short MinLineIntensity
        {
            get { return minLineIntensity; }
            set { minLineIntensity = value; }
        }
        /// <summary>
        /// Gets the maximum line intensity.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
        }
        /// <summary>
        /// Gets or sets the radius search for the local peak value.
        /// </summary>
        public int LocalPeakRadius
        {
            get { return localPeakRadius; }
            set { localPeakRadius = Math.Max(1, Math.Min(10, value)); }
        }
        /// <summary>
        /// Gets the count of lines found.
        /// </summary>
        public int LinesCount
        {
            get { return lines.Count; }
        }
        /// <summary>
        /// Returns an array of lines with absolute intensity.
        /// </summary>
        /// <param name="count">Count</param>
        /// <returns>Array</returns>
        public HoughLine[] GetMostIntensiveLines(int count)
        {
            // lines count
            int n = Math.Min(count, lines.Count);

            // result array
            HoughLine[] dst = new HoughLine[n];
            lines.CopyTo(0, dst, 0, n);

            return dst;
        }
        /// <summary>
        /// Returns an array of lines with relative intensity.
        /// </summary>
        /// <param name="minRelativeIntensity">Minimum relative intensity</param>
        /// <returns>Array</returns>
        public HoughLine[] GetLinesByRelativeIntensity(double minRelativeIntensity)
        {
            int count = 0, n = lines.Count;

            while ((count < n) && (((HoughLine)lines[count]).RelativeIntensity >= minRelativeIntensity))
                count++;

            return GetMostIntensiveLines(count);
        }
        #endregion

        #region Image processing
        /// <summary>
        /// Returns the Hough synogram.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap ToBitmap()
        {
            // check if Hough transformation was made already
            if (houghMap == null)
            {
                throw new ApplicationException("Hough conversion was not made");
            }

            int height = houghMap.GetLength(0);
            int width = houghMap.GetLength(1);

            // create new image
            Bitmap image = new Bitmap(width, height, PixelFormat.Format8bppIndexed);

            // get palette
            ColorPalette cp = image.Palette;
            // init palette
            for (int i = 0; i < 256; i++)
            {
                cp.Entries[i] = Color.FromArgb(i, i, i);
            }
            // set palette back
            image.Palette = cp;


            // lock destination bitmap data
            BitmapData imageData = BitmapConverter.Lock8bpp(image);
            float scale = 255.0f / maxMapIntensity;

            // do the job
            unsafe
            {
                byte* dst = (byte*)imageData.Scan0.ToPointer();
                int y, x;

                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, dst++)
                    {
                        *dst = Maths.Byte(scale * houghMap[y, x]);
                    }
                }
            }

            // unlock destination images
            image.UnlockBits(imageData);

            return image;
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            // lock source image
            BitmapData imageData = BitmapConverter.Lock8bpp(Data);

            try
            {
                // process the image
                Apply(imageData);
            }
            finally
            {
                // unlock image
                Data.UnlockBits(imageData);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format8bppIndexed)
            {
                throw new Exception("Image format should be: 8bppIndexed");
            }

            // get source image size
            int width = bmData.Width;
            int height = bmData.Height;
            Rectangle rect = new Rectangle(0, 0, width, height);
            int halfWidth = width / 2;
            int halfHeight = height / 2;

            // make sure the specified rectangle recides with the source image
            rect.Intersect(new Rectangle(0, 0, width, height));

            int startX = -halfWidth + rect.Left;
            int startY = -halfHeight + rect.Top;
            int stopX = width - halfWidth - (width - rect.Right);
            int stopY = height - halfHeight - (height - rect.Bottom);

            int offset = bmData.Stride - rect.Width;

            // calculate Hough map's width
            int halfHoughWidth = (int)Math.Sqrt(halfWidth * halfWidth + halfHeight * halfHeight);
            int houghWidth = halfHoughWidth * 2;

            houghMap = new short[houghHeight, houghWidth];

            // do the job
            unsafe
            {
                byte* src = (byte*)bmData.Scan0.ToPointer() + rect.Top * bmData.Stride + rect.Left;
                int y, x, theta, radius;

                // for each row
                for (y = startY; y < stopY; y++)
                {
                    // for each pixel
                    for (x = startX; x < stopX; x++, src++)
                    {
                        if (*src != 0)
                        {
                            // for each Theta value
                            for (theta = 0; theta < houghHeight; theta++)
                            {
                                radius = (int)Math.Round(cosMap[theta] * x - sinMap[theta] * y) + halfHoughWidth;

                                if ((radius < 0) || (radius >= houghWidth))
                                    continue;

                                houghMap[theta, radius]++;
                            }
                        }
                    }
                    src += offset;
                }
            }

            // find max value in Hough map
            maxMapIntensity = 0;
            int i, j;

            for (i = 0; i < houghHeight; i++)
            {
                for (j = 0; j < houghWidth; j++)
                {
                    if (houghMap[i, j] > maxMapIntensity)
                    {
                        maxMapIntensity = houghMap[i, j];
                    }
                }
            }

            CollectLines();
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Collect lines with intesities greater or equal then specified.
        /// </summary>
        private void CollectLines()
        {
            int maxTheta = houghMap.GetLength(0);
            int maxRadius = houghMap.GetLength(1);

            short intensity;
            bool foundGreater;

            int halfHoughWidth = maxRadius >> 1;
            int theta, radius, cycledTheta, cycledRadius, tt, tr, ttMax, trMax;

            // clean lines collection
            lines.Clear();

            // for each Theta value
            for (theta = 0; theta < maxTheta; theta++)
            {
                // for each Radius value
                for (radius = 0; radius < maxRadius; radius++)
                {
                    // get current value
                    intensity = houghMap[theta, radius];

                    if (intensity < minLineIntensity)
                        continue;

                    foundGreater = false;

                    // check neighboors
                    for (tt = theta - localPeakRadius, ttMax = theta + localPeakRadius; tt < ttMax; tt++)
                    {
                        // break if it is not local maximum
                        if (foundGreater == true)
                            break;

                        cycledTheta = tt;
                        cycledRadius = radius;

                        // check limits
                        if (cycledTheta < 0)
                        {
                            cycledTheta = maxTheta + cycledTheta;
                            cycledRadius = maxRadius - cycledRadius;
                        }
                        if (cycledTheta >= maxTheta)
                        {
                            cycledTheta -= maxTheta;
                            cycledRadius = maxRadius - cycledRadius;
                        }

                        for (tr = cycledRadius - localPeakRadius, trMax = cycledRadius + localPeakRadius; tr < trMax; tr++)
                        {
                            // skip out of map values
                            if (tr < 0)
                                continue;
                            if (tr >= maxRadius)
                                break;

                            // compare the neighboor with current value
                            if (houghMap[cycledTheta, tr] > intensity)
                            {
                                foundGreater = true;
                                break;
                            }
                        }
                    }

                    // was it local maximum ?
                    if (!foundGreater)
                    {
                        // we have local maximum
                        lines.Add(new HoughLine((double)theta / stepsPerDegree, (short)(radius - halfHoughWidth), intensity, (double)intensity / maxMapIntensity));
                    }
                }
            }

            lines.Sort();
        }
        #endregion
    }
    /// <summary>
    /// Defines the Hough circle transform filter.
    /// </summary>
    public class HoughCircleTransform : IBitmapFilter
    {
        #region Private data
        private int radiusToDetect;

        private short[,] houghMap;
        private short maxMapIntensity = 0;

        private int width;
        private int height;

        private int localPeakRadius = 4;
        private short minCircleIntensity = 10;
        private ArrayList circles = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes the Hough circle transform filter.
        /// </summary>
        /// <param name="radiusToDetect">Radius</param>
        public HoughCircleTransform(int radiusToDetect)
        {
            this.radiusToDetect = radiusToDetect;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Gets the minimum circle intensity.
        /// </summary>
        public short MinCircleIntensity
        {
            get { return minCircleIntensity; }
            set { minCircleIntensity = value; }
        }
        /// <summary>
        /// Gets or sets the radius search for the local peak value.
        /// </summary>
        public int LocalPeakRadius
        {
            get { return localPeakRadius; }
            set { localPeakRadius = Math.Max(1, Math.Min(10, value)); }
        }
        /// <summary>
        /// Gets the maximum line intensity.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
        }
        /// <summary>
        /// Gets the count of lines found.
        /// </summary>
        public int CirclesCount
        {
            get { return circles.Count; }
        }
        /// <summary>
        /// Returns an array of lines with absolute intensity.
        /// </summary>
        /// <param name="count">Count</param>
        /// <returns>Array</returns>
        public HoughCircle[] GetMostIntensiveCircles(int count)
        {
            // lines count
            int n = Math.Min(count, circles.Count);

            // result array
            HoughCircle[] dst = new HoughCircle[n];
            circles.CopyTo(0, dst, 0, n);

            return dst;
        }
        /// <summary>
        /// Returns an array of lines with relative intensity.
        /// </summary>
        /// <param name="minRelativeIntensity">Minimum relative intensity</param>
        /// <returns>Array</returns>
        public HoughCircle[] GetCirclesByRelativeIntensity(double minRelativeIntensity)
        {
            int count = 0, n = circles.Count;

            while ((count < n) && (((HoughCircle)circles[count]).RelativeIntensity >= minRelativeIntensity))
                count++;

            return GetMostIntensiveCircles(count);
        }
        #endregion

        #region Image processing
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        public void Apply(Bitmap Data)
        {
            // lock source image
            BitmapData imageData = BitmapConverter.Lock8bpp(Data);

            try
            {
                // process the image
                Apply(imageData);
            }
            finally
            {
                // unlock image
                Data.UnlockBits(imageData);
            }
        }
        /// <summary>
        /// Apply filter.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format8bppIndexed)
            {
                throw new Exception("Image format should be: 8bppIndexed");
            }

            // get source image size
            width = bmData.Width;
            height = bmData.Height;

            int srcOffset = bmData.Stride - width;

            // allocate Hough map of the same size like image
            houghMap = new short[height, width];

            // do the job
            unsafe
            {
                byte* src = (byte*)bmData.Scan0.ToPointer();

                // for each row
                for (int y = 0; y < height; y++)
                {
                    // for each pixel
                    for (int x = 0; x < width; x++, src++)
                    {
                        if (*src != 0)
                        {
                            DrawHoughCircle(x, y);
                        }
                    }
                    src += srcOffset;
                }
            }

            // find max value in Hough map
            maxMapIntensity = 0;
            for (int i = 0; i < height; i++)
            {
                for (int j = 0; j < width; j++)
                {
                    if (houghMap[i, j] > maxMapIntensity)
                    {
                        maxMapIntensity = houghMap[i, j];
                    }
                }
            }

            CollectCircles();
        }
        /// <summary>
        /// Returns the Hough synogram.
        /// </summary>
        /// <returns>Bitmap</returns>
        public Bitmap ToBitmap()
        {
            // check if Hough transformation was made already
            if (houghMap == null)
            {
                throw new ApplicationException("Hough conversion was not made");
            }

            int width = houghMap.GetLength(1);
            int height = houghMap.GetLength(0);

            // create new image
            Bitmap image = new Bitmap(width, height, PixelFormat.Format8bppIndexed);

            // get palette
            ColorPalette cp = image.Palette;
            // init palette
            for (int i = 0; i < 256; i++)
            {
                cp.Entries[i] = Color.FromArgb(i, i, i);
            }
            // set palette back
            image.Palette = cp;

            // lock destination bitmap data
            BitmapData imageData = BitmapConverter.Lock8bpp(image);
            float scale = 255.0f / maxMapIntensity;

            // do the job
            unsafe
            {
                byte* dst = (byte*)imageData.Scan0.ToPointer();
                int y, x;

                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, dst++)
                    {
                        *dst = Maths.Byte(scale * houghMap[y, x]);
                    }
                }
            }

            // unlock destination images
            image.UnlockBits(imageData);

            return image;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Collect circles with intesities greater or equal then specified.
        /// </summary>
        private void CollectCircles()
        {
            short intensity;
            bool foundGreater;

            // clean circles collection
            circles.Clear();
            int y, x, ty, tx, txMax, tyMax;

            // for each Y coordinate
            for (y = 0; y < height; y++)
            {
                // for each X coordinate
                for (x = 0; x < width; x++)
                {
                    // get current value
                    intensity = houghMap[y, x];

                    if (intensity < minCircleIntensity)
                        continue;

                    foundGreater = false;

                    // check neighboors
                    for (ty = y - localPeakRadius, tyMax = y + localPeakRadius; ty < tyMax; ty++)
                    {
                        // continue if the coordinate is out of map
                        if (ty < 0)
                            continue;
                        // break if it is not local maximum or coordinate is out of map
                        if ((foundGreater == true) || (ty >= height))
                            break;

                        for (tx = x - localPeakRadius, txMax = x + localPeakRadius; tx < txMax; tx++)
                        {
                            // continue or break if the coordinate is out of map
                            if (tx < 0)
                                continue;
                            if (tx >= width)
                                break;

                            // compare the neighboor with current value
                            if (houghMap[ty, tx] > intensity)
                            {
                                foundGreater = true;
                                break;
                            }
                        }
                    }

                    // was it local maximum ?
                    if (!foundGreater)
                    {
                        // we have local maximum
                        circles.Add(new HoughCircle(x, y, radiusToDetect, intensity, (double)intensity / maxMapIntensity));
                    }
                }
            }

            circles.Sort();
        }
        /// <summary>
        /// Draw Hough circle:
        /// http://www.cs.unc.edu/~mcmillan/comp136/Lecture7/circle.html
        /// TODO: more optimizations of circle drawing could be done.
        /// </summary>
        /// <param name="xCenter">Co. X</param>
        /// <param name="yCenter">Co. Y</param>
        private void DrawHoughCircle(int xCenter, int yCenter)
        {
            int x = 0;
            int y = radiusToDetect;
            int p = (5 - radiusToDetect * 4) / 4;

            SetHoughCirclePoints(xCenter, yCenter, x, y);

            while (x < y)
            {
                x++;
                if (p < 0)
                {
                    p += 2 * x + 1;
                }
                else
                {
                    y--;
                    p += 2 * (x - y) + 1;
                }
                SetHoughCirclePoints(xCenter, yCenter, x, y);
            }
        }
        /// <summary>
        /// Set circle points.
        /// </summary>
        /// <param name="cx">Cx</param>
        /// <param name="cy">Cy</param>
        /// <param name="x">Co. X</param>
        /// <param name="y">Co. Y</param>
        private void SetHoughCirclePoints(int cx, int cy, int x, int y)
        {
            if (x == 0)
            {
                SetHoughPoint(cx, cy + y);
                SetHoughPoint(cx, cy - y);
                SetHoughPoint(cx + y, cy);
                SetHoughPoint(cx - y, cy);
            }
            else if (x == y)
            {
                SetHoughPoint(cx + x, cy + y);
                SetHoughPoint(cx - x, cy + y);
                SetHoughPoint(cx + x, cy - y);
                SetHoughPoint(cx - x, cy - y);
            }
            else if (x < y)
            {
                SetHoughPoint(cx + x, cy + y);
                SetHoughPoint(cx - x, cy + y);
                SetHoughPoint(cx + x, cy - y);
                SetHoughPoint(cx - x, cy - y);
                SetHoughPoint(cx + y, cy + x);
                SetHoughPoint(cx - y, cy + x);
                SetHoughPoint(cx + y, cy - x);
                SetHoughPoint(cx - y, cy - x);
            }
        }
        /// <summary>
        /// Set point.
        /// </summary>
        /// <param name="x">Co. X</param>
        /// <param name="y">Co. Y</param>
        private void SetHoughPoint(int x, int y)
        {
            if ((x >= 0) && (y >= 0) && (x < width) && (y < height))
            {
                houghMap[y, x]++;
            }
        }
        #endregion
    }
    #endregion

    #region Mathematics
    /// <summary>
    /// Uses to work with brightness represented as a value belonging to the interval [0, 1].
    /// </summary>
    public static class Intensity
    {
        #region Nonlinear methods components
        /// <summary>
        /// Implements the Single Scale Retinex algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double SingleScaleRetinex(double x, double xlow, double nbase, double a, double b)
        {
            // Singe scale retinex modified algorithm
            // by Asiryan Valeriy
            // 
            return Math.Exp(a * Math.Log(x / xlow, nbase) + b) - 0.5f;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="nbase">Logarithm base</param>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] SingleScaleRetinex(double nbase, double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.SingleScaleRetinex(w, v, nbase, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the local contrast inversion algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double LocalContrastInversion(double x, double xlow, double a, double b)
        {
            return a * x / (xlow + b);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor (0, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LocalContrastInversion(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LocalContrastInversion(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the local contrast enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="a">Factor [-1, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double LocalContrastEnhancement(double x, double xlow, double a)
        {
            return x + a * (x - xlow);

            //return (1.0 + a) * (x - xlow) + xlow;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LocalContrastEnhancement(double a, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LocalContrastEnhancement(w, v, a);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the homomorphic enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double HomomorphicEnhancement(double x, double mu, double a, double b)
        {
            return Math.Exp(Math.Log(x) - a * Math.Log(mu + b));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset (0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] HomomorphicEnhancement(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.HomomorphicEnhancement(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the ξ-contrast enhancement algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double KsiContrastEnchancement(double x, double mu, double a, double b)
        {
            // x ∈ [0, 1], μ ∈ [0, 1] - mean of x.
            // σ - variance, ξ - coefficient.
            //
            // σ = x - μ
            // ξ = σ / μ
            // result value:
            // x' = x + α * ξ + β, where α, β ∈ [-1, 1].

            double sigma = x - mu;
            double ksi = sigma / mu;
            return x + a * ksi + b;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Contrast [-1, 1]</param>
        /// <param name="b">Offset [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] KsiContrastEnchancement(double a, double b, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.KsiContrastEnchancement(w, v, a, b);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the SAUCE algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="d">Degree of difference [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double SAUCE(double x, double mu, double d)
        {
            // Ravimal Bandara algorithm
            // implementation:

            double a = (mu - d / 2.0);
            double b = (mu + d / 2.0);

            if (x < a)
            {
                return 0;
            }
            else if (x >= b)
            {
                return 1;
            }
            return (x - a) / d;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Factor [-1, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] SAUCE(double a, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.SAUCE(w, v, a);
                }
            }

            return table;
        }
        /// <summary>
        /// Implements the Bradley threshold correction.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="xlow">Filter brightness</param>
        /// <param name="difference">Difference [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bradley(double x, double xlow, double difference = 0.15)
        {
            // Bradley local threshold void.
            // Derek Bradley, Gerhard Roth (2005). Adaptive Thresholding Using the Integral Image.
            // Retrieved from http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf

            double z = 1.0 - difference;
            return (x < xlow * z) ? 0 : 1;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="difference">Difference [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] Bradley(double difference, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.Bradley(w, v, difference);
                }
            }

            return table;
        }
        #endregion

        #region LogStretrch Histogram method components
        /// <summary>
        /// Logarithm of 0.5.
        /// </summary>
        public const double log05 = -0.693147180559945;
        /// <summary>
        /// Logarithmic epsilon.
        /// </summary>
        public const double logEpsilon = 1e-9;
        /// <summary>
        /// Implements the logarithmic stretch algorithm.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="mu">Filter brightness</param>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogStretch(double x, double mu, double s, double l)
        {
            return Intensity.LogPow(x, Maths.Range(Intensity.log05 / Math.Log(mu), s, l));
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="s">Shadows</param>
        /// <param name="l">Highlights</param>
        /// <param name="length">Length</param>
        /// <returns>Matrix</returns>
        public static double[,] LogStretch(double s, double l, int length)
        {
            double[,] table = new double[length, length];
            double w, v;
            int x, y;

            for (x = 0; x < length; x++)
            {
                w = x / (double)length;

                for (y = 0; y < length; y++)
                {
                    v = y / (double)length;

                    table[x, y] = Intensity.LogStretch(w, v, s, l);
                }
            }

            return table;
        }
        /// <summary>
        /// Returns the number raised to the logarithmic power.
        /// </summary>
        /// <param name="a">Number</param>
        /// <param name="power">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogPow(double a, double power)
        {
            return Math.Exp(Math.Log(a) * power);
        }
        #endregion

        #region Linear methods components
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="g">Gamma</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Gamma(double g, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Gamma(x / (double)length, g);
            }
            return table;
        }
        /// <summary>
        /// Implements the gamma correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="g">Gamma</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gamma(double x, double g)
        {
            return Math.Pow(x, g);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Shift(double b, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Shift(x / (double)length, b);
            }
            return table;
        }
        /// <summary>
        /// Implements the shift correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="b">Offset (-0.5, 0.5)</param>
        /// <returns>Double precision floating point number</returns>
        public static double Shift(double x, double b)
        {
            double v = log05 / Math.Log(0.5 - b);
            return LogPow(x, v);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Bin(double threshold, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Bin(x / (double)length, threshold);
            }
            return table;
        }
        /// <summary>
        /// Implements the threshold correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="threshold">Threshold [0, 1]</param>
        /// <returns>Double precision floating point number</returns>
        public static double Bin(double x, double threshold)
        {
            return (x > threshold) ? 1 : 0;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="average">Average</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Exposure(double average, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Exposure(x / (double)length, average);
            }
            return table;
        }
        /// <summary>
        /// Implements the exposure correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="average">Average</param>
        /// <returns>Double precision floating point number</returns>
        public static double Exposure(double x, double average)
        {
            double T = 255.0 / average;
            return 1 - Math.Exp(-T * x);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Sin(double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Sin(x / (double)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the sine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Sin(double x, double delta)
        {
            return 0.5 * Math.Sin((3.14 * x) - (3.14 / 2)) + 0.5 + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Cos(double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Cos(x / (double)length, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the cosine correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Cos(double x, double delta)
        {
            return 0.5 * Math.Cos((3.14 * x) - 3.14) + 0.5 + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Log(double a, double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Log(x / (double)length, a, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the logarithmic correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="a">Logarithm base</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Log(double x, double a, double delta)
        {
            return Math.Log((1.0 + (x + delta) / 0.5), a);
        }
        /// <summary>
        /// Returns the correction mask for formula: Y = (X + V).
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Add(double value, int length)
        {
            double[] table = new double[length];
            for (int x = 0; x < length; x++)
            {
                table[x] = x / (double)length + value;
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="value">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Contrast(double value, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Contrast(x / (double)length, value);
            }
            return table;
        }
        /// <summary>
        /// Implements the contrast correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="value">Contrast</param>
        /// <returns>Double precision floating point number</returns>
        public static double Contrast(double x, double value)
        {
            value = (1 + value);
            double xc = x;
            xc -= 0.5f;
            xc *= value;
            xc += 0.5f;
            return xc;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="power">Value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] LogContrast(double power, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = LogContrast(x / (double)length, power);
            }
            return table;
        }
        /// <summary>
        /// Implements the log-contrast correction.
        /// </summary>
        /// <param name="x">Brightness</param>
        /// <param name="power">Power</param>
        /// <returns>Double precision floating point number</returns>
        public static double LogContrast(double x, double power)
        {
            if (x <= 0.5)
            {
                return Intensity.LogPow(x * 2, power) * 0.5;
            }
            return 1.0 - Intensity.LogPow((1 - x) * 2, power) * 0.5;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Invert(int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Invert(x / (double)length);
            }
            return table;
        }
        /// <summary>
        /// Negates the value.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Invert(double x)
        {
            return 1.0 - x;
        }
        /// <summary>
        /// Equalizes a value relative to the {min, max} range.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <returns>Double precision floating point number</returns>
        public static double Equalize(double x, double min, double max)
        {
            double a = max - min;
            double b = x - min;
            double c = (a != 0) ? b / a : x;
            return Maths.Double(c);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="max">Maximum value</param>
        /// <param name="min">Minimum value</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Equalize(double min, double max, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Equalize(x / (double)length, min, max);
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="range">Pair of numbers Max и Min</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Linear(RangeDouble range, double delta, int length)
        {
            return Intensity.Linear(range.Max, range.Min, delta, length);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="xmax">Maximum value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="delta">Delta</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Linear(double xmax, double xmin, double delta, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Linear(x / (double)length, xmax, xmin, delta);
            }
            return table;
        }
        /// <summary>
        /// Implements the linear correction.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="xmax">Maximum value</param>
        /// <param name="xmin">Minimum value</param>
        /// <param name="delta">Delta</param>
        /// <returns>Double precision floating point number</returns>
        public static double Linear(double x, double xmax, double xmin, double delta)
        {
            return (x - xmin) / (xmax - xmin) + delta;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="xmin">Minimum value of the input range</param>
        /// <param name="xmax">Maximum value of the input range</param>
        /// <param name="ymin">Minimum value of the output range</param>
        /// <param name="ymax">Maximum value of the output range</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Levels(double xmin, double xmax, double ymin, double ymax, int length)
        {
            double[] table = new double[length];
            double k = 0, b = 0;
            double v;

            if (xmax != xmin)
            {
                k = (ymax - ymin) / (xmax - xmin);
                b = (ymin) - k * xmin;
            }

            for (int i = 0; i < length; i++)
            {
                v = i / (double)length;

                if (v >= xmax)
                { v = ymax; }
                else if (v <= xmin)
                { v = ymin; }
                else
                { v = k * v + b; }

                table[i] = v;
            }
            return table;
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="input">Input values</param>
        /// <param name="output">Output values</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Levels(RangeDouble input, RangeDouble output, int length)
        {
            return Intensity.Levels(input.Min, input.Max, output.Min, output.Max, length);
        }
        /// <summary>
        /// Returns the correction mask.
        /// </summary>
        /// <param name="levels">Number of levels</param>
        /// <param name="length">Length</param>
        /// <returns>Array</returns>
        public static double[] Quantize(int levels, int length)
        {
            if (levels > length)
                throw new Exception("Number of levels cannot be greater than length");

            int interval = length / levels + 1;
            double[] table = new double[length];
            double min = double.MaxValue;
            double max = double.MinValue;
            double v;
            int q, i;

            // calculating
            for (i = 0; i < length; i++)
            {
                q = (i / interval) * interval;
                v = q / (length - 1.0);

                if (v < min)
                    min = v;

                if (v > max)
                    max = v;

                table[i] = v;
            }

            // normalizing
            v = (max - min);

            for (i = 0; i < length; i++)
            {
                table[i] = (table[i] - min) / v;
            }

            return table;
        }
        #endregion
    }
    /// <summary>
    /// Used to blending layers.
    /// <remarks>
    /// More information can be found on the website:
    /// http://www.pegtop.net/delphi/articles/blendmodes/index.htm
    /// </remarks>
    /// </summary>
    public static class BlendMode
    {
        #region Blend mode components
        /// <summary>
        /// Implements the averaging function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Average(double a, double b)
        {
            return (a + b) / 2.0;
        }
        /// <summary>
        /// Implements the screening function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Screen(double a, double b)
        {
            return 1 - (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the difference function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Difference(double a, double b)
        {
            return Math.Abs(a - b);
        }
        /// <summary>
        /// Implements the negation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Negation(double a, double b)
        {
            return 1 - Math.Abs(1 - a - b);
        }
        /// <summary>
        /// Implements the exclusion function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Exclusion(double a, double b)
        {
            return a + b - 2 * a * b;
        }
        /// <summary>
        /// Implements the overlaying function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Overlay(double a, double b)
        {
            if (a < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the "hard light" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double HardLight(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Implements the "dodge" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Dodge(double a, double b)
        {
            return a / (1 - b);
        }
        /// <summary>
        /// Implements the "soft dodge" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double SoftDodge(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * a / (1 - b);
            }
            return 1 - 0.5 * (1 - b) / a;
        }
        /// <summary>
        /// Implements the "burn" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Burn(double a, double b)
        {
            return 1 - (1 - a) / b;
        }
        /// <summary>
        /// Implements the "soft burn" function".
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double SoftBurn(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * b / (1 - a);
            }
            return 1 - 0.5 * (1 - a) / b;
        }
        /// <summary>
        /// Implements the reflection function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Reflect(double a, double b)
        {
            return a * a / (1 - b);
        }
        /// <summary>
        /// Implements the glow function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Glow(double a, double b)
        {
            return b * b / (1 - a);
        }
        /// <summary>
        /// Implements the stamp function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Stamp(double a, double b)
        {
            return a + 2 * b - 1;
        }
        /// <summary>
        /// Implements the "freeze" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Freeze(double a, double b)
        {
            double x = 1 - a;
            return 1 - x * x / b;
        }
        /// <summary>
        /// Implements the "heat" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Heat(double a, double b)
        {
            double x = 1 - b;
            return 1 - x * x / a;
        }
        /// <summary>
        /// Implements the interpolation function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Interpolation(double a, double b)
        {
            return 0.5 - 0.25 * Math.Cos(Math.PI * a) - 0.25 * Math.Cos(Math.PI * b);
        }
        /// <summary>
        /// Implements the function of "soft light" (Adobe Photoshop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Photoshop(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b + a * a * (1 - 2 * a);
            }
            return 2 * a * (1 - b) + Math.Sqrt(a) * (2 * b - 1);
        }
        /// <summary>
        /// Implements the function of "soft light" (Illusions.hu).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Illusions(double a, double b)
        {
            double x = 2 * (0.5 - b);
            double y = Math.Pow(2, x);
            return Math.Pow(a, y);
        }
        /// <summary>
        /// Implements the function of "soft light" (Pegtop).
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Pegtop(double a, double b)
        {
            return (1 - 2 * b) * a * a + 2 * b * a;
        }
        /// <summary>
        /// Implements the "Cairo" function.
        /// </summary>
        /// <param name="a">First layer</param>
        /// <param name="b">Second layer</param>
        /// <returns>Double precision floating point number</returns>
        public static double Fw3c(double a, double b)
        {
            if (b <= 0.5)
            {
                return a - (1 - 2 * b) * a * (1 - a);
            }
            return a + (2 * b - 1) * (BlendMode.Gw3c(a) - a);
        }
        /// <summary>
        /// Implements the "Cairo" function.
        /// </summary>
        /// <param name="a">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public static double Gw3c(double a)
        {
            if (a <= 0.25)
            {
                return ((16 * a - 12) * a + 4) * a;
            }
            return Math.Sqrt(a);
        }
        #endregion
    }
    /// <summary>
    /// Uses to work with the statistical characteristics of the image.
    /// </summary>
    public static class Statistics
    {
        #region Histogram
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Array</returns>
        public unsafe static int[] Histogram(BitmapData bmData)
        {
            int[] rgb = new int[256];
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;
            byte brightness;

            for (y = 0; y < height; y++)
            {
                for (x = 0; x < width; x++, p += 4)
                {
                    brightness = (byte)RGB.Average(p[2], p[1], p[0]);
                    rgb[brightness]++;
                }
            }
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="channel">Channel of RGBA model</param>
        /// <returns>Array</returns>
        public unsafe static int[] Histogram(BitmapData bmData, RGBA channel)
        {
            int[] rgb = new int[256];
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int y, x, width = bmData.Width, height = bmData.Height;

            if (channel == RGBA.Red)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[2])]++;
                    }
                }
            }
            else if (channel == RGBA.Green)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[1])]++;
                    }
                }
            }
            else if (channel == RGBA.Blue)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[0])]++;
                    }
                }
            }
            else if (channel == RGBA.Alpha)
            {
                for (y = 0; y < height; y++)
                {
                    for (x = 0; x < width; x++, p += 4)
                    {
                        rgb[(int)(p[3])]++;
                    }
                }
            }
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Array</returns>
        public static int[] Histogram(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            int[] rgb = Histogram(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Gets a histogram of the image.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <param name="channel">Channel of RGBA model</param>
        /// <returns>Array</returns>
        public static int[] Histogram(Bitmap Data, RGBA channel)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            int[] rgb = Histogram(bmData, channel);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        #endregion

        #region Density
        /// <summary>
        /// Gets an array of values of the density function.
        /// </summary>
        /// <param name="H">Histogram</param>
        /// <returns>Array</returns>
        public static int[] CDF(int[] H)
        {
            int length = H.Length;
            int[] cdf = new int[length];
            cdf[0] = H[0];

            for (int i = 1; i < length; i++)
            {
                cdf[i] = H[i] + cdf[i - 1];
            }
            return cdf;
        }
        /// <summary>
        /// Gets an array of equalized histogram values by recalculating the brightness density function.
        /// </summary>
        /// <param name="H">Histogram</param>
        /// <returns>Array</returns>
        public static double[] Equalize(int[] H)
        {
            // CDF calculating
            int length = H.Length;
            int[] cdf = Statistics.CDF(H);
            int min = Statistics.MinIndex(cdf);
            int max = Statistics.MaxIndex(cdf);

            // table
            double[] table = new double[length];
            double factor = cdf[max] - cdf[min];

            // scaling
            for (int i = 0; i < length; i++)
            {
                table[i] = (cdf[i] - cdf[min]) / factor;
            }

            return table;
        }
        #endregion

        #region Otsu's threshold
        /// <summary>
        /// Calculates the optimal threshold using the Otsu method for the original bitmap.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Integer number</returns>
        public static int OtsuThreshold(BitmapData bmData)
        {
            double[] v = new double[256];
            int[] h = Statistics.Histogram(bmData);
            double diff, w1, w2;

            for (int k = 1; k < 256; k++)
            {
                w1 = Omega(0, k, h);
                w2 = Omega(k + 1, 255, h);
                diff = Mu(0, k, h) * w2 - Mu(k + 1, 255, h) * w1;
                v[k] = (diff * diff) / (w1 * w2);
            }

            return MaxIndex(v);
        }
        /// <summary>
        /// Calculates the optimal threshold using the Otsu method for the original bitmap.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Integer number</returns>
        public static int OtsuThreshold(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            int threshold = OtsuThreshold(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return threshold;
        }
        #endregion

        #region SIS threshold
        /// <summary>
        /// Calculates the optimal threshold for the original bitmap.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Integer number</returns>
        public unsafe static int SISThreshold(BitmapData bmData)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int width1 = width - 1, height1 = height - 1, offset = stride - width;
            double ex, ey, weight, weightTotal = 0, total = 0;
            byte* p = (byte*)bmData.Scan0.ToPointer(); p += stride;
            int y, x;

            for (y = 1; y < height1; y++)
            {
                p++;

                for (x = 1; x < width1; x++, p++)
                {
                    ex = Math.Abs(p[1] - p[-1]);
                    ey = Math.Abs(p[stride] - p[-stride]);
                    weight = (ex > ey) ? ex : ey;
                    weightTotal += weight;
                    total += weight * (*p);
                }
                p += offset + 1;
            }

            return (weightTotal == 0) ? 0 : (byte)(total / weightTotal);
        }
        /// <summary>
        /// Calculates the optimal threshold for the original bitmap.
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Integer number</returns>
        public static int SISThreshold(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            int threshold = SISThreshold(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return threshold;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Gets the index of the maximum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        private static int MaxIndex(double[] data)
        {
            int index = 0, length = data.Length - 1, i;
            double maximum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] > maximum)
                {
                    maximum = data[i];
                    index = i;
                }
            }

            return index;
        }
        /// <summary>
        /// Gets the index of the maximum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        private static int MaxIndex(int[] data)
        {
            int index = 0, length = data.Length - 1, i;
            int maximum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] > maximum)
                {
                    maximum = data[i];
                    index = i;
                }
            }

            return index;
        }
        /// <summary>
        /// Gets the index of the minimum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        private static int MinIndex(double[] data)
        {
            int index = 0, length = data.Length - 1, i;
            double minimum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] < minimum)
                {
                    minimum = data[i];
                    index = i;
                }
            }

            return index;
        }
        /// <summary>
        /// Gets the index of the minimum element of the array.
        /// </summary>
        /// <param name="data">Array</param>
        /// <returns>Integer number</returns>
        private static int MinIndex(int[] data)
        {
            int index = 0, length = data.Length - 1, i;
            int minimum = data[index];

            for (i = 1; i < length; i++)
            {
                if (data[i] < minimum)
                {
                    minimum = data[i];
                    index = i;
                }
            }

            return index;
        }
        /// <summary>
        /// Omega.
        /// </summary>
        /// <param name="init">Init</param>
        /// <param name="end">End</param>
        /// <param name="h">Histogram</param>
        /// <returns>Double precision floating point number</returns>
        private static double Omega(int init, int end, int[] h)
        {
            int sum = 0, i;

            for (i = init; i <= end; i++)
            {
                sum += h[i];
            }

            return sum;
        }
        /// <summary>
        /// Mean.
        /// </summary>
        /// <param name="init">Init</param>
        /// <param name="end">End</param>
        /// <param name="h">Histogram</param>
        /// <returns>Double precision floating point number</returns>
        private static double Mu(int init, int end, int[] h)
        {
            int sum = 0, i;

            for (i = init; i <= end; i++)
            {
                sum += i * h[i];
            }

            return sum;
        }
        #endregion
    }
    /// <summary>
    /// Uses to work with point matrices.
    /// </summary>
    public static class PointMatrix
    {
        #region Filters
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] FlipY(int width, int height)
        {
            PointInt[,] matrix = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y;

                for (y = 0; y < height; y++)
                {

                    matrix[x, y].X = x;
                    matrix[x, y].Y = height - y - 1;
                }
            }
            );

            return matrix;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] FlipX(int width, int height)
        {
            PointInt[,] matrix = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y;

                for (y = 0; y < height; y++)
                {

                    matrix[x, y].X = width - x - 1;
                    matrix[x, y].Y = y;
                }
            }
            );

            return matrix;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Offset</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] ShiftX(int width, int height, int value)
        {
            PointInt[,] matrix = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y;

                for (y = 0; y < height; y++)
                {

                    matrix[x, y].X = Maths.Mod(x + value, width);
                    matrix[x, y].Y = y;
                }
            }
            );

            return matrix;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Offset</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] ShiftY(int width, int height, int value)
        {
            PointInt[,] matrix = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y;

                for (y = 0; y < height; y++)
                {

                    matrix[x, y].X = x;
                    matrix[x, y].Y = Maths.Mod(y + value, height);
                }
            }
            );

            return matrix;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Value [0, 100]</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] Noise(int width, int height, int value)
        {
            PointInt[,] noise = new PointInt[width, height];
            int y, x, newX, newY;
            int nHalf = value / 2;
            Random rnd = new Random();

            for (x = 0; x < width; x++)
            {
                for (y = 0; y < height; y++)
                {
                    newX = rnd.Next(value) - nHalf;

                    if (x + newX > 0 && x + newX < width)
                        noise[x, y].X = newX;
                    else
                        noise[x, y].X = 0;

                    newY = rnd.Next(value) - nHalf;

                    if (y + newY > 0 && y + newY < width)
                        noise[x, y].Y = newY;
                    else
                        noise[x, y].Y = 0;
                }
            }
            return noise;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Value [0, 100]</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] Pixelate(int width, int height, int value)
        {
            PointInt[,] pixelate = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y, newX, newY;

                for (y = 0; y < height; y++)
                {
                    newX = (int)(value - x % value);
                    newY = (int)(value - y % value);

                    pixelate[x, y].X = newX - value;
                    pixelate[x, y].Y = newY - value;
                }
            }
            );

            return pixelate;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Value [0, 100]</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] Grid(int width, int height, int value)
        {
            PointInt[,] grid = new PointInt[width, height];

            Parallel.For(0, width, x =>
            {
                int y, newX, newY;
                
                for (y = 0; y < height; y++)
                {
                    newX = (int)(value - x % value);
                    newY = (int)(value - y % value);

                    if (newX == value)
                        grid[x, y].X = -x;
                    else
                        grid[x, y].X = newX - value;

                    if (newY == value)
                        grid[x, y].Y = -y;
                    else
                        grid[x, y].Y = newY - value;
                }
            }
            );

            return grid;
        }
        /// <summary>
        /// Returns the point matrix.
        /// </summary>
        /// <param name="width">Image width</param>
        /// <param name="height">Image height</param>
        /// <param name="value">Value [0, 100]</param>
        /// <returns>Array of ordered pairs of X and Y</returns>
        public static PointInt[,] Water(int width, int height, int value)
        {
            PointInt[,] water = new PointInt[width, height];
            double pix = 2.0 * Math.PI / 127.5;

            Parallel.For(0, width, x =>
            {
                int y;
                double x0, y0;

                y0 = value * Math.Cos(pix * x);

                for (y = 0; y < height; y++)
                {
                    x0 = value * Math.Sin(pix * y);

                    water[x, y].X = (int)Maths.Mod(x + x0, width);
                    water[x, y].Y = (int)Maths.Mod(y + y0, height);
                }
            }
            );

            return water;
        }
        #endregion
    }
    /// <summary>
    /// Defines the Perlin noise.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Perlin_noise
    /// </remarks>
    /// </summary>
    public class PerlinNoise
    {
        #region Private data
        private double initFrequency = 1.0;
        private double initAmplitude = 1.0;
        private double persistence = 0.65;
        private int octaves = 4;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes the Perlin noise.   
        /// </summary>
        /// <param name="octaves">Octaves[1, 32]</param>
        /// <param name="persistence">Persistence</param>
        /// <param name="frequency">Frequency</param>
        /// <param name="amplitude">Amplitude</param>
        public PerlinNoise(int octaves = 4, double persistence = 0.65, double frequency = 1, double amplitude = 1)
        {
            Octaves = octaves; Persistence = persistence; Frequency = frequency; Amplitude = amplitude;
        }
        /// <summary>
        /// Gets or sets the frequency.
        /// </summary>
        public double Frequency
        {
            get { return initFrequency; }
            set { initFrequency = value; }
        }
        /// <summary>
        /// Gets or sets the amplitude value.
        /// </summary>
        public double Amplitude
        {
            get { return initAmplitude; }
            set { initAmplitude = value; }
        }
        /// <summary>
        /// Gets or sets the persistence value.
        /// </summary>
        public double Persistence
        {
            get { return persistence; }
            set { persistence = value; }
        }
        /// <summary>
        /// Gets or sets the number of octaves[1, 32].
        /// </summary>
        public int Octaves
        {
            get { return octaves; }
            set { octaves = System.Math.Max(1, System.Math.Min(32, value)); }
        }
        /// <summary>
        /// One-dimensional Perlin noise function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            double frequency = initFrequency;
            double amplitude = initAmplitude;
            double sum = 0;

            // octaves
            for (int i = 0; i < octaves; i++)
            {
                sum += PerlinNoise.SmoothedNoise(x * frequency) * amplitude;

                frequency *= 2;
                amplitude *= persistence;
            }
            return sum;
        }
        /// <summary>
        /// Two-dimensional Perlin noise function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="y">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function2D(double x, double y)
        {
            double frequency = initFrequency;
            double amplitude = initAmplitude;
            double sum = 0;

            // octaves
            for (int i = 0; i < octaves; i++)
            {
                sum += PerlinNoise.SmoothedNoise(x * frequency, y * frequency) * amplitude;

                frequency *= 2;
                amplitude *= persistence;
            }
            return sum;
        }
        #endregion

        #region Private static voids
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double Noise(int x)
        {
            int n = (x << 13) ^ x;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private static double Noise(int x, int y)
        {
            int n = x + y * 57;
            n = (n << 13) ^ n;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        private static double SmoothedNoise(double x)
        {
            int xInt = (int)x;
            double xFrac = x - xInt;

            return PerlinNoise.CosineInterpolate(Noise(xInt), Noise(xInt + 1), xFrac);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        private static double SmoothedNoise(double x, double y)
        {
            // params
            int xInt = (int)x;
            int yInt = (int)y;
            double xFrac = x - xInt;
            double yFrac = y - yInt;

            // get four noise values
            double x0y0 = PerlinNoise.Noise(xInt, yInt);
            double x1y0 = PerlinNoise.Noise(xInt + 1, yInt);
            double x0y1 = PerlinNoise.Noise(xInt, yInt + 1);
            double x1y1 = PerlinNoise.Noise(xInt + 1, yInt + 1);

            // x interpolation
            double v1 = PerlinNoise.CosineInterpolate(x0y0, x1y0, xFrac);
            double v2 = PerlinNoise.CosineInterpolate(x0y1, x1y1, xFrac);

            // y interpolation
            return PerlinNoise.CosineInterpolate(v1, v2, yFrac);
        }
        ///
        private static double CosineInterpolate(double x1, double x2, double a)
        {
            double f = (1 - Math.Cos(a * Math.PI)) * 0.5;

            return x1 * (1 - f) + x2 * f;
        }
        #endregion
    }
    #endregion

    #region Enums
    /// <summary>
    /// Defines the color space.
    /// </summary>
    public enum Space
    {
        #region Space
        /// <summary>
        /// Color space RGB.
        /// </summary>
        RGB,
        /// <summary>
        /// Color space HSB.
        /// </summary>
        HSB,
        /// <summary>
        /// Color space HSB.
        /// </summary>
        HSL,
        /// <summary>
        /// Color space YCbCr.
        /// </summary>
        YCbCr,
        /// <summary>
        /// Grayscale.
        /// </summary>
        Grayscale,
        #endregion
    }
    /// <summary>
    /// Defines the channel of RGBA model.
    /// </summary>
    public enum RGBA
    {
        #region RGBA
        /// <summary>
        /// Alpha.
        /// </summary>
        Alpha = 3,
        /// <summary>
        /// Red.
        /// </summary>
        Red = 2,
        /// <summary>
        /// Green.
        /// </summary>
        Green = 1,
        /// <summary>
        /// Blue.
        /// </summary>
        Blue = 0,
        #endregion
    }
    /// <summary>
    /// Defines the direction of the gradient vector.
    /// </summary>
    public enum Gradient
    {
        #region Gradient
        /// <summary>
        /// North direction.
        /// </summary>
        North = 0,
        /// <summary>
        /// North-West direction.
        /// </summary>
        NorthWest = 1,
        /// <summary>
        /// West direction.
        /// </summary>
        West = 2,
        /// <summary>
        /// South-West direction.
        /// </summary>
        SouthWest = 3,
        /// <summary>
        /// South direction.
        /// </summary>
        South = 4,
        /// <summary>
        /// South-East direction.
        /// </summary>
        SouthEast = 5,
        /// <summary>
        /// East direction.
        /// </summary>
        East = 6,
        /// <summary>
        /// North-East direction.
        /// </summary>
        NorthEast = 7,
        #endregion
    }
    #endregion
}
