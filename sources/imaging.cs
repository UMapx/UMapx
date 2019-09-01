// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Collections;
using System.Collections.Generic;
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
    /// Определяет общий интерфейс фильтров двух изображений.
    /// </summary>
    public interface IBitmapFilter2
    {
        #region Filter components
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        void Apply(BitmapData bmData, BitmapData bmSrc);
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        void Apply(Bitmap Data, Bitmap Src);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс фильтров.
    /// </summary>
    public interface IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        void Apply(BitmapData bmData);
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        void Apply(Bitmap Data);
        #endregion
    }
    /// <summary>
    /// Определяет общий вид холстов.
    /// </summary>
    public interface ICanvas
    {
        #region Interface
        /// <summary>
        /// Получает ширину холста.
        /// </summary>
        int Width { get; set; }
        /// <summary>
        /// Получает высоту холста.
        /// </summary>
        int Height { get; set; }
        /// <summary>
        /// Создает холст.
        /// </summary>
        /// <returns>Точечный рисунок</returns>
        Bitmap Create();
        #endregion
    }
    #endregion

    #region Abstract classes
    /// <summary>
    /// Определяет абстрактный класс перестроения данных.
    /// </summary>
    public abstract class Rebuilder
    {
        #region Protected components
        /// <summary>
        /// Использовать перестроение данных класса или нет.
        /// </summary>
        protected bool rebuild = false;
        /// <summary>
        /// Реализует перестроение данных класса.
        /// </summary>
        protected abstract void Rebuild();
        #endregion
    }
    #endregion

    #region Rebuilder filters
    /// <summary>
    /// Определяет фильтр масочной коррекции.
    /// </summary>
    public class Correction : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Значения.
        /// </summary>
        protected double[] values;
        /// <summary>
        /// Цветовое пространство.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр масочной коррекции.
        /// </summary>
        /// <param name="values">Одномерная маска</param>
        /// <param name="space">Цветовое пространство</param>
        public Correction(double[] values, Space space)
        {
            Values = values;
            Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр масочной коррекции.
        /// </summary>
        public Correction()
        {
            Values = new double[256];
            Space = Imaging.Space.RGB;
        }
        /// <summary>
        /// Получает или задает табулированную маску.
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
                    throw new Exception("Размер маски должен быть равен 256");

                this.values = value;
            }
        }
        /// <summary>
        /// Получает или задает цветовое пространство.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
    /// Определяет фильтр масочной коррекции каналов.
    /// </summary>
    public class RGBACorrection : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Канал.
        /// </summary>
        protected RGBA channel;
        /// <summary>
        /// Значения.
        /// </summary>
        protected double[] values;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр масочной коррекции каналов.
        /// </summary>
        /// <param name="values">Одномерная маска</param>
        /// <param name="channel">Канал</param>
        public RGBACorrection(double[] values, RGBA channel)
        {
            Values = values;
            Channel = channel;
        }
        /// <summary>
        /// Инициализирует фильтр масочной коррекции каналов.
        /// </summary>
        public RGBACorrection()
        {
            Values = new double[256];
        }
        /// <summary>
        /// Получает или задает табулированную маску.
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
                    throw new Exception("Размер маски должен быть равен 256");

                this.values = value;
            }
        }
        /// <summary>
        /// Получает или задает цветовой канал модели RGBA.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр локальной масочной коррекции.
    /// </summary>
    public class LocalCorrection : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Фильтр локального усреднения.
        /// </summary>
        protected BoxBlur gb;
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double[,] values;
        /// <summary>
        /// Цветовое пространство.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции.
        /// </summary>
        public LocalCorrection()
        {
            gb = new BoxBlur(3);
            Space = Space.RGB;
            Values = new double[256, 256];
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="values">Двумерная маска</param>
        public LocalCorrection(int radius, double[,] values, Space space)
        {
            gb = new BoxBlur(radius);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="values">Двумерная маска</param>
        public LocalCorrection(int width, int height, double[,] values, Space space)
        {
            gb = new BoxBlur(width, height);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="values">Двумерная маска</param>
        public LocalCorrection(SizeInt size, double[,] values, Space space)
        {
            gb = new BoxBlur(size);
            Space = space;
            Values = values;
        }
        /// <summary>
        /// Получает или задает размер фильтра.
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
        /// Получает или задает цветовое пространство.
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
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
                    throw new Exception("Размер маски должен быть равен 256");

                this.values = value;
            }
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
    /// Определяет фильтр локальной масочной коррекции каналов.
    /// </summary>
    public class RGBALocalCorrection : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Фильтр локального усреднения.
        /// </summary>
        private BoxBlur gb;
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double[,] values;
        /// <summary>
        /// Канал.
        /// </summary>
        protected RGBA channel;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции каналов.
        /// </summary>
        public RGBALocalCorrection()
        {
            gb = new BoxBlur(3);
            Channel = RGBA.Red;
            Values = new double[256, 256];
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции каналов.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="channel">Канал</param>
        /// <param name="values">Двумерная маска</param>
        public RGBALocalCorrection(int radius, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(radius);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции каналов.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="channel">Канал</param>
        /// <param name="values">Двумерная маска</param>
        public RGBALocalCorrection(int width, int height, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(width, height);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Инициализирует фильтр локальной масочной коррекции каналов.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        /// <param name="channel">Канал</param>
        /// <param name="values">Двумерная маска</param>
        public RGBALocalCorrection(SizeInt size, double[,] values, RGBA channel)
        {
            gb = new BoxBlur(size);
            Channel = channel; 
            Values = values;
        }
        /// <summary>
        /// Получает или задает размер фильтра.
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
        /// Получает или задает цветовой канал модели RGBA.
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
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
                    throw new Exception("Размер маски должен быть равен 256");

                this.values = value;
            }
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр суммирования точек.
    /// </summary>
    public class PointAddition : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Массив точек.
        /// </summary>
        protected PointInt[,] points = new PointInt[0, 0];
        /// <summary>
        /// Ширина изображения.
        /// </summary>
        protected int width;
        /// <summary>
        /// Высота изображения.
        /// </summary>
        protected int height;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр суммирования точек.
        /// </summary>
        /// <param name="points">Массив упорядоченных пар чисел X и Y</param>
        public PointAddition(PointInt[,] points)
        {
            Points = points;
        }
        /// <summary>
        /// Инициализирует фильтр суммирования точек.
        /// </summary>
        public PointAddition() { }
        /// <summary>
        /// Получает или задает массив упорядоченных пар чисел X и Y.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр умножения точек.
    /// </summary>
    public class PointMultiplication : Rebuilder, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Массив точек.
        /// </summary>
        protected PointInt[,] points = new PointInt[0, 0];
        /// <summary>
        /// Ширина изображения.
        /// </summary>
        protected int width;
        /// <summary>
        /// Высота изображения.
        /// </summary>
        protected int height;
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр умножения точек.
        /// </summary>
        /// <param name="points">Массив упорядоченных пар чисел X и Y</param>
        public PointMultiplication(PointInt[,] points)
        {
            Points = points;
        }
        /// <summary>
        /// Инициализирует фильтр умножения точек.
        /// </summary>
        public PointMultiplication() { }
        /// <summary>
        /// Получает или задает массив упорядоченных пар чисел X и Y.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр инверсии яркости.
    /// </summary>
    public class InvertChannels : Correction, IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Инициализирует фильтр инверсии яркости.
        /// </summary>
        public InvertChannels(Space space)
        {
            this.values = Intensity.Invert(256);
        }
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild() { }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр бинариазции.
    /// </summary>
    public class Threshold : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Порог.
        /// </summary>
        protected double threshold;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует новый фильтр.
        /// </summary>
        /// <param name="threshold">Пороговое значение [0, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public Threshold(double threshold, Space space)
        {
            this.Value = threshold;
            this.Space = space;
        }
        /// <summary>
        /// Инициализирует новый фильтр.
        /// </summary>
        public Threshold()
        {
            Value = 0.5;
        }
        /// <summary>
        /// Получает или задает значение порогового значения [0, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bin(this.threshold, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции уровней.
    /// <remarks>
    /// Пример использования фильтра:
    /// https://digital-photography-school.com/using-levels-photoshop-image-correct-color-contrast/
    /// </remarks>
    /// </summary>
    public class LevelsCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Входные значения каналов.
        /// </summary>
        protected RangeDouble input;
        /// <summary>
        /// Выходные значения каналов.
        /// </summary>
        protected RangeDouble output;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции уровней.
        /// </summary>
        /// <param name="input">Входные значения каналов</param>
        /// <param name="output">Выходные значения каналов</param>
        /// <param name="space">Цветовое пространство</param>
        public LevelsCorrection(RangeDouble input, RangeDouble output, Space space)
        {
            Input = input; Output = output;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции уровней.
        /// </summary>
        public LevelsCorrection()
        {
            Input = new RangeDouble(0, 1);
            Output = new RangeDouble(0, 1);
        }
        /// <summary>
        /// Получает или задает входные значения каналов.
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
        /// Получает или задает выходные значения каналов.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Levels(input, output, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр экспозиционной коррекции.
    /// <remarks>
    /// Более подробное описание фильтра можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Exposure_(photography)
    /// </remarks>
    /// </summary>
    public class ExposureCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Среднее число.
        /// </summary>
        protected double average;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр экспозиционной коррекции.
        /// </summary>
        /// <param name="average">Среднее число [0, 2500]</param>
        /// <param name="space">Цветовое пространство</param>
        public ExposureCorrection(double average, Space space)
        {
            Average = average;
        }
        /// <summary>
        /// Инициализирует фильтр экспозиционной коррекции.
        /// </summary>
        public ExposureCorrection()
        {
            Average = 128;
        }
        /// <summary>
        /// Получает или задает значение среднего числа [0, 2500].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Exposure(this.average, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр линейной коррекции.
    /// </summary>
    public class LinearCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Граничные значения каналов.
        /// </summary>
        protected RangeDouble range;
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр линейной коррекции.
        /// </summary>
        /// <param name="range">Граничные значения каналов</param>
        /// <param name="delta">Дельта [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public LinearCorrection(RangeDouble range, double delta, Space space)
        {
            Range = range; Delta = delta;
        }
        /// <summary>
        /// Инициализирует фильтр линейной коррекции.
        /// </summary>
        /// <param name="delta">Дельта [-100, 100]</param>
        /// <param name="space">Цветовое пространство</param>
        public LinearCorrection(double delta, Space space)
        {
            Range = new RangeDouble(0, 1); Delta = delta;
        }
        /// <summary>
        /// Инициализирует фильтр линейной коррекции.
        /// </summary>
        public LinearCorrection()
        {
            Range = new RangeDouble(0, 1); Delta = 0.5;
        }
        /// <summary>
        /// Получает или задает граничные значения каналов.
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
        /// Получает или задает значение дельты [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Linear(range, delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр синусоидальной коррекции.
    /// </summary>
    public class SinCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Дельта.
        /// </summary>
        protected double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр синусоидальной коррекции.
        /// </summary>
        /// <param name="delta">Дельта [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public SinCorrection(double delta, Space space)
        {
            Delta = delta;
        }
        /// <summary>
        /// Инициализирует фильтр синусоидальной коррекции.
        /// </summary>
        public SinCorrection()
        {
            Delta = 20;
        }
        /// <summary>
        /// Получает или задает значение дельты [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Sin(delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр косинусоидальной коррекции.
    /// </summary>
    public class CosCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Дельта.
        /// </summary>
        protected double delta;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр косинусоидальной коррекции.
        /// </summary>
        /// <param name="delta">Дельта [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public CosCorrection(double delta, Space space)
        {
            Delta = delta;
        }
        /// <summary>
        /// Инициализирует фильтр синусоидальной коррекции.
        /// </summary>
        public CosCorrection()
        {
            Delta = 0.5;
        }
        /// <summary>
        /// Получает или задает значение дельты.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Cos(delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр логарифмической коррекции.
    /// </summary>
    public class LogCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Основание логарифма.
        /// </summary>
        protected double nbase = 3.14f;
        /// <summary>
        /// Дельта.
        /// </summary>
        protected double delta = 0.5;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр логарифмической коррекции.
        /// </summary>
        /// <param name="nbase">Основание логарифма</param>
        /// <param name="delta">Дельта [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public LogCorrection(double nbase, double delta, Space space)
        {
            Base = nbase; Delta = delta;
        }
        /// <summary>
        /// Инициализирует фильтр логарифмической коррекции.
        /// </summary>
        public LogCorrection() { }
        /// <summary>
        /// Получает или задает значение основания логарифма.
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
        /// Получает или задает значение дельты.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Log(nbase, delta / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции яркости.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    /// </summary>
    public class BrightnessCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Яркость.
        /// </summary>
        protected double brightness;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции якрости.
        /// </summary>
        /// <param name="brightness">Яркость [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public BrightnessCorrection(double brightness, Space space)
        {
            Brightness = brightness;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции якрости.
        /// </summary>
        public BrightnessCorrection()
        {
            Brightness = 0.5;
        }
        /// <summary>
        /// Получает или задает значение якрости [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.brightness / 2.0, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции контрастности.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// http://esate.ru/uroki/OpenGL/image_processing/_p4106/
    /// </remarks>
    /// </summary>
    public class ContrastCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Контрастность.
        /// </summary>
        protected double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции контрастности.
        /// </summary>
        /// <param name="value">Контрастность [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public ContrastCorrection(double value, Space space)
        {
            Contrast = value;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции контрастности.
        /// </summary>
        public ContrastCorrection()
        {
            Contrast = 0.5;
        }
        /// <summary>
        /// Получает или задает значение контрастности.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Contrast(this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр гамма-коррекции.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gamma_correction
    /// </remarks>
    /// </summary>
    public class GammaCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Гамма.
        /// </summary>
        protected double g;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр гамма-коррекции.
        /// </summary>
        /// <param name="g">Гамма [0, 20]</param>
        /// <param name="space">Цветовое пространство</param>
        public GammaCorrection(double g, Space space)
        {
            Gamma = g; Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр гамма-коррекции.
        /// </summary>
        public GammaCorrection()
        {
            Gamma = 2.2;
        }
        /// <summary>
        /// Получает или задает значение гаммы [0, 20].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Gamma(this.g, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции ко смещением.
    /// </summary>
    public class ShiftCorrection : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double offset;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции ко смещением.
        /// </summary>
        /// <param name="offset">Смещение (-0.5, 0.5)</param>
        /// <param name="space">Цветовое пространство</param>
        public ShiftCorrection(double offset, Space space)
        {
            Offset = offset; Space = space;
        }
        /// <summary>
        /// Получает или задает значение смещения (-0.5, 0.5).
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Shift(this.offset, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр глобального улучшения контраста.
    /// </summary>
    public class ContrastEnhancement : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр глобального улучшения контраста.
        /// </summary>
        /// <param name="contrast">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public ContrastEnhancement(double contrast, Space space)
        {
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LogContrast(1 + this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр квантования.
    /// <remarks>
    /// Примеры работы фильтра можно найти на сайте: 
    /// http://en.wikipedia.org/wiki/Posterization
    /// </remarks>
    /// </summary>
    public class Quantization : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Количество уровней.
        /// </summary>
        protected int levels;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр квантования.
        /// </summary>
        /// <param name="levels">Количество уровней</param>
        /// <param name="space">Цветовое пространство</param>
        public Quantization(int levels, Space space)
        {
            Levels = levels; Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр квантования.
        /// </summary>
        public Quantization()
        {
            Levels = 4;
        }
        /// <summary>
        /// Получает или задает количество уровней.
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
        /// Реализует перестроение данных фильтра.
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
    /// Определяет фильтр локальной бинаризации Брэдли.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf
    /// </remarks>
    /// </summary>
    public class LocalThreshold : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Предел разницы яркости между пикселем обработки и средним значением локальных пикселей.
        /// </summary>
        protected double difference;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр локальной бинаризации Брэдли.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="difference">Предел разницы яркости между пикселем обработки и средним значением локальных пикселей [0, 1]</param>
        public LocalThreshold(int radius, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(radius);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр локальной бинаризации Брэдли.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="difference">Предел разницы яркости между пикселем обработки и средним значением локальных пикселей [0, 1]</param>
        public LocalThreshold(int width, int height, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(width, height);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр локальной бинаризации Брэдли.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="difference">Предел разницы яркости между пикселем обработки и средним значением локальных пикселей [0, 1]</param>
        public LocalThreshold(SizeInt size, Space space, double difference = 0.15)
        {
            gb = new BoxBlur(size);
            Difference = difference;
            Space = space;
        }
        /// <summary>
        /// Получает или задает предел разницы яркости между пикселем обработки и средним значением локальных пикселей.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Bradley(this.difference, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр локального улучшения контраста.
    /// <remarks>
    /// Данный фильтр также известен под названием "Unsharp Masking". 
    /// Подробное описание алгоритма можно найти на сайте:
    /// http://www.cambridgeincolour.com/tutorials/local-contrast-enhancement.htm
    /// Примеры использования:
    /// http://www.knowhowtransfer.com/photoshop-professional-plugins/alce-local-contrast-enhancer/
    /// </remarks>
    /// </summary>
    public class LocalContrastEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double contrast;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр локального улучшения контраста.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [-1, 1]</param>
        public LocalContrastEnhancement(int radius, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(radius);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Инициализирует фильтр локального улучшения контраста.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [-1, 1]</param>
        public LocalContrastEnhancement(int width, int height, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(width, height);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Инициализирует фильтр локального улучшения контраста.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [-1, 1]</param>
        public LocalContrastEnhancement(SizeInt size, Space space, double contrast = 0.75)
        {
            gb = new BoxBlur(size);
            Space = space; Contrast = contrast;
        }
        /// <summary>
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastEnhancement(this.contrast, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр инверсии локального контраста.
    /// <remarks>
    /// Данный фильтр используется для выравнивания освещенности изображений путем усреднения яркости.
    /// </remarks>
    /// </summary>
    public class LocalContrastInversion : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double a;
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр инверсии локального контраста.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public LocalContrastInversion(int radius, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(radius);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Инициализирует фильтр инверсии локального контраста.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public LocalContrastInversion(int width, int height, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(width, height);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Инициализирует фильтр инверсии локального контраста.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public LocalContrastInversion(SizeInt size, Space space, double a = 0.75, double b = 0.05)
        {
            this.gb = new BoxBlur(size);
            Space = space;
            A = a;
            B = b;
        }
        /// <summary>
        /// Получает или задает значение коэффициента сжатия контраста (0, 1].
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
        /// Получает или задает значение смещения (0, 1].
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
                    throw new Exception("Неверное значение аргумента");

                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.LocalContrastInversion(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр повышения контраста.
    /// </summary>
    public class KsiContrastEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double a;
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр повышения контраста.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение [-1, 1]</param>
        public KsiContrastEnhancement(int radius, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Инициализирует фильтр повышения контраста.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение [-1, 1]</param>
        public KsiContrastEnhancement(int width, int height, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Инициализирует фильтр повышения контраста.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение [-1, 1]</param>
        public KsiContrastEnhancement(SizeInt size, Space space, double a = 0.75, double b = 0.05)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
        /// Получает или задает значение смещения [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.KsiContrastEnchancement(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр гомоморфной обработки.
    /// <remarks>
    /// Гомоморфный фильтр чаще всего используется для выравнивания освещенности изображений. 
    /// Он одновременно нормализует яркость изображения и увеличивает контрастность.
    /// 
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Homomorphic_filtering
    /// </remarks>
    /// </summary>
    public class HomomorphicEnhancement : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double a;
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double b;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр гомоморфной обработки.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public HomomorphicEnhancement(int radius, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Инициализирует фильтр гомоморфной обработки.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public HomomorphicEnhancement(int width, int height, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Инициализирует фильтр гомоморфной обработки.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        public HomomorphicEnhancement(SizeInt size, Space space, double a = 0.5, double b = 0.05)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b;
        }
        /// <summary>
        /// Получает или задает значение коэффициента сжатия контраста [-1, 1].
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
        /// Получает или задает значение смещения (0, 1].
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
                    throw new Exception("Неверное значение аргумента");

                this.b = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.HomomorphicEnhancement(this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр Single Scale Retinex.
    /// <remarks>
    /// Алгоритмы группы Single Scale Retinex предназначены для коррекции неравномерно освещенных изображений. 
    /// Фильтр данной группы выравнивает освещенность изображения, сохраняя локальный контраст в плохо и ярко освещенных областях.
    /// 
    /// Более подробную информацию можно найти в статье:
    /// https://dragon.larc.nasa.gov/background/pubabs/papers/gspx1.pdf
    /// </remarks>
    /// </summary>
    public class SingleScaleRetinex : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double a;
        /// <summary>
        /// Смещение.
        /// </summary>
        protected double b;
        /// <summary>
        /// Основание логарифма.
        /// </summary>
        protected double nbase;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр Single Scale Retinex.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="nbase">Основание логарифма</param>
        public SingleScaleRetinex(int radius, Space space, double a = 1, double b = 0, double nbase = Math.PI)
        {
            gb = new BoxBlur(radius);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Инициализирует фильтр Single Scale Retinex.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="nbase">Основание логарифма</param>
        public SingleScaleRetinex(int width, int height, Space space, double a = 1, double b = 0, double nbase = Math.PI)
        {
            gb = new BoxBlur(width, height);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Инициализирует фильтр Single Scale Retinex.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="nbase">Основание логарифма</param>
        public SingleScaleRetinex(SizeInt size, Space space, double a = 1, double b = 0, double nbase = Math.PI)
        {
            gb = new BoxBlur(size);
            Space = space; A = a; B = b; Base = nbase;
        }
        /// <summary>
        /// Получает или задает коэффициент сжатия контраста [-1, 1].
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
        /// Получает или задает значение смещения (0, 1].
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
        /// Получает или задает основание логарифма.
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
                    throw new Exception("Основание логарифма должно быть больше 0");

                this.nbase = value;
                this.rebuild = true;
            }
        }
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.SingleScaleRetinex(this.nbase, this.a, this.b, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции теней и светов.
    /// <remarks>
    /// Коррекция теней и светов (Shadows-Highlights correction) используется для коррекции неравномерно освещенных изображений. В отличие от других локальных алгоритмов
    /// (например, Single Scale Retinex, Homomorphic Enhancement, Flat-Field Correction) фильтр позволяет регулировать значения яркости отдельно в темных и светлых областях
    /// изображения.
    /// </remarks>
    /// </summary>
    public class ShadowsHighlightsCorrection : LocalCorrection, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Тени.
        /// </summary>
        protected double shadows;
        /// <summary>
        /// Света.
        /// </summary>
        protected double lights;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции теней и светов.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="shadows">Тени [0, 1]</param>
        /// <param name="highlights">Света [0, 1]</param>
        public ShadowsHighlightsCorrection(int radius, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(radius);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции теней и светов.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="shadows">Тени [0, 1]</param>
        /// <param name="highlights">Света [0, 1]</param>
        public ShadowsHighlightsCorrection(int width, int height, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(width, height);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции теней и светов.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="shadows">Тени [0, 1]</param>
        /// <param name="highlights">Света [0, 1]</param>
        public ShadowsHighlightsCorrection(SizeInt size, Space space, double shadows = 0.4, double highlights = 0.4)
        {
            gb = new BoxBlur(size);
            Space = space; Shadows = shadows; Highlights = highlights;
        }
        /// <summary>
        /// Получает или задает значение теней [0, 1].
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
        /// Получает или задает значение светов [0, 1].
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
        /// Реализует перестроение данных фильтра.
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
        /// Возвращает значение гамма-параметра.
        /// </summary>
        /// <param name="v">Значение фильтра [0, 1]</param>
        /// <returns>Гамма</returns>
        private double value2gamma(double v)
        {
            // Вычисление гамма-параметра:
            return (v - Intensity.logEpsilon) / 2.0;
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр Flat-Field-коррекции.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
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
        /// Инициализирует фильтр "Flat-Field" коррекции.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public FlatFieldCorrection(int radius = 15)
        {
            gb = new BoxBlur(radius);
        }
        /// <summary>
        /// Инициализирует фильтр "Flat-Field" коррекции.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public FlatFieldCorrection(int width, int height)
        {
            gb = new BoxBlur(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр "Flat-Field" коррекции.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        public FlatFieldCorrection(SizeInt size)
        {
            gb = new BoxBlur(size);
        }
        /// <summary>
        /// Получает или задает размер фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            gb.Apply(bmSrc);
            flatfield(bmData, bmSrc);
            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Получает значения медиан изображения.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Одномерный массив</returns>
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
    /// Определяет фильтр на основе структуры HSB.
    /// <remarks>
    /// Фильтр обесцвечивает указанную часть изображения.
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
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="hue">Диапазон оттенков [0, 359]</param>
        public HSBGrayscale(RangeInt hue)
        {
            Hue = hue;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="min">Нижний предел [0, 359]</param>
        /// <param name="max">Верхний предел [0, 359]</param>
        public HSBGrayscale(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Получает или задает диапазон оттенков [0, 359].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр на основе структуры HSL.
    /// <remarks>
    /// Фильтр обесцвечивает указанную часть изображения.
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
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="hue">Диапазон оттенков [0, 359]</param>
        public HSLGrayscale(RangeInt hue)
        {
            Hue = hue;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="min">Нижний предел [0, 359]</param>
        /// <param name="max">Верхний предел [0, 359]</param>
        public HSLGrayscale(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Получает или задает диапазон оттенков [0, 359].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр градаций серого.
    /// </summary>
    public class Grayscale : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// C(r).
        /// </summary>
        protected double cr;
        /// <summary>
        /// C(g).
        /// </summary>
        protected double cg;
        /// <summary>
        /// C(b).
        /// </summary>
        protected double cb;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр градаций серого.
        /// </summary>
        /// <param name="cr">Коэффициент красного канала</param>
        /// <param name="cg">Коэффициент зеленого канала</param>
        /// <param name="cb">Коэффициент синего канала</param>
        public Grayscale(double cr, double cg, double cb)
        {
            Cr = cr;
            Cg = cg;
            Cb = cb;
        }
        /// <summary>
        /// Инициализирует фильтр градаций серого.
        /// </summary>
        public Grayscale()
        {
            Cr = 0.333f;
            Cg = 0.333f;
            Cb = 0.333f;
        }
        /// <summary>
        /// Получает или задает значение коэффициента красного канала.
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
        /// Получает или задает значение коэффициента зеленого канала.
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
        /// Получает или задает значение коэффициента синего канала.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static voids
        /// <summary>
        /// Инициализирует фильтр градаций серого (BT709).
        /// </summary>
        public static Grayscale BT709
        {
            get
            {
                return new Grayscale(0.212f, 0.715f, 0.072f);
            }
        }
        /// <summary>
        /// Определяет фильтр градаций серого (R-Y).
        /// </summary>
        public static Grayscale RY
        {
            get
            {
                return new Grayscale(0.5f, 0.419f, 0.081f);
            }
        }
        /// <summary>
        /// Определяет фильтр градаций серого Y.
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
        /// Проверяет является ли точечный рисунок изображением в оттенках серого.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Логическое значение</returns>
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
    /// Определяет цветовой фотофильтр.
    /// </summary>
    public class PhotoFilter : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Сила фильтра.
        /// </summary>
        protected double s;
        /// <summary>
        /// Цвет.
        /// </summary>
        protected Color color;
        /// <summary>
        /// Функция смешивания.
        /// </summary>
        protected IDoubleMesh blendf;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует цветовой фотофильтр.
        /// </summary>
        /// <param name="blendf">Функция смешивания</param>
        /// <param name="color">Цвет фильтра</param>
        /// <param name="strength">Сила фильтра [0, 1]</param>
        public PhotoFilter(IDoubleMesh blendf, Color color, double strength = 0.5)
        {
            BlendFunction = blendf;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Инициализирует цветовой фотофильтр.
        /// </summary>
        /// <param name="color">Цвет фильтра</param>
        /// <param name="strength">Сила фильтра [0, 1]</param>
        public PhotoFilter(Color color, double strength = 0.5)
        {
            BlendFunction = BlendMode.Pegtop;
            Color = color;
            Strength = strength;
        }
        /// <summary>
        /// Инициализирует цветовой фотофильтр.
        /// </summary>
        public PhotoFilter()
        {
            BlendFunction = BlendMode.Pegtop;
            Color = Color.White;
            Strength = 0.5;
        }
        /// <summary>
        /// Получает или задает цвет фильтра.
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
        /// Получает или задает цвет фильтра.
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
        /// Получает или задает силу фильтра [0, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Инициализирует холодный фильтр (82).
        /// </summary>
        public static PhotoFilter Cold82
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 181, 255));
            }
        }
        /// <summary>
        /// Инициализирует холодный фильтр LBB.
        /// </summary>
        public static PhotoFilter ColdLBB
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 93, 255));
            }
        }
        /// <summary>
        /// Инициализирует теплый фильтр (81).
        /// </summary>
        public static PhotoFilter Warm81
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(235, 177, 19));
            }
        }
        /// <summary>
        /// Инициализирует теплый фильтр LBA.
        /// </summary>
        public static PhotoFilter WarmLBA
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(250, 150, 0));
            }
        }
        /// <summary>
        /// Инициализирует фильтр сепии.
        /// </summary>
        public static PhotoFilter Sepia
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Инициализирует красный фильтр.
        /// </summary>
        public static PhotoFilter Red
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(234, 26, 26));
            }
        }
        /// <summary>
        /// Инициализирует синий фильтр.
        /// </summary>
        public static PhotoFilter Blue
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(29, 53, 234));
            }
        }
        /// <summary>
        /// Инициализирует зеленый фильтр.
        /// </summary>
        public static PhotoFilter Green
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(25, 201, 25));
            }
        }
        /// <summary>
        /// Инициализирует подводный фильтр.
        /// </summary>
        public static PhotoFilter Underwater
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(0, 194, 177));
            }
        }
        /// <summary>
        /// Инициализирует пурпурный фильтр.
        /// </summary>
        public static PhotoFilter Purple
        {
            get
            {
                return new PhotoFilter(Color.FromArgb(227, 24, 227));
            }
        }
        /// <summary>
        /// Инициализирует оранжевый фильтр.
        /// </summary>
        public static PhotoFilter Orange
        {
            get
            {
                return new PhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Инициализирует желтый фильтр.
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
    /// Определяет фильтр коррекции температуры.
    /// <remarks>
    /// Фильтр использует аппроксимацию кривой Планка.
    /// </remarks>
    /// </summary>
    public class TemperatureCorrection : PhotoFilter, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Температура.
        /// </summary>
        protected double temperature;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции температуры.
        /// </summary>
        /// <param name="temperature">Температура [1E3К, 1E4К]</param>
        /// <param name="strength">Сила фильтра [0, 1]</param>
        public TemperatureCorrection(double temperature, double strength = 0.5)
        {
            Temperature = temperature; Strength = strength;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции температуры.
        /// </summary>
        public TemperatureCorrection()
        {
            Temperature = 1000; Strength = 0.5;
        }
        /// <summary>
        /// Получает или задает значение температуры [1E3К, 1E4К].
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
                this.color = RGB.Temp2RGB(this.temperature);
            }
        }
        #endregion
    }
    /// <summary>
    /// Определяет цветовой фотофильтр на основе модели YUV.
    /// </summary>
    public class YUVPhotoFilter : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Сила фильтра.
        /// </summary>
        protected double s;
        /// <summary>
        /// Цвет.
        /// </summary>
        protected Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует цветовой фотофильтр на основе модели YUV.
        /// </summary>
        /// <param name="color">Цвет фильтра</param>
        /// <param name="strength">Сила фильтра [0, 1]</param>
        public YUVPhotoFilter(Color color, double strength = 0.5)
        {
            Color = color; Strength = strength;
        }
        /// <summary>
        /// Инициализирует цветовой фотофильтр на основе модели YUV.
        /// </summary>
        public YUVPhotoFilter()
        {
            Color = Color.White; Strength = 0.5;
        }
        /// <summary>
        /// Получает или задает цвет фильтра.
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
        /// Получает или задает силу фильтра [0, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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

                    // Исходные значения каналов:
                    iR = p[k + 2]; iG = p[k + 1]; iB = p[k];

                    //Вычисление нового цвета:
                    luma = RGB.HDTV(iR, iG, iB);
                    nYUV = YUV.FromRGB(luma, luma, luma);
                    iYUV = YUV.FromRGB(color.R, color.G, color.B);
                    rgb  = AddColor(nYUV, iYUV).ToRGB;
                    nR   = rgb.Red; nG = rgb.Green; nB = rgb.Blue;

                    // Фиксация значений в соответствии с силой фильтра:
                    p[k + 2] = Maths.Byte(nR * s + iR * z);
                    p[k + 1] = Maths.Byte(nG * s + iG * z);
                    p[k    ] = Maths.Byte(nB * s + iB * z);
                }
            }
            );

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        /// <summary>
        /// Проверяет является ли цвет оттенком серого.
        /// </summary>
        /// <param name="color">Цвет в терминах красного, зеленого и синего каналов</param>
        /// <returns>Логическое значение</returns>
        public static bool IsGrayColor(Color color)
        {
            if (color.R == color.G && color.G == color.B)
            {
                return true;
            }
            return false;
        }
        /// <summary>
        /// Смешивает два цвета в пространстве YUV.
        /// </summary>
        /// <param name="yuv1">Первый цвет</param>
        /// <param name="yuv2">Второй цвет</param>
        /// <returns>Результат</returns>
        public static YUV AddColor(YUV yuv1, YUV yuv2)
        {
            return new YUV(yuv1.Y, yuv1.U + yuv2.U, yuv1.V + yuv2.V);
        }
        #endregion

        #region Complete filters
        /// <summary>
        /// Инициализирует фильтр сепии.
        /// </summary>
        public static YUVPhotoFilter Sepia
        {
            get
            {
                return new YUVPhotoFilter(Color.FromArgb(172, 122, 51));
            }
        }
        /// <summary>
        /// Инициализирует оранжевый фильтр.
        /// </summary>
        public static YUVPhotoFilter Orange
        {
            get
            {
                return new YUVPhotoFilter(Color.Orange);
            }
        }
        /// <summary>
        /// Инициализирует желтый фильтр.
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
    /// Определяет фильтр коррекции насыщенности.
    /// </summary>
    public class SaturationCorrection : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Насыщенность.
        /// </summary>
        protected double saturation;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции насыщенности.
        /// </summary>
        /// <param name="saturation">Насыщенность [-100, 100]</param>
        public SaturationCorrection(double saturation)
        {
            Saturation = saturation;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции насыщенности.
        /// </summary>
        public SaturationCorrection()
        {
            Saturation = 20;
        }
        /// <summary>
        /// Получает или задает значение насыщенности [-100, 100].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр замены цвета.
    /// </summary>
    public class ColorReplace : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Входной цвет.
        /// </summary>
        protected Color input = Color.Transparent;
        /// <summary>
        /// Выходной цвет.
        /// </summary>
        protected Color output = Color.Transparent;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр замены цвета.
        /// </summary>
        /// <param name="input">Входной цвет в терминах каналов красного, зеленого, синего и альфа каналов</param>
        /// <param name="output">Выходной цвет в терминах каналов красного, зеленого, синего и альфа каналов</param>
        public ColorReplace(Color input, Color output)
        {
            Input = input;
            Output = output;
        }
        /// <summary>
        /// Инициализирует фильтр замены цвета.
        /// </summary>
        public ColorReplace() { }
        /// <summary>
        /// Входной цвет в терминах каналов красного, зеленого, синего и альфа каналов.
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
        /// Выходной цвет в терминах каналов красного, зеленого, синего и альфа каналов.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр тонального дитеринга.
    /// <remarks>
    /// Примеры работы данного фильтра можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    public class ToneDiffusionDithering : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Матрица диффузии.
        /// </summary>
        protected double[,] matrix;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр тонального дитеринга.
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
        /// Инициализирует фильтр тонального дитеринга.
        /// </summary>
        /// <param name="matrix">Двумерная матрица</param>
        public ToneDiffusionDithering(double[,] matrix)
        {
            Matrix = matrix;
        }
        /// <summary>
        /// Получает или задает двумерную матрицу.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Генератор случайных чисел.
        /// </summary>
        private static System.Random rand = new System.Random();
        /// <summary>
        /// Инициализирует фильтр заказного дитеринга.
        /// <remarks>
        /// Подробное описание фильтра можно найти на сайте:
        /// http://en.wikipedia.org/wiki/Ordered_dithering
        /// Примеры работы данного фильтра можно найти на сайте:
        /// https://en.wikipedia.org/wiki/Dither
        /// </remarks>
        /// </summary>
        /// <param name="radius">Радиус матрицы [0, 255]</param>
        /// <returns>Фильтр тонального дитеринга</returns>
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
        /// Инициализирует фильтр случайного дитеринга.
        /// </summary>
        /// <param name="radius">Радиус матрицы [0, 255]</param>
        /// <returns>Фильтр тонального дитеринга</returns>
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
        /// Инициализирует фильтр классического дитеринга.
        /// </summary>
        /// <returns>Фильтр тонального дитеринга</returns>
        public static ToneDiffusionDithering Basic()
        {
            return new ToneDiffusionDithering(new double[4, 4] {
			                            {  15, 143,  47, 175 },
			                            { 207,  79, 239, 111 },
			                            {  63, 191,  31, 159 },
			                            { 255, 127, 223,  95 }});
        }
        /// <summary>
        /// Инициализирует фильтр дитеринга по алгоритму Байера.
        /// </summary>
        /// <returns>Фильтр тонального дитеринга</returns>
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
    /// Определяет фильтр дитеринга диффузии ошибки.
    /// <remarks>
    /// Примеры работы данного фильтра можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Dither
    /// </remarks>
    /// </summary>
    public class ErrorDiffusionDithering : Rebuilder, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Координата X.
        /// </summary>
        protected int x;
        /// <summary>
        /// Координата Y.
        /// </summary>
        protected int y;
        /// <summary>
        /// Ширина изображения.
        /// </summary>
        protected int width;
        /// <summary>
        /// Высота изображения.
        /// </summary>
        protected int height;
        /// <summary>
        /// Шаг изображения.
        /// </summary>
        protected int stride;
        /// <summary>
        /// Матрица.
        /// </summary>
        protected double[][] matrix;
        /// <summary>
        /// Суммарное значение.
        /// </summary>
        private double summary;
        /// <summary>
        /// Количество уровней квантования.
        /// </summary>
        protected int levels;
        /// <summary>
        /// Таблица квантования.
        /// </summary>
        private double[] table;
        #endregion

        #region Class components
        /// <summary>
        /// Инициализирует фильтр дитеринга диффузии ошибки.
        /// </summary>
        /// <param name="levels">Количество уровней квантования</param>
        /// <param name="matrix">Матрица</param>
        public ErrorDiffusionDithering(int levels, double[][] matrix)
        {
            this.Levels = levels;
            this.Matrix = matrix;
        }
        /// <summary>
        /// Получает или задает количество уровней квантования.
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
        /// Получает или задает матрицу преобразования.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return;
        }
        /// <summary>
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.table = Intensity.Quantize(this.levels, 256).Mul(255);
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Выполняет диффузию ошибки.
        /// </summary>
        /// <param name="rError">Ошибка красного канала</param>
        /// <param name="gError">Ошибка зеленого канала</param>
        /// <param name="bError">Ошибка синего канала</param>
        /// <param name="ptr">Текущий пиксель</param>
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
        /// Вычисляет ближайший цвет для заданных значений каналов.
        /// </summary>
        /// <param name="red">Красный</param>
        /// <param name="green">Зеленый</param>
        /// <param name="blue">Синий</param>
        /// <param name="table">Таблица цветов</param>
        /// <returns>Цвет в терминах красного, зеленого и синего каналов</returns>
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
        /// Инициализирует фильтр дитеринга Эткинсона.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Буркеса.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Фэна.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Сиерры.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Сиерры.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Сиерры.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Флойда-Стеинберга.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Джарвиса-Джадиса-Нинке.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Стивенсона.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Шиау.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
        /// Инициализирует фильтр дитеринга Стаки.
        /// </summary>
        /// <returns>Фильтр дитеринга диффузии ошибки</returns>
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
    /// Определяет фильтр аддитивного Гауссова шума.
    /// <remarks>
    /// Гауссовский шум представляет собой статистический шум, имеющий функцию плотности вероятности (PDF), равную функции нормального распределения, 
    /// которая также известна как распределение Гаусса. Другими словами, значения, которые может принимать шум, распределены по Гауссу.
    /// Примеры работы данного фильтра можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gaussian_noise
    /// </remarks>
    /// </summary>
    public class AdditiveNoise : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Сила.
        /// </summary>
        protected int amount = 10;
        /// <summary>
        /// Генератор случайных чисел.
        /// </summary>
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр аддитивного Гауссова шума.
        /// </summary>
        public AdditiveNoise() { }
        /// <summary>
        /// Инициализирует фильтр аддитивного Гауссова шума.
        /// </summary>
        /// <param name="amount">Процентная доля шума [0, 100]</param>
        public AdditiveNoise(int amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Получает или задает процентную долю шума [0, 100].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр шума "соль и перец".
    /// <remarks>
    /// Salt and pepper (англ. "соль с черным перцем", то есть чередование серых и белых частиц) - одна из форм шума, который как правило встречается 
    /// на графических и видео изображениях. Он представляет собой случайно возникающие белые и черные пиксели. Очень часто для проверки видео фильтров 
    /// данный шум используют как тестовый, добавляя к сигналу. В обычных же условиях шум Salt and Pepper возникает в изображения при быстрых переходных 
    /// процессах, таких как неправильная коммутация. Эффективным способом подавления этого типа шума является использование медианного фильтра.
    /// Примеры работы данного фильтра можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Salt-and-pepper_noise
    /// </remarks>
    /// </summary>
    public class SaltAndPepper : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Сила.
        /// </summary>
        protected double amount = 10;
        /// <summary>
        /// Генератор случайных чисел.
        /// </summary>
        private Random generator = new Random();
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр шума "соль и перец".
        /// </summary>
        public SaltAndPepper() { }
        /// <summary>
        /// Инициализирует фильтр шума "соль и перец".
        /// </summary>
        /// <param name="amount">Процентная доля шума [0, 100].</param>
        public SaltAndPepper(double amount)
        {
            Amount = amount;
        }
        /// <summary>
        /// Получает или задает процентную долю шума [0, 100].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр коррекции прозрачности.
    /// </summary>
    public class TransparencyCorrection : RGBACorrection, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Прозрачность.
        /// </summary>
        protected double transparency;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции прозрачности.
        /// </summary>
        /// <param name="transparency">Прозрачность [-1, 1]</param>
        public TransparencyCorrection(double transparency)
        {
            this.Channel = RGBA.Alpha;
            Transparency = transparency;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции прозрачности.
        /// </summary>
        public TransparencyCorrection()
        {
            this.Channel = RGBA.Alpha;
            Transparency = 0;
        }
        /// <summary>
        /// Получает или задает значение прозрачности [-1, 1].
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Add(this.transparency, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр коррекции уровней канала.
    /// <remarks>
    /// Пример использования фильтра:
    /// https://digital-photography-school.com/using-levels-photoshop-image-correct-color-contrast/
    /// </remarks>
    /// </summary>
    public class LevelsChannelCorrection : RGBACorrection, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Входные значения каналов.
        /// </summary>
        protected RangeDouble input;
        /// <summary>
        /// Выходные значения каналов.
        /// </summary>
        protected RangeDouble output;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр коррекции уровней.
        /// </summary>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        /// <param name="input">Входные значения каналов</param>
        /// <param name="output">Выходные значения каналов</param>
        public LevelsChannelCorrection(RGBA channel, RangeDouble input, RangeDouble output)
        {
            this.Channel = channel; Input = input; Output = output;
        }
        /// <summary>
        /// Инициализирует фильтр коррекции уровней.
        /// </summary>
        public LevelsChannelCorrection()
        {
            Input = new RangeDouble(0, 255);
            Output = new RangeDouble(0, 255);
        }
        /// <summary>
        /// Получает или задает входные значения каналов.
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
        /// Получает или задает выходные значения каналов.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Levels(input, output, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр поворота каналов.
    /// <remarks>
    /// Каналы изображения меняются местами: B = R, G = B, R = G.
    /// </remarks>
    /// </summary>
    public class RotateChannel : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Инициализирует фильтр поворота каналов.
        /// </summary>
        public RotateChannel() { }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр эквализации каналов.
    /// <remarks>
    /// Каналы изображения эквализируются в соотвествии с выбранным каналом C: R = G = B = C.
    /// </remarks>
    /// </summary>
    public class EqualizeChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр эквализации каналов.
        /// </summary>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        public EqualizeChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Инициализирует фильтр эквализации каналов.
        /// </summary>
        public EqualizeChannel() { }
        /// <summary>
        /// Получает или задает цветовой канал модели RGBA.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр скрытия канала.
    /// </summary>
    public class HideChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр скрытия канала.
        /// </summary>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        public HideChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Инициализирует фильтр скрытия канала.
        /// </summary>
        public HideChannel() { }
        /// <summary>
        /// Получает или задает цветовой канал модели RGBA.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр отображения каналов.
    /// <remarks>
    /// Каналы изображения скрываются в соотвествии с выбранным каналом.
    /// </remarks>
    /// </summary>
    public class ShowChannel : IBitmapFilter
    {
        #region Private data
        private RGBA channel = RGBA.Red;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр отображения каналов.
        /// </summary>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        public ShowChannel(RGBA channel)
        {
            Channel = channel;
        }
        /// <summary>
        /// Инициализирует фильтр отображения каналов.
        /// </summary>
        public ShowChannel() { }
        /// <summary>
        /// Получает или задает цветовой канал модели RGBA.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет класс цветного холста.
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
        /// Инициализирует класс цветного холста.
        /// </summary>
        public CanvasColor() { }
        /// <summary>
        /// Инициализирует класс цветного холста.
        /// </summary>
        /// <param name="width">Ширина холста</param>
        /// <param name="height">Высота холста</param>
        /// <param name="color">Цвет в терминах каналов красного, зеленого и синего</param>
        public CanvasColor(int width, int height, Color color)
        {
            Width = width; Height = height; Color = color;
        }
        /// <summary>
        /// Получает или задает цвет холста.
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
        /// Получает или задает ширину холста.
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
        /// Получает или задает высоту холста.
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
        /// Создает холст.
        /// </summary>
        /// <returns>Точечный рисунок</returns>
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
    /// Определяет класс градиентного холста.
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
        /// Инициализирует класс градиентного холста.
        /// </summary>
        public CanvasGradient() { }
        /// <summary>
        /// Инициализирует класс градиентного холста.
        /// </summary>
        /// <param name="width">Ширина холста</param>
        /// <param name="height">Высота холста</param>
        /// <param name="angle">Направление линейного градиента</param>
        /// <param name="color1">Первый цвет</param>
        /// <param name="color2">Второй цвет</param>
        public CanvasGradient(int width, int height, double angle, Color color1, Color color2)
        {
            Width = width; Height = height; Angle = angle; Color1 = color1; Color2 = color2;
        }
        /// <summary>
        /// Получает или задает первый цвет холста.
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
        /// Получает или задает второй цвет холста.
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
        /// Получает или задает направление линейного градиента.
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
        /// Получает или задает ширину холста.
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
        /// Получает или задает высоту холста.
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
        /// Создает холст.
        /// </summary>
        /// <returns>Точечный рисунок</returns>
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
    /// Определяет фильтр пикселизации.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Pixelate : PointAddition, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Параметр фильтра.
        /// </summary>
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр пикселизации.
        /// </summary>
        /// <param name="value">Глубина [0, 100]</param>
        public Pixelate(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Инициализирует фильтр пикселизации.
        /// </summary>
        public Pixelate() { }
        /// <summary>
        /// Получает или задает значение глубины пикселизации.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Pixelate(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр разбиения на сетку.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Grid : PointAddition, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Параметр фильтра.
        /// </summary>
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр разбиения на сетку.
        /// </summary>
        /// <param name="value">Глубина [0, 100]</param>
        public Grid(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Инициализирует фильтр разбиения на сетку.
        /// </summary>
        public Grid() { }
        /// <summary>
        /// Получает или задает значение глубины разбиения.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Grid(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр дрожания.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Jitter : PointAddition, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Параметр фильтра.
        /// </summary>
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр дрожания.
        /// </summary>
        /// <param name="value">Глубина [0, 100]</param>
        public Jitter(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Инициализирует фильтр дрожания.
        /// </summary>
        public Jitter() { }
        /// <summary>
        /// Получает или задает значение глубины дрожания.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Noise(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Определяет водный фильтр.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Water : PointMultiplication, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Параметр фильтра.
        /// </summary>
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует водный фильтр.
        /// </summary>
        /// <param name="value">Глубина [0, 100]</param>
        public Water(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Инициализирует водный фильтр.
        /// </summary>
        public Water() { }
        /// <summary>
        /// Получает или задает значение глубины фильтра.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.points = PointMatrix.Water(this.width, this.height, this.value);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр шума.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// https://www.codeproject.com/Articles/2122/Image-Processing-for-Dummies-with-C-and-GDI-Part
    /// </remarks>
    /// </summary>
    public class Noise : PointMultiplication, IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Параметр фильтра.
        /// </summary>
        private int value = 20;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует водный фильтр.
        /// </summary>
        /// <param name="value">Глубина [0, 100]</param>
        public Noise(int value)
        {
            Value = value;
        }
        /// <summary>
        /// Инициализирует водный фильтр.
        /// </summary>
        public Noise() { }
        /// <summary>
        /// Получает или задает значение глубины фильтра.
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
        /// Реализует перестроение данных фильтра.
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
    /// Определяет фильтр сдвига.
    /// </summary>
    public class Shift : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Смещение по оси X.
        /// </summary>
        protected int x;
        /// <summary>
        /// Смещение по оси Y.
        /// </summary>
        protected int y;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр сдвига.
        /// </summary>
        /// <param name="x">Значение смещения по оси X</param>
        /// <param name="y">Значение смещения по оси Y</param>
        public Shift(int x, int y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Инициализирует фильтр сдвига.
        /// </summary>
        /// <param name="point">Пара целых чисел, представляющих упорядоченную пару координат X и Y</param>
        public Shift(PointInt point)
        {
            X = point.X;
            Y = point.Y;
        }
        /// <summary>
        /// Инициализирует фильтр сдвига по оси X.
        /// </summary>
        public Shift() { }
        /// <summary>
        /// Получает или задает значение смещения по оси X.
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
        /// Получает или задает значение смещения по оси Y.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // image properties:
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int stride = bmSrc.Stride;

            // exception!
            if (bmData.Width  != width ||
                bmData.Height != height)
                throw new Exception("Размеры входных изображений должны быть одинаковыми");

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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="stride">Ширина шага по индексу</param>
        private unsafe void ShiftY(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1; // Смещения.
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j + y; // Смещение по оси Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i; // Смещение по оси X

                    if (sY < height && sY >= 0)
                    {
                        // Смещение в пределах рисунка
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front    ];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sY < 0) ? -1 : 1; // Определение направления смещения.
                        // Смещение за пределами рисунка
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="stride">Ширина шага по индексу</param>
        private unsafe void ShiftX(BitmapData bmData, BitmapData bmSrc, int width, int height, int stride)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            byte* pSrc = (byte*)(void*)bmSrc.Scan0.ToPointer();
            int sX = 0, sY = 0, front = 0, reverse = 0, f = 1; // Смещения.
            int i, j;

            for (j = 0; j < height; j++)
            {
                sY = j; // Смещение по оси Y

                for (i = 0; i < width; i++, p += 4)
                {
                    sX = i + x; // Смещение по оси X

                    if (sX < width && sX >= 0)
                    {
                        // Смещение в пределах рисунка
                        front = Math.Abs(sX * 4 + sY * stride);
                        p[0] = pSrc[front];
                        p[1] = pSrc[front + 1];
                        p[2] = pSrc[front + 2];
                        p[3] = pSrc[front + 3];
                    }
                    else
                    {
                        f = (sX < 0) ? -1 : 1; // Определение направления смещения.
                        // Смещение за пределами рисунка
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
    /// Определяет фильтр отображения.
    /// </summary>
    public class Flip : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// По оси X.
        /// </summary>
        protected bool x = true;
        /// <summary>
        /// По оси Y.
        /// </summary>
        protected bool y = true;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр отображения.
        /// </summary>
        public Flip() { }
        /// <summary>
        /// Инициализирует фильтр отображения.
        /// </summary>
        /// <param name="x">Отображение по оси X</param>
        /// <param name="y">Отображение по оси Y</param>
        public Flip(bool x, bool y)
        {
            X = x;
            Y = y;
        }
        /// <summary>
        /// Получает или задает отображение по оси X.
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
        /// Получает или задает отображение по оси Y.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        public unsafe void Apply(BitmapData bmData)
        {
            #region Data
            int pixel = 4, start = 0;
            int width = bmData.Width, height = bmData.Height;
            int offset = width * pixel, stride = bmData.Stride;
            #endregion

            #region FlipY
            // Отображение по оси X:
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
            // Отображение по оси Y:
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр поворота.
    /// </summary>
    public class Rotate : IBitmapFilter2
    {
        #region Private data
        private double angle;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр поворота.
        /// </summary>
        /// <param name="angle">Угол поворота</param>
        /// <param name="color">Цвет фона</param>
        public Rotate(double angle, Color color)
        {
            Angle = angle; Color = color;
        }
        /// <summary>
        /// Получает или задает угол поворота.
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
        /// Получает или задает цвет фона.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // exception!
            if (bmSrc.Width != bmData.Width || bmSrc.Height != bmData.Height)
                throw new Exception("Размеры входных изображений должны быть одинаковыми");

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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр вырезания.
    /// </summary>
    public class Crop : IBitmapFilter2
    {
        #region Private data
        Rectangle rectangle;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр вырезания.
        /// </summary>
        /// <param name="rectangle">Прямоугольная область</param>
        public Crop(Rectangle rectangle)
        {
            Rectangle = rectangle;
        }
        /// <summary>
        /// Инициализирует фильтр вырезания.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        /// <param name="width">Ширина</param>
        /// <param name="height">Высота</param>
        public Crop(int x, int y, int width, int height)
        {
            Rectangle = new Rectangle(x, y, width, height);
        }
        /// <summary>
        /// Получает или задает прямоугольную область.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
    /// Определяет фильтр изменения размеров изображения.
    /// </summary>
    public class Resize : IBitmapFilter2
    {
        #region Private data
        int newWidth;
        int newHeight;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр изменения размеров изображения.
        /// </summary>
        /// <param name="width">Ширина</param>
        /// <param name="height">Высота</param>
        public Resize(int width = 512, int height = 512)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр изменения размеров изображения.
        /// </summary>
        /// <param name="size">Ширина и высота</param>
        public Resize(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Получает или задает новый размер изображения.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
    /// Определяет фильтр наложения.
    /// </summary>
    public class Merge : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Прозрачность.
        /// </summary>
        protected int transparency = 255;
        /// <summary>
        /// Координаты X и Y.
        /// </summary>
        protected PointInt point;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр наложения.
        /// </summary>
        /// <param name="point">Пара целых чисел, представляющих упорядоченную пару координат X и Y</param>
        /// <param name="transparency">Прозрачность [0, 255]</param>
        public Merge(PointInt point, int transparency = 128)
        {
            Point = point;
            Transparency = transparency;
        }
        /// <summary>
        /// Инициализирует фильтр наложения.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        /// <param name="transparency">Прозрачность [0, 255]</param>
        public Merge(int x, int y, int transparency)
        {
            Point = new PointInt(x, y);
            Transparency = transparency;
        }
        /// <summary>
        /// Инициализирует фильтр наложения.
        /// </summary>
        public Merge() { }
        /// <summary>
        /// Получает или задает значение прозрачности [0, 255]
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
        /// Получает или задает пару целых чисел, представляющих упорядоченную пару координат X и Y.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Реализует слияние двух случайных величин с заданными параметрами.
        /// </summary>
        /// <param name="x">Первая случайная величина</param>
        /// <param name="y">Вторая случайная величина</param>
        /// <param name="a0">Первый параметр</param>
        /// <param name="a1">Второй параметр</param>
        /// <returns>Целое число без знака</returns>
        private static byte merge(byte x, byte y, int a0, int a1)
        {
            return Maths.Byte((x * a0 + y * a1) / 255);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр ректификации.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Image_rectification
    /// </remarks>
    /// </summary>
    public class Rectification : IBitmapFilter2
    {
        #region Private data
        private double[,] H;
        private Color color;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр ректификации.
        /// </summary>
        /// <param name="homography">Матрица гомографии [3, 3]</param>
        /// <param name="color">Цвет фона</param>
        public Rectification(double[,] homography, Color color)
        {
            this.H = homography; Color = color;
        }
        /// <summary>
        /// Получает или задает матрицу гомографии [3, 3], используемую для сопоставления изображения, 
        /// переданного фильтру, на накладываемое изображение.
        /// </summary>
        public double[,] Homography
        {
            get
            {
                return this.H;
            }
            set
            {
                if (value.GetLength(0) != 3 && value.GetLength(1) != 3)
                    throw new Exception("Размерность матрицы должна быть [3, 3]");

                this.H = value;
            }
        }
        /// <summary>
        /// Получает или задает цвет фона.
        /// </summary>
        public Color Color
        {
            get { return color; }
            set { color = value; }
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            // get source image size
            int width = bmSrc.Width;
            int height = bmSrc.Height;
            int srcStride = bmSrc.Stride;

            // get destination image size
            int nWidth = bmData.Width;
            int nHeight = bmData.Height;
            int dstStride = bmData.Stride;

            // fill values
            byte R = color.R;
            byte G = color.G;
            byte B = color.B;
            byte A = color.A;

            // Retrieve homography matrix as float array
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();

            int y, x, ox, oy, c;
            double cx, cy, hw, hx, hy;

            // Project the second image
            for (y = 0; y < nHeight; y++)
            {
                for (x = 0; x < nWidth; x++, dst += 4)
                {
                    cx = x;
                    cy = y;

                    // projection using homogenous coordinates
                    hw = (H[2, 0] * cx + H[2, 1] * cy + H[2, 2]);
                    hx = (H[0, 0] * cx + H[0, 1] * cy + H[0, 2]) / hw;
                    hy = (H[1, 0] * cx + H[1, 1] * cy + H[1, 2]) / hw;

                    // coordinate of the nearest point
                    ox = (int)(hx);
                    oy = (int)(hy);

                    // validate source pixel's coordinates
                    if ((ox >= 0) && (oy >= 0) && (ox < width) && (oy < height))
                    {
                        c = oy * srcStride + ox * 4;

                        // 32bpp
                        dst[0] = src[c + 0];
                        dst[1] = src[c + 1];
                        dst[2] = src[c + 2];
                        dst[3] = src[c + 3];
                    }
                    else
                    {
                        // 32bpp
                        dst[0] = B;
                        dst[1] = G;
                        dst[2] = R;
                        dst[3] = A;
                    }
                }
            }

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтры текстуризации.
    /// </summary>
    public class Texturer : IBitmapFilter
    {
        #region Private data
        private double[,] texture = null;
        private double depth = 0.5;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтры текстуризации.
        /// </summary>
        /// <param name="texture">Матрица-текстура</param>
        public Texturer(double[,] texture)
        {
            Texture = texture;
        }
        /// <summary>
        /// Инициализирует фильтры текстуризации.
        /// </summary>
        /// <param name="texture">Матрица-текстура</param>
        /// <param name="depth">Глубина фильтра [0, 1]</param>
        public Texturer(double[,] texture, double depth = 1.0)
        {
            Texture = texture; Depth = depth;
        }
        /// <summary>
        /// Получает или задает матрицу-текстуру.
        /// </summary>
        public double[,] Texture
        {
            get { return texture; }
            set { texture = value; }
        }
        /// <summary>
        /// Получает или задает глубину фильтра [0, 1].
        /// </summary>
        public double Depth
        {
            get { return this.depth; }
            set { this.depth = Maths.Double(value); }
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Генератор случайных чисел.
        /// </summary>
        private static Random rand = new Random();
        /// <summary>
        /// Случайное значение.
        /// </summary>
        private static int r;
        /// <summary>
        /// Реализует построение текстуры дерева заданного размера.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="rings">Количество колец</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение текстуры ткани заданного размера.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение текстуры мрамора заданного размера.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="yPeriod">Период по высоте</param>
        /// <param name="xPeriod">Период по ширине</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение текстуры лабиринта заданного размера.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
        /// Реализует построение текстуры облаков заданного размера.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матрица</returns>
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
    /// Определяет фильтр линейной операции.
    /// <remarks>
    /// Данный фильтр работает по следующему алгоритму: C(x,y) = a * A(x,y) + b * B(x,y), где A, B - исходные изображения, 
    /// a, b - коэффициенты.
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
        /// Инициализирует фильтр линейной операции.
        /// </summary>
        /// <param name="a">Коэффициент первого изображения</param>
        /// <param name="b">Коэффициент второго изображения</param>
        public Operation(double a, double b)
        {
            A = a; B = b; 
        }
        /// <summary>
        /// Получает или задает коэффициент первого изображения.
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
        /// Получает или задает коэффициент второго изображения.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Возвращает фильтр операции сложения.
        /// </summary>
        public static Operation Addition
        {
            get
            {
                return new Operation(1, 1);
            }
        }
        /// <summary>
        /// Возвращает фильтр операции вычитания.
        /// </summary>
        public static Operation Subtraction
        {
            get
            {
                return new Operation(1, -1);
            }
        }
        /// <summary>
        /// Возвращает фильтр операции усреднения.
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
    /// Определяет фильтр создания стереоэффекта для пары изображений.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
    /// </remarks>
    /// </summary>
    public class StereoAnaglyph : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Алгоритм.
        /// </summary>
        protected Anaglyph algorithm = Anaglyph.Gray;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр создания стереоэффекта для пары изображений.
        /// </summary>
        /// <param name="algorithm">Алгоритм</param>
        public StereoAnaglyph(Anaglyph algorithm)
        {
            this.algorithm = algorithm;
        }
        /// <summary>
        /// Получает или задает алгоритм.
        /// </summary>
        public Anaglyph Algorithm
        {
            get { return algorithm; }
            set { algorithm = value; }
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения (левое изображение)</param>
        /// <param name="bmSrc">Атрибуты точечного изображения (правое изображение)</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Определяет алгоритм создания стереоэффекта.
        /// </summary>
        /// <remarks>
        /// Более подробную информацию можно найти на сайте:
        /// http://www.3dtv.at/Knowhow/AnaglyphComparison_en.aspx
        /// </remarks>
        public enum Anaglyph
        {
            /// <summary>
            /// Создает стереоэффект для пары изображений согласно следующим вычислениям:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            True,
            /// <summary>
            /// Создает стереоэффект для пары изображений согласно следующим вычислениям:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=0.299*R<sub>r</sub>+0.587*G<sub>r</sub>+0.114*B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Gray,
            /// <summary>
            /// Создает стереоэффект для пары изображений согласно следующим вычислениям:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=R<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            Color,
            /// <summary>
            /// Создает стереоэффект для пары изображений согласно следующим вычислениям:
            /// <list type="bullet">
            /// <item>R<sub>a</sub>=0.299*R<sub>l</sub>+0.587*G<sub>l</sub>+0.114*B<sub>l</sub>;</item>
            /// <item>G<sub>a</sub>=G<sub>r</sub>;</item>
            /// <item>B<sub>a</sub>=B<sub>r</sub>.</item>
            /// </list>
            /// </summary>
            HalfColor,
            /// <summary>
            /// Создает стереоэффект для пары изображений согласно следующим вычислениям:
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
    /// Определяет фильтр рисунка маслом.
    /// <remarks>
    /// Данный фильтр был портирован с языка C++.
    /// Подробное описание алгоритма можно найти на сайте:
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
        /// Инициализирует фильтр рисунка маслом.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="depth">Глубина [0, 1]</param>
        public OilPainting(int radius = 3, double depth = 1.0)
        {
            Size = new SizeInt(radius, radius);
            Depth = depth;
        }
        /// <summary>
        /// Инициализирует фильтр рисунка маслом.
        /// </summary>
        /// <param name="height">Радиус по высоте</param>
        /// <param name="width">Радиус по ширине</param>
        /// <param name="depth">Глубина [0, 1]</param>
        public OilPainting(int width, int height, double depth = 1.0)
        {
            Size = new SizeInt(width, height);
            Depth = depth;
        }
        /// <summary>
        /// Инициализирует фильтр рисунка маслом.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="depth">Глубина [0, 1]</param>
        public OilPainting(SizeInt size, double depth = 1.0)
        {
            Size = size;
            Depth = depth;
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Получает или задает значение глубины [0, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр на основе матрицы свертки.
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
        private bool twoside;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр на основе двумерной матрицы свертки.
        /// </summary>
        /// <param name="m">Матричный оператор</param>
        /// <param name="offset">Значение смещения</param>
        /// <param name="twoside">Двусторонняя обработка или нет</param>
        public Convolution(double[,] m, double offset = 0, bool twoside = false)
        {
            Matrix = m; Offset = offset; Twoside = twoside;
        }
        /// <summary>
        /// Инициализирует фильтр на основе двумерной матрицы свертки.
        /// </summary>
        public Convolution()
        {
            Matrix = Matrice.One(3, 3); Offset = 0; Twoside = false;
        }
        /// <summary>
        /// Получает или задает матричный оператор.
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
        /// Получает или задает значение смещения.
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
        /// Двусторонняя обработка или нет.
        /// </summary>
        public bool Twoside
        {
            get
            {
                return this.twoside;
            }
            set
            {
                this.twoside = value;
            }
        }
        /// <summary>
        /// Устанавливает значения параметров свертки.
        /// </summary>
        /// <param name="m">Матричный оператор</param>
        private void Data(double[,] m)
        {
            this.l0 = m.GetLength(0);
            this.l1 = m.GetLength(1);
            this.radius0 = l0 >> 1;
            this.radius1 = l1 >> 1;
            this.kernel = Jagged.ToJagged(m);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public unsafe void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            #region Data
            if (this.l0 != this.l1 && this.twoside == true)
                throw new Exception("Матрица должна быть квадратной");

            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            #endregion

            if (!this.twoside)
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Получает значение оператора градиента.
        /// </summary>
        /// <param name="Gx">Случайная величина</param>
        /// <param name="Gy">Случайная величина</param>
        /// <returns>Число двойной точности</returns>
        public static double G(double Gx, double Gy)
        {
            return Math.Sqrt(Gx * Gx + Gy * Gy);
        }
        /// <summary>
        /// Получает направление оператора градиента.
        /// </summary>
        /// <param name="Gx">Случайная величина</param>
        /// <param name="Gy">Случайная величина</param>
        /// <returns>Число двойной точности</returns>
        public static double Tetta(double Gx, double Gy)
        {
            return Math.Atan(Gx / Gy);
        }
        #endregion

        #region Public static methods
        #region Radius matrix
        /// <summary>
        /// Реализует построение перевернутого оператора Гаусса.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="sigma">Отклонение (!=0)</param>
        /// <returns>Матричный оператор</returns>
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
        /// Реализует построение оператора Гаусса.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="sigma">Отклонение (!=0)</param>
        /// <returns>Матричный оператор</returns>
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
        /// Реализует построение оператора "нерезкого маскирования".
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="sigma">Отклонение (!=0)</param>
        /// <returns>Матричный оператор</returns>
        public static Convolution Unsharp(int m, int l, double sigma)
        {
            // построение оператора Гаусса:
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

            // Вычисление матрицы:
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

            // Вычисление цетрального значения:
            invG[m / 2, l / 2] = 2 * summary - invG[m / 2, l / 2];
            return new Convolution(invG);
        }
        /// <summary>
        /// Реализует построение оператора выделения верхних частот.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <param name="boost">Усиление</param>
        /// <returns>Матричный оператор</returns>
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
        /// Реализует построение оператора выделения нижних частот.
        /// </summary>
        /// <param name="m">Высота</param>
        /// <param name="l">Ширина</param>
        /// <returns>Матричный оператор</returns>
        public static Convolution LowPass(int m, int l)
        {
            return Convolution.HighPass(m, l, 1);
        }
        /// <summary>
        /// Реализует построение оператора выделения рельефа.
        /// </summary>
        /// <param name="n">Размер матрицы</param>
        /// <returns>Матричный оператор</returns>
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
        /// Реализует построение оператора Робертса [2 x 2].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution Roberts()
        {
            return new Convolution(new double[2, 2] { { 1, 0 }, { 0, -1 } });
        }
        /// <summary>
        /// Реализует построение оператора Прюитта [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution Prewitt()
        {
            return new Convolution(new double[3, 3] { { -1, -1, -1 }, { 0, 0, 0 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Реализует построение оператора Собеля [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution Sobel()
        {
            return new Convolution(new double[3, 3] { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } });
        }
        /// <summary>
        /// Реализует построение оператора Щарра [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution Scharr()
        {
            return new Convolution(new double[3, 3] { { 3, 10, 3 }, { 0, 0, 0 }, { -3, -10, -3 } });
        }
        /// <summary>
        /// Реализует построение оператора Лапласа [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution Laplacian()
        {
            return new Convolution(new double[3, 3] { { 0, 1, 0 }, { 1, -4, 1 }, { 0, 1, 0 } });
        }
        /// <summary>
        /// Реализует построение диагонального оператора Лапласа [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution LaplacianDiagonal()
        {
            return new Convolution(new double[3, 3] { { 1, 1, 1 }, { 1, -8, 1 }, { 1, 1, 1 } });
        }
        /// <summary>
        /// Реализует построение обратного оператора Лапласа [3 x 3].
        /// </summary>
        /// <returns>Матричный оператор</returns>
        public static Convolution LaplacianInvert()
        {
            return new Convolution(new double[3, 3] { { -1, 0, -1 }, { 0, 4, 0 }, { -1, 0, -1 } });
        }
        #endregion

        #region Fixed radius compass matrix
        /// <summary>
        /// Реализует построение оператора-компаса Кирша [3 x 3].
        /// </summary>
        /// <param name="direction">Направление оператора градиента</param>
        /// <returns>Матричный оператор</returns>
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
        /// Реализует построение оператора-компаса Робертса [2 x 2].
        /// </summary>
        /// <param name="direction">Направление оператора градиента</param>
        /// <returns>Матричный оператор</returns>
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
    /// Определяет фильтр морфологических операций.
    /// </summary>
    public class Morphology : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Пороговое значение.
        /// </summary>
        private int threshold;
        /// <summary>
        /// Ширина фильтра.
        /// </summary>
        private int rw;
        /// <summary>
        /// Высота фильтра.
        /// </summary>
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="threshold">Пороговое значение</param>
        public Morphology(int radius = 3, int threshold = 0)
        {
            Size = new SizeInt(radius, radius);
            Threshold = threshold;
        }
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="threshold">Пороговое значение</param>
        public Morphology(int width, int height, int threshold)
        {
            Size = new SizeInt(width, height);
            Threshold = threshold;
        }
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        /// <param name="threshold">Пороговое значение</param>
        public Morphology(SizeInt size, int threshold = 0)
        {
            Size = size;
            Threshold = threshold;
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Получает или задает пороговое значение.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Инициализирует фильтр медианы.
        /// </summary>
        /// <param name="radius">Радиус</param>
        public static Morphology Median(int radius)
        {
            int threshold = radius / 2;
            return new Morphology(radius, radius, threshold);
        }
        /// <summary>
        /// Инициализирует фильтр эрозии.
        /// </summary>
        /// <param name="radius">Радиус</param>
        public static Morphology Erosion(int radius)
        {
            return new Morphology(radius, radius, 0);
        }
        /// <summary>
        /// Инициализирует фильтр расширения.
        /// </summary>
        /// <param name="radius">Радиус</param>
        public static Morphology Dilatation(int radius)
        {
            return new Morphology(radius, radius, radius - 1);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр расширения.
    /// </summary>
    public class Dilatation : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Ширина фильтра.
        /// </summary>
        private int rw;
        /// <summary>
        /// Высота фильтра.
        /// </summary>
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр расширения.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public Dilatation(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Инициализирует фильтр расширения.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public Dilatation(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр расширения.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public Dilatation(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
    /// Определяет фильтр эрозии.
    /// </summary>
    public class Erosion : IBitmapFilter2
    {
        #region Private data
        /// <summary>
        /// Ширина фильтра.
        /// </summary>
        private int rw;
        /// <summary>
        /// Высота фильтра.
        /// </summary>
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр эрозии.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public Erosion(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Инициализирует фильтр эрозии.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public Erosion(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр эрозии.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public Erosion(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
    /// Определяет фильтр Top-Hat.
    /// </summary>
    public class TopHat : IBitmapFilter2
    {
        #region Private data
        private Opening opening = new Opening();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр Top-Hat.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public TopHat(int radius = 3)
        {
            opening = new Opening(radius);
        }
        /// <summary>
        /// Инициализирует фильтр Top-Hat.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public TopHat(int width, int height)
        {
            opening = new Opening(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр Top-Hat.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public TopHat(SizeInt size)
        {
            opening = new Opening(size);
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр Bottom-Hat.
    /// </summary>
    public class BottomHat : IBitmapFilter2
    {
        #region Private data
        private Closing closing = new Closing();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр Bottom-Hat.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public BottomHat(int radius = 3)
        {
            closing = new Closing(radius);
        }
        /// <summary>
        /// Инициализирует фильтр Bottom-Hat.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public BottomHat(int width, int height)
        {
            closing = new Closing(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр Bottom-Hat.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public BottomHat(SizeInt size)
        {
            closing = new Closing(size);
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр закрытия.
    /// </summary>
    public class Closing : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Dilatation dilatation = new Dilatation();
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр закрытия.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public Closing(int radius = 3)
        {
            erosion = new Erosion(radius);
            dilatation = new Dilatation(radius);
        }
        /// <summary>
        /// Инициализирует фильтр закрытия.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public Closing(int width, int height)
        {
            erosion = new Erosion(width, height);
            dilatation = new Dilatation(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр закрытия.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public Closing(SizeInt size)
        {
            erosion = new Erosion(size);
            dilatation = new Dilatation(size);
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            dilatation.Apply(bmSrc, bmData);
            erosion.Apply(bmData, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр открытия.
    /// </summary>
    public class Opening : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Dilatation dilatation = new Dilatation();
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр открытия.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public Opening(int radius = 3)
        {
            erosion = new Erosion(radius);
            dilatation = new Dilatation(radius);
        }
        /// <summary>
        /// Инициализирует фильтр открытия.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public Opening(int width, int height)
        {
            erosion = new Erosion(width, height);
            dilatation = new Dilatation(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр открытия.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public Opening(SizeInt size)
        {
            erosion = new Erosion(size);
            dilatation = new Dilatation(size);
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData, BitmapData bmSrc)
        {
            erosion.Apply(bmSrc, bmData);
            dilatation.Apply(bmData, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр свечения краев.
    /// </summary>
    public class EdgeGlow : IBitmapFilter2
    {
        #region Private data
        private Erosion erosion = new Erosion();
        private Operation subtraction = Operation.Subtraction;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр свечения краев.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public EdgeGlow(int radius = 3)
        {
            erosion = new Erosion(radius);
        }
        /// <summary>
        /// Инициализирует фильтр свечения краев.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public EdgeGlow(int width, int height)
        {
            erosion = new Erosion(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр свечения краев.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public EdgeGlow(SizeInt size)
        {
            erosion = new Erosion(size);
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
        public void Apply(Bitmap Data, Bitmap Src)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmSrc = BitmapConverter.Lock32bpp(Src);
            Apply(bmData, bmSrc);
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Src, bmSrc);
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр квадратного размытия.
    /// </summary>
    public class BoxBlur : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Ширина фильтра.
        /// </summary>
        private int rw;
        /// <summary>
        /// Высота фильтра.
        /// </summary>
        private int rh;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        public BoxBlur(int radius = 3)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public BoxBlur(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр квадратного размытия.
        /// </summary>
        /// <param name="size">Размеры фильтра</param>
        public BoxBlur(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Получает или задает размеры фильтра.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
    /// Определяет фильтр на основе структуры HSB.
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
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="hue">Оттенок [0, 359]</param>
        /// <param name="saturation">Насыщенность [-1, 1]</param>
        /// <param name="brightness">Яркость [-1, 1]</param>
        public HSBFilter(int hue, double saturation, double brightness)
        {
            Hue = hue;
            Saturation = saturation;
            Brightness = brightness;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        public HSBFilter()
        {
            new HSBFilter(0, 0, 0);
        }
        /// <summary>
        /// Оттенок [0, 359].
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
        /// Насыщенность [-1, 1].
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
        /// Яркость [-1, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр на основе структуры HSL.
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
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        /// <param name="hue">Оттенок [0, 359]</param>
        /// <param name="saturation">Насыщенность [-1, 1]</param>
        /// <param name="lightness">Световая интенсивность [-1, 1]</param>
        public HSLFilter(int hue, double saturation, double lightness)
        {
            Hue = hue;
            Saturation = saturation;
            Lightness = lightness;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры HSB.
        /// </summary>
        public HSLFilter()
        {
            new HSLFilter(0, 0, 0);
        }
        /// <summary>
        /// Оттенок [0, 359].
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
        /// Насыщенность [-1, 1].
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
        /// Световая интенсивность [-1, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            Apply(bmData);
            BitmapConverter.Unlock(Data, bmData);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр на основе структуры YCbCr.
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
        /// Инициализирует фильтр на основе структуры YCbCr.
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
        /// Инициализирует фильтр на основе структуры YCbCr.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр на основе структуры CMYK.
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
        /// Инициализирует фильтр на основе структуры CMYK.
        /// </summary>
        /// <param name="cyan">Голубой [-1, 1]</param>
        /// <param name="magenta">Пурпурный [-1, 1]</param>
        /// <param name="yellow">Желтый [-1, 1]</param>
        /// <param name="keycolor">Черный [-1, 1]</param>
        public CMYKFilter(double cyan, double magenta, double yellow, double keycolor)
        {
            Cyan = cyan;
            Magenta = magenta;
            Yellow = yellow;
            Keycolor = keycolor;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры CMYK.
        /// </summary>
        public CMYKFilter()
        {
            new CMYKFilter(0, 0, 0, 0);
        }
        /// <summary>
        /// Голубой [-1, 1].
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
        /// Пурпурный [-1, 1].
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
        /// Желтый [-1, 1].
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
        /// Черный [-1, 1].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр на основе структуры RGB.
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
        /// Инициализирует фильтр на основе структуры RGB.
        /// </summary>
        /// <param name="red">Красный [-255, 255]</param>
        /// <param name="green">Зеленый [-255, 255]</param>
        /// <param name="blue">Синий [-255, 255]</param>
        public RGBFilter(int red, int green, int blue)
        {
            Red = red;
            Green = green;
            Blue = blue;
        }
        /// <summary>
        /// Инициализирует фильтр на основе структуры RGB.
        /// </summary>
        public RGBFilter()
        {
            new RGBFilter(0, 0, 0);
        }
        /// <summary>
        /// Красный [-255, 255].
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
        /// Зеленый [-255, 255].
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
        /// Синий [-255, 255].
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр глобальной эквализации гистограммы.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
    /// http://www.cromwell-intl.com/3d/histogram/
    /// </remarks>
    /// </summary>
    public class HistogramEqualization : IBitmapFilter
    {
        #region Filter components
        /// <summary>
        /// Инициализирует фильтр глобальной эквализации гистограммы.
        /// </summary>
        public HistogramEqualization() { }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр локальной эквализации гистограммы.
    /// <remarks>
    /// Подробное описание алгоритма можно найти на сайте:
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
        /// Инициализирует фильтр локальной эквализации гистограммы.
        /// </summary>
        /// <param name="radius">Радиус матрицы</param>
        public LocalHistogramEqualization(int radius = 10)
        {
            Size = new SizeInt(radius, radius);
        }
        /// <summary>
        /// Инициализирует фильтр локальной эквализации гистограммы.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        public LocalHistogramEqualization(int width, int height)
        {
            Size = new SizeInt(width, height);
        }
        /// <summary>
        /// Инициализирует фильтр локальной эквализации гистограммы.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        public LocalHistogramEqualization(SizeInt size)
        {
            Size = size;
        }
        /// <summary>
        /// Получает или задает размер фильтра.
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
        /// Устанавливает значения параметров свертки.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        private void Data(SizeInt size)
        {
            this.l0 = size.Width;
            this.l1 = size.Height;
            this.rw = l0 >> 1;
            this.rh = l1 >> 1;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmSrc">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="Src">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
    /// Определяет фильтр глобального сжатия гистограммы.
    /// </summary>
    public class HistogramStretch : Correction, IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// { Min, Max }.
        /// </summary>
        protected RangeDouble range;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр глобального сжатия гистограммы.
        /// </summary>
        /// <param name="min">Минимальное значение яркости [0, 1]</param>
        /// <param name="max">Максимально значение яркости [0, 1]</param>
        /// <param name="space">Цветовое пространство</param>
        public HistogramStretch(double min, double max, Space space)
        {
            Range = new RangeDouble(min, max);
            Space = space;
        }
        /// <summary>
        /// Инициализирует фильтр глобального сжатия гистограммы.
        /// </summary>
        /// <param name="range">Диапазон значений яркости</param>
        /// <param name="space">Цветовое пространство</param>
        public HistogramStretch(RangeDouble range, Space space)
        {
            Range = range;
            Space = space;
        }
        /// <summary>
        /// Получает или задает диапазон значений, которые принимает сигнал.
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
        /// Реализует перестроение данных фильтра.
        /// </summary>
        protected override void Rebuild()
        {
            this.values = Intensity.Equalize(range.Min, range.Max, 256);
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс фильтра локального сжатия гистограммы.
    /// <remarks>
    /// Фильтр используется для коррекции неравномерно освещенных изображений и повышения детализации.
    /// Более подробную информацию можно найти на сайте:
    /// http://www.academia.edu/7629047/Image_enhancement_by_local_histogram_stretching
    /// </remarks>
    /// </summary>
    public class LocalHistogramStretch : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Фильтр локального усреднения.
        /// </summary>
        private BoxBlur gb = new BoxBlur();
        /// <summary>
        /// Фильтр эрозии.
        /// </summary>
        private Erosion er = new Erosion();
        /// <summary>
        /// Фильтр приращения.
        /// </summary>
        private Dilatation di = new Dilatation();
        /// <summary>
        /// Цветовое пространство.
        /// </summary>
        protected Space space;
        /// <summary>
        /// Коэффициент сжатия контраста.
        /// </summary>
        protected double contrast;
        /// <summary>
        /// Сглаживание.
        /// </summary>
        protected bool smoothing;
        #endregion

        #region Filter componetns
        /// <summary>
        /// Иницилазирует класс фильтра локального сжатия гистограммы.
        /// </summary>
        /// <param name="radius">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [0, 1]</param>
        /// <param name="smoothing">Использовать сглаживание или нет</param>
        public LocalHistogramStretch(int radius, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = new SizeInt(radius, radius);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Иницилазирует класс фильтра локального сжатия гистограммы.
        /// </summary>
        /// <param name="width">Ширина фильтра</param>
        /// <param name="height">Высота фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [0, 1]</param>
        /// <param name="smoothing">Использовать сглаживание или нет</param>
        public LocalHistogramStretch(int width, int height, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = new SizeInt(width, height);
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Иницилазирует класс фильтра локального сжатия гистограммы.
        /// </summary>
        /// <param name="size">Размер фильтра</param>
        /// <param name="space">Цветовое пространство</param>
        /// <param name="contrast">Коэффициент сжатия контраста [0, 1]</param>
        /// <param name="smoothing">Использовать сглаживание или нет</param>
        public LocalHistogramStretch(SizeInt size, Space space, double contrast = 0.5, bool smoothing = true)
        {
            Size = size;
            Space = space;
            Smoothing = smoothing;
            Contrast = contrast;
        }
        /// <summary>
        /// Получает или задает размер фильтра.
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
        /// Получает или задает цветовое пространство.
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
        /// Получает или задает значение коэффициента сжатия контраста [0, 1].
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
        /// Использовать сглаживание или нет.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData)
        {
            // Создание копий изображения:
            Bitmap Max = (Bitmap)BitmapConverter.Bitmap(bmData).Clone();
            Bitmap Min = (Bitmap)Max.Clone();

            // Применение фильтров морфологии:
            di.Apply(Max); er.Apply(Min);

            // Получение атрибутов:
            BitmapData bmMax = BitmapConverter.Lock32bpp(Max);
            BitmapData bmMin = BitmapConverter.Lock32bpp(Min);

            // Сглаживание ФНЧ:
            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            // Применение фильтра сжатия гистограммы:
            Apply(bmData, bmMax, bmMin);

            // Высвобождение данных:
            BitmapConverter.Unlock(Max, bmMax);
            BitmapConverter.Unlock(Min, bmMin);

            // Удаление объектов:
            Max.Dispose(); Min.Dispose();
            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        public void Apply(Bitmap Data)
        {
            // Создание копий изображения:
            Bitmap Max = (Bitmap)Data.Clone();
            Bitmap Min = (Bitmap)Data.Clone();

            // Применение фильтров морфологии:
            di.Apply(Max); er.Apply(Min);

            // Получение атрибутов:
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            BitmapData bmMax = BitmapConverter.Lock32bpp(Max);
            BitmapData bmMin = BitmapConverter.Lock32bpp(Min);

            // Сглаживание ФНЧ:
            if (smoothing)
                gb.Apply(bmMax); gb.Apply(bmMin);

            // Применение фильтра сжатия гистограммы:
            Apply(bmData, bmMax, bmMin);

            // Высвобождение данных:
            BitmapConverter.Unlock(Data, bmData);
            BitmapConverter.Unlock(Max, bmMax);
            BitmapConverter.Unlock(Min, bmMin);

            // Удаление объектов:
            Max.Dispose(); Min.Dispose();
            return;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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
                        // Смещение по пикселям:
                        v = k + q;

                        // Получение яркостных компонент:
                        mag = p[v] / 255.0;
                        max = pMax[v] / 255.0;
                        min = pMin[v] / 255.0;

                        // Вычисление функции контраста:
                        num1 = max - min;

                        if (num1 < required)
                        {
                            num2 = min + (required - num1) * min / (num1 - 1f);
                            min = Maths.Double(num2);
                            max = Maths.Double(num2 + required);
                        }

                        // Локальное сжатие гистограммы:
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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

                    // Преобразование в другое цветовое пространство:
                    imag = HSB.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSB.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSB.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    // Получение яркостных компонент:
                    mag = imag.Brightness;
                    max = imax.Brightness;
                    min = imin.Brightness;

                    // Вычисление функции контраста:
                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    // Локальное сжатие гистограммы:
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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

                    // Преобразование в другое цветовое пространство:
                    imag = HSL.FromRGB(p[k2], p[k1], p[k]);
                    imax = HSL.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = HSL.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    // Получение яркостных компонент:
                    mag = imag.Lightness;
                    max = imax.Lightness;
                    min = imin.Lightness;

                    // Вычисление функции контраста:
                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    // Локальное сжатие гистограммы:
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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

                    // Преобразование в другое цветовое пространство:
                    imag = YCbCr.FromRGB(p[k2], p[k1], p[k]);
                    imax = YCbCr.FromRGB(pMax[k2], pMax[k1], pMax[k]);
                    imin = YCbCr.FromRGB(pMin[k2], pMin[k1], pMin[k]);

                    // Получение яркостных компонент:
                    mag = imag.Y;
                    max = imax.Y;
                    min = imin.Y;

                    // Вычисление функции контраста:
                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    // Локальное сжатие гистограммы:
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="bmMax">Атрибуты точечного изображения</param>
        /// <param name="bmMin">Атрибуты точечного изображения</param> 
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

                    // Получение яркостных компонент:
                    mag = RGB.Average(p[v], p[v + 1], p[v + 2]) / 255.0;
                    max = RGB.Average(pMax[v], pMax[v + 1], pMax[v + 2]) / 255.0;
                    min = RGB.Average(pMin[v], pMin[v + 1], pMin[v + 2]) / 255.0;

                    // Вычисление функции контраста:
                    num1 = max - min;

                    if (num1 < required)
                    {
                        num2 = min + (required - num1) * min / (num1 - 1f);
                        min = Maths.Double(num2);
                        max = Maths.Double(num2 + required);
                    }

                    // Локальное сжатие гистограммы:
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
    /// Определяет фильтр обработки изображений.
    /// </summary>
    public class BitmapFilter : IBitmapFilter
    {
        #region Private data
        /// <summary>
        /// Направленный фильтр.
        /// </summary>
        protected IFilter filter;
        /// <summary>
        /// Цветовое пространство.
        /// </summary>
        protected Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр обработки изображений.
        /// </summary>
        /// <param name="filter">Фильтр</param>
        /// <param name="space">Цветовое пространство</param>
        public BitmapFilter(IFilter filter, Space space = Space.RGB)
        {
            this.filter = filter;
            this.space = space;
        }
        /// <summary>
        /// Получает или задает фильтр.
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
        /// Получает или задает цветовое пространство.
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        private unsafe void ApplyRGB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] r = new double[width, height];
            double[,] g = new double[width, height];
            double[,] b = new double[width, height];

            // Получение яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    r[i, j] = p[k + 2] / 255.0;
                    g[i, j] = p[k + 1] / 255.0;
                    b[i, j] = p[k] / 255.0;
                }
            }
            );

            // Применение фильтра:
            this.filter.Apply(r);
            this.filter.Apply(g);
            this.filter.Apply(b);

            // Указание новых яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k + 2] = Maths.Byte(r[i, j] * 255);
                    p[k + 1] = Maths.Byte(g[i, j] * 255);
                    p[k] = Maths.Byte(b[i, j] * 255);
                }
            }
            );

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        private unsafe void ApplyHSB(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] input = new double[width, height];

            // Получение яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                HSB lumI;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    lumI = HSB.FromRGB(p[k + 2], p[k + 1], p[k]);
                    input[i, j] = lumI.Brightness;
                }
            }
            );

            // Применение фильтра:
            this.filter.Apply(input);

            // Указание новых яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                HSB lumI; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    lumI = HSB.FromRGB(p[k + 2], p[k + 1], p[k]);
                    lumI.Brightness = input[i, j];

                    rgb = lumI.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        private unsafe void ApplyHSL(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] input = new double[width, height];

            // Получение яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                HSL lumI;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    lumI = HSL.FromRGB(p[k + 2], p[k + 1], p[k]);
                    input[i, j] = lumI.Lightness;
                }
            }
            );

            // Применение фильтра:
            this.filter.Apply(input);

            // Указание новых яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                HSL lumI; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    lumI = HSL.FromRGB(p[k + 2], p[k + 1], p[k]);
                    lumI.Lightness = input[i, j];

                    rgb = lumI.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        private unsafe void ApplyYCbCr(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] input = new double[width, height];

            // Получение яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                YCbCr lumI;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    lumI = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k]);
                    input[i, j] = lumI.Y;
                }
            }
            );

            // Применение фильтра:
            this.filter.Apply(input);

            // Указание новых яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                YCbCr lumI; RGB rgb;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;

                    lumI = YCbCr.FromRGB(p[k + 2], p[k + 1], p[k]);
                    lumI.Y = input[i, j];

                    rgb = lumI.ToRGB;
                    p[k + 2] = rgb.Red; p[k + 1] = rgb.Green; p[k] = rgb.Blue;
                }
            }
            );

            return;
        }
        /// <summary>
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        private unsafe void ApplyGrayscale(BitmapData bmData)
        {
            byte* p = (byte*)bmData.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            double[,] r = new double[width, height];

            // Получение яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    r[i, j] = RGB.Average(p[k], p[k + 1], p[k + 2]) / 255.0;
                }
            }
            );

            // Применение фильтра:
            this.filter.Apply(r);

            // Указание новых яркостных характеристик изображения:
            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    k = jstride + i * 4;
                    p[k] = p[k + 1] = p[k + 2] = Maths.Byte(r[i, j] * 255);
                }
            }
            );

            return;
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр смешивания изображений.
    /// <remarks>
    /// Данный фильтр усредняет несколько входных изображений одинакового размера, получая тем самым единственное
    /// выходное изображение. Фильтр может использоваться для реализации эффектов HDR.
    /// Более подробную информацию можно найти на сайте:
    /// https://web.stanford.edu/class/cs231m/lectures/lecture-5-stitching-blending.pdf
    /// (стр. 65-75)
    /// </remarks>
    /// </summary>
    public class BitmapBlender
    {
        #region Private data
        /// <summary>
        /// Фильтр смешивания.
        /// </summary>
        private IBlendFilter filter;
        /// <summary>
        /// Цветовое пространство.
        /// </summary>
        private Space space;
        #endregion

        #region Filter components
        /// <summary>
        /// Инициализирует фильтр смешивания изображений.
        /// </summary>
        /// <param name="filter">Фильтр смешивания</param>
        /// <param name="space">Цветовое пространство</param>
        public BitmapBlender(IBlendFilter filter, Space space = Space.RGB)
        {
            this.filter = filter;
            this.space = space;
        }
        /// <summary>
        /// Получает или задает фильтр смешивания.
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
        /// Получает или задает цветовое пространство.
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="images">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="images">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="Data">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву атрибутов точечных рисунков.
        /// </summary>
        /// <param name="bmData">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="Data">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
        private Bitmap ApplyGrayscale(Bitmap[] Data)
        {
            int length = Data.Length;
            double[][,] r = new double[length][,];

            for (int i = 0; i < length; i++)
            {
                r[i] = BitmapConverter.FromBitmap(Data[i]);
            }

            return BitmapConverter.ToBitmap(this.filter.Apply(r));
        }
        /// <summary>
        /// Применяет фильтр к массиву атрибутов точечных рисунков.
        /// </summary>
        /// <param name="bmData">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
        private Bitmap ApplyGrayscale(BitmapData[] bmData)
        {
            int length = bmData.Length;
            double[][,] r = new double[length][,];

            for (int i = 0; i < length; i++)
            {
                r[i] = BitmapConverter.FromBitmap(bmData[i]);
            }

            return BitmapConverter.ToBitmap(this.filter.Apply(r));
        }
        /// <summary>
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="Data">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву атрибутов точечных рисунков.
        /// </summary>
        /// <param name="bmData">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="Data">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву атрибутов точечных рисунков.
        /// </summary>
        /// <param name="bmData">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву точечных рисунков.
        /// </summary>
        /// <param name="Data">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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
        /// Применяет фильтр к массиву атрибутов точечных рисунков.
        /// </summary>
        /// <param name="bmData">Массив точечных рисунков</param>
        /// <returns>Точечный рисунок</returns>
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

    #region Hough transforms
    /// <summary>
    /// Определяет линию Хафа.
    /// </summary>
    public struct HoughLine : IComparable
    {
        #region Struct components
        /// <summary>
        /// Наклон линии - угол между полярной осью и радиусом линии θ ∈ [0, 180).
        /// </summary>
        public readonly double Theta;
        /// <summary>
        /// Расстояние линии от центра изображения (-inf, +inf).
        /// <remarks>
        /// Отрицательный радиус линии означает, что линия находится в нижней части системы полярных координат. Поэтому 
        /// угол θ должен быть увеличен на 180 градусов, а радиус должен быть положительным.
        /// </remarks>
        /// </summary>
        public readonly short Radius;
        /// <summary>
        /// Абсолютная интенсивность линии (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Относительная интенсивность линии (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Инициализирует линию Хафа.
        /// </summary>
        /// <param name="theta">Наклон линии θ ∈ [0, 180)</param>
        /// <param name="radius">Расстояние линии от центра изображения (-inf, +inf)</param>
        /// <param name="intensity">Абсолютная интенсивность линии (0, +inf)</param>
        /// <param name="relativeIntensity">Относительная интенсивность линии (0, 1]</param>
        public HoughLine(double theta, short radius, short intensity, double relativeIntensity)
        {
            Theta = theta;
            Radius = radius;
            Intensity = intensity;
            RelativeIntensity = relativeIntensity;
        }
        /// <summary>
        /// Сравнивает объект с другим экземпляром этого класса.
        /// </summary>
        /// <param name="value">Объект</param>
        /// <returns>Целое число со знаком</returns>
        public int CompareTo(object value)
        {
            return (-Intensity.CompareTo(((HoughLine)value).Intensity));
        }
        #endregion
    }
    /// <summary>
    /// Определяет окружность Хафа.
    /// </summary>
    public struct HoughCircle : IComparable
    {
        #region Struct components
        /// <summary>
        /// Координата X.
        /// </summary>
        public readonly int X;
        /// <summary>
        /// Координата Y.
        /// </summary>
        public readonly int Y;
        /// <summary>
        /// Радиус окружности.
        /// </summary>
        public readonly int Radius;
        /// <summary>
        /// Абсолютная интенсивность линии (0, +inf).
        /// </summary>
        public readonly short Intensity;
        /// <summary>
        /// Относительная интенсивность линии (0, 1].
        /// </summary>
        public readonly double RelativeIntensity;
        /// <summary>
        /// Инициализирует окружность Хафа.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        /// <param name="radius">Радиус окружности</param>
        /// <param name="intensity">Абсолютная интенсивность линии (0, +inf)</param>
        /// <param name="relativeIntensity">Относительная интенсивность линии (0, 1]</param>
        public HoughCircle(int x, int y, int radius, short intensity, double relativeIntensity)
        {
            X = x;
            Y = y;
            Radius = radius;
            Intensity = intensity;
            RelativeIntensity = relativeIntensity;
        }
        /// <summary>
        /// Сравнивает объект с другим экземпляром этого класса.
        /// </summary>
        /// <param name="value">Объект</param>
        /// <returns>Целое число со знаком</returns>
        public int CompareTo(object value)
        {
            return (-Intensity.CompareTo(((HoughCircle)value).Intensity));
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр преобразования Хафа.
    /// </summary>
    public class HoughLineTransform : IBitmapFilter
    {
        #region Private data
        // Hough transformation quality settings
        private int stepsPerDegree;
        private int houghHeight;
        private double thetaStep;

        // precalculated Sine and Cosine values
        private double[] sinMap;
        private double[] cosMap;
        // Hough map
        private short[,] houghMap;
        private short maxMapIntensity = 0;

        private int localPeakRadius = 4;
        private short minLineIntensity = 10;
        private ArrayList lines = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Иницилизирует фильтр преобразования Хафа.
        /// </summary>
        public HoughLineTransform()
        {
            StepsPerDegree = 1;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Шаг в грудасах.
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
        /// Получает минимальную интенсивность линий.
        /// </summary>
        public short MinLineIntensity
        {
            get { return minLineIntensity; }
            set { minLineIntensity = value; }
        }
        /// <summary>
        /// Получает максимальную интенсивность линий.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
        }
        /// <summary>
        /// Получает или задает радиус поиска локального пикового значения.
        /// </summary>
        public int LocalPeakRadius
        {
            get { return localPeakRadius; }
            set { localPeakRadius = Math.Max(1, Math.Min(10, value)); }
        }
        /// <summary>
        /// Получает количество найденных линий.
        /// </summary>
        public int LinesCount
        {
            get { return lines.Count; }
        }
        /// <summary>
        /// Возвращает массив линий с абсолютной интенсивностью.
        /// </summary>
        /// <param name="count">Количество</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает массив линий с относительной интенсивностью.
        /// </summary>
        /// <param name="minRelativeIntensity">Минимальная относительная интесивность</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает синограмму Хафа.
        /// </summary>
        /// <returns>Точечный рисунок</returns>
        public Bitmap ToBitmap()
        {
            // check if Hough transformation was made already
            if (houghMap == null)
            {
                throw new ApplicationException("Преобразование Хафа не было произведено");
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format8bppIndexed)
            {
                throw new Exception("Формат изображения должен быть: 8bppIndexed");
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
    /// Определяет фильтр преобразования Хафа.
    /// </summary>
    public class HoughCircleTransform : IBitmapFilter
    {
        #region Private data
        // circle radius to detect
        private int radiusToDetect;

        // Hough map
        private short[,] houghMap;
        private short maxMapIntensity = 0;

        // Hough map's width and height
        private int width;
        private int height;

        private int localPeakRadius = 4;
        private short minCircleIntensity = 10;
        private ArrayList circles = new ArrayList();
        #endregion

        #region Initialize
        /// <summary>
        /// Иницилизирует фильтр преобразования Хафа.
        /// </summary>
        /// <param name="radiusToDetect">Радиус поиска</param>
        public HoughCircleTransform(int radiusToDetect)
        {
            this.radiusToDetect = radiusToDetect;
        }
        #endregion

        #region Class components
        /// <summary>
        /// Получает или задает минимальную интенсивность окружности.
        /// </summary>
        public short MinCircleIntensity
        {
            get { return minCircleIntensity; }
            set { minCircleIntensity = value; }
        }
        /// <summary>
        /// Получает или задает радиус поиска локального пикового значения.
        /// </summary>
        public int LocalPeakRadius
        {
            get { return localPeakRadius; }
            set { localPeakRadius = Math.Max(1, Math.Min(10, value)); }
        }
        /// <summary>
        /// Получает максимальную интенсивность линий.
        /// </summary>
        public short MaxIntensity
        {
            get { return maxMapIntensity; }
        }
        /// <summary>
        /// Получает количество найденных окружностей.
        /// </summary>
        public int CirclesCount
        {
            get { return circles.Count; }
        }
        /// <summary>
        /// Возвращает массив окружностей с абсолютной интенсивностью.
        /// </summary>
        /// <param name="count">Количество</param>
        /// <returns>Одномерный массив</returns>
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
        /// Возвращает массив окружностей с относительной интенсивностью.
        /// </summary>
        /// <param name="minRelativeIntensity">Минимальная относительная интесивность</param>
        /// <returns>Одномерный массив</returns>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
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
        /// Применяет фильтр к точечному рисунку.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        public void Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format8bppIndexed)
            {
                throw new Exception("Формат изображения должен быть: 8bppIndexed");
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
        /// Возвращает синограмму Хафа.
        /// </summary>
        /// <returns>Точечный рисунок</returns>
        public Bitmap ToBitmap()
        {
            // check if Hough transformation was made already
            if (houghMap == null)
            {
                throw new ApplicationException("Преобразование Хафа не было произведено");
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

    #region Bitmap conversions
    /// <summary>
    /// Используется для работы с точечными изображениями.
    /// </summary>
    public static class BitmapConverter
    {
        #region Bitmap convert components
        /// <summary>
        /// Конвертирует точечный рисунок в файл значка.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <param name="size">Размер (ширина = высота)</param>
        /// <returns>Файл значка</returns>
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
        /// Конвертирует точечный рисунок в JPEG формат.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap ToJpeg(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Jpeg);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Конвертирует точечный рисунок в BMP формат.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap ToBmp(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Bmp);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Конвертирует точечный рисунок в GIF формат.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap ToGif(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Gif);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Конвертирует точечный рисунок в Png формат.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap ToPng(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Png);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Конвертирует точечный рисунок в Tiff формат.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap ToTiff(this Bitmap b)
        {
            Bitmap bmp = new Bitmap(b);
            MemoryStream stream = new MemoryStream();
            bmp.Save(stream, ImageFormat.Tiff);

            return new Bitmap(stream);
        }
        /// <summary>
        /// Получает точечный рисунок из атрибутов точечного изображения.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap Bitmap(this BitmapData bmData)
        {
            return new Bitmap(bmData.Width, bmData.Height, bmData.Stride, bmData.PixelFormat, bmData.Scan0);
        }
        /// <summary>
        /// Конвертирует точечный рисунок в определенный формат
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <param name="pixelformat">Формат данных о цвете для каждой точки точечного рисунка</param>
        /// <returns>Точечный рисунок</returns>
        public static Bitmap Bitmap(this Bitmap b, PixelFormat pixelformat)
        {
            return b.Clone(new Rectangle(0, 0, b.Width, b.Height), pixelformat);
        }
        #endregion

        #region BitmapData voids
        /// <summary>
        /// Блокирует точечный рисунок в системной памяти.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Атрибуты точечного изображения</returns>
        public static BitmapData Lock32bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format32bppArgb);
        }
        /// <summary>
        /// Блокирует точечный рисунок в системной памяти.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <returns>Атрибуты точечного изображения</returns>
        public static BitmapData Lock8bpp(this Bitmap b)
        {
            return b.LockBits(new Rectangle(0, 0, b.Width, b.Height), ImageLockMode.ReadWrite, PixelFormat.Format8bppIndexed);
        }
        /// <summary>
        /// Разблокирует точечный рисунок из системной памяти.
        /// </summary>
        /// <param name="b">Точечный рисунок</param>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        public static void Unlock(this Bitmap b, BitmapData bmData)
        {
            b.UnlockBits(bmData);
            return;
        }
        #endregion

        #region RGB
        /// <summary>
        /// Преобразовывает точечный рисунок в RGBA-структуру.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>RGBA-структура</returns>
        public static double[][,] ToRGB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToRGB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в RGBA-структуру.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>RGBA-структура</returns>
        public unsafe static double[][,] ToRGB(this BitmapData bmData, bool alpha = false)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int total = 4, length = (alpha) ? total : total - 1;
            double[][,] array = BitmapConverter.Dim(length, height, width);
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * total;

                    // color
                    array[0][j, i] = p[k] / 255.0;
                    array[1][j, i] = p[k + 1] / 255.0;
                    array[2][j, i] = p[k + 2] / 255.0;

                    // if exist?
                    if (alpha)
                    {
                        array[3][j, i] = p[k + 3] / 255.0;
                    }
                }
            });
            return array;
        }
        /// <summary>
        /// Преобразовывает RGBA-структуру в цветное изображение.
        /// </summary>
        /// <param name="array">RGBA-структура</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap FromRGB(this double[][,] array)
        {
            int total = 4, length = array.Length;
            int width = array[0].GetLength(1), height = array[0].GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * total;

                    // recording model:
                    p[k] = Maths.Byte(array[0][j, i] * 255.0);
                    p[k + 1] = Maths.Byte(array[1][j, i] * 255.0);
                    p[k + 2] = Maths.Byte(array[2][j, i] * 255.0);
                    p[k + 3] = (byte)((length < 4) ? 255 : array[3][j, i] * 255.0);
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        #endregion

        #region HSB
        /// <summary>
        /// Преобразовывает точечный рисунок в HSB-структуру.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>HSB-структура</returns>
        public static double[][,] ToHSB(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSB(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в HSB-структуру.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>HSB-структура</returns>
        public unsafe static double[][,] ToHSB(this BitmapData bmData, bool alpha = false)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int total = 4, length = (alpha) ? total : total - 1;
            double[][,] array = BitmapConverter.Dim(length, height, width);
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; HSB mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift
                    k = jstride + i * total;

                    // color transform:
                    rgb = new RGB(p[k + 2], p[k + 1], p[k + 0]);
                    mod = HSB.FromRGB(rgb);

                    // recording model:
                    array[0][j, i] = mod.Hue;
                    array[1][j, i] = mod.Saturation;
                    array[2][j, i] = mod.Brightness;

                    // if exist?
                    if (alpha)
                    {
                        array[3][j, i] = p[k + 3] / 255.0;
                    }
                }
            });
            return array;
        }
        /// <summary>
        /// Преобразовывает HSB-структуру в цветное изображение.
        /// </summary>
        /// <param name="array">HSB-структура</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap FromHSB(this double[][,] array)
        {
            int total = 4, length = array.Length;
            int width = array[0].GetLength(1), height = array[0].GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; HSB mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * total;

                    // color transform:
                    mod = new HSB(array[0][j, i], array[1][j, i], array[2][j, i]);
                    rgb = mod.ToRGB;

                    // recording model:
                    p[k] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                    p[k + 3] = (byte)((length < 4) ? 255 : array[3][j, i] * 255.0);
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        #endregion

        #region HSL
        /// <summary>
        /// Преобразовывает точечный рисунок в HSL-структуру.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>HSL-структура</returns>
        public static double[][,] ToHSL(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToHSL(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в HSL-структуру.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>HSL-структура</returns>
        public unsafe static double[][,] ToHSL(this BitmapData bmData, bool alpha = false)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int total = 4, length = (alpha) ? total : total - 1;
            double[][,] array = BitmapConverter.Dim(length, height, width);
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; HSL mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift
                    k = jstride + i * total;

                    // color transform:
                    rgb = new RGB(p[k + 2], p[k + 1], p[k + 0]);
                    mod = HSL.FromRGB(rgb);

                    // recording model:
                    array[0][j, i] = mod.Hue;
                    array[1][j, i] = mod.Saturation;
                    array[2][j, i] = mod.Lightness;

                    // if exist?
                    if (alpha)
                    {
                        array[3][j, i] = p[k + 3] / 255.0;
                    }
                }
            });
            return array;
        }
        /// <summary>
        /// Преобразовывает HSL-структуру в цветное изображение.
        /// </summary>
        /// <param name="array">HSL-структура</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap FromHSL(this double[][,] array)
        {
            int total = 4, length = array.Length;
            int width = array[0].GetLength(1), height = array[0].GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; HSL mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * total;

                    // color transform:
                    mod = new HSL(array[0][j, i], array[1][j, i], array[2][j, i]);
                    rgb = mod.ToRGB;

                    // recording model:
                    p[k] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                    p[k + 3] = (byte)((length < 4) ? 255 : array[3][j, i] * 255.0);
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        #endregion

        #region YCbCr
        /// <summary>
        /// Преобразовывает точечный рисунок в YCbCr-структуру.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>YCbCr-структура</returns>
        public static double[][,] ToYCbCr(this Bitmap Data, bool alpha = false)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[][,] rgb = BitmapConverter.ToYCbCr(bmData, alpha);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в YCbCr-структуру.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="alpha">Учитывать альфа-канал или нет</param>
        /// <returns>YCbCr-структура</returns>
        public unsafe static double[][,] ToYCbCr(this BitmapData bmData, bool alpha = false)
        {
            int width = bmData.Width, height = bmData.Height, stride = bmData.Stride;
            int total = 4, length = (alpha) ? total : total - 1;
            double[][,] array = BitmapConverter.Dim(length, height, width);
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; YCbCr mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift
                    k = jstride + i * total;

                    // color transform:
                    rgb = new RGB(p[k + 2], p[k + 1], p[k + 0]);
                    mod = YCbCr.FromRGB(rgb);

                    // recording model:
                    array[0][j, i] = mod.Y;
                    array[1][j, i] = mod.Cb;
                    array[2][j, i] = mod.Cr;

                    // if exist?
                    if (alpha)
                    {
                        array[3][j, i] = p[k + 3] / 255.0;
                    }
                }
            });
            return array;
        }
        /// <summary>
        /// Преобразовывает YCbCr-структуру в цветное изображение.
        /// </summary>
        /// <param name="array">YCbCr-структура</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap FromYCbCr(this double[][,] array)
        {
            int total = 4, length = array.Length;
            int width = array[0].GetLength(1), height = array[0].GetLength(0);
            Bitmap bitmap = new Bitmap(width, height);
            BitmapData bmData = BitmapConverter.Lock32bpp(bitmap);
            int stride = bmData.Stride;
            byte* p = (byte*)bmData.Scan0.ToPointer();

            Parallel.For(0, height, j =>
            {
                RGB rgb; YCbCr mod;
                int i, k, jstride = j * stride;

                for (i = 0; i < width; i++)
                {
                    // shift:
                    k = jstride + i * total;

                    // color transform:
                    mod = new YCbCr(array[0][j, i], array[1][j, i], array[2][j, i]);
                    rgb = mod.ToRGB;

                    // recording model:
                    p[k] = rgb.Blue;
                    p[k + 1] = rgb.Green;
                    p[k + 2] = rgb.Red;
                    p[k + 3] = (byte)((length < 4) ? 255 : array[3][j, i] * 255.0);
                }
            });

            BitmapConverter.Unlock(bitmap, bmData);
            return bitmap;
        }
        #endregion

        #region Private voids
        /// <summary>
        /// Возвращает набор матриц в соответствии с заданными параметрами.
        /// </summary>
        /// <param name="length">Длина набора</param>
        /// <param name="height">Высота матриц</param>
        /// <param name="width">Ширина матриц</param>
        /// <returns>Набор матриц</returns>
        internal static double[][,] Dim(int length, int height, int width)
        {
            double[][,] matx = new double[length][,];

            for (int k = 0; k < length; k++)
            {
                matx[k] = new double[height, width];
            }
            return matx;
        }
        #endregion

        #region Bitmap matrix voids
        /// <summary>
        /// Преобразовывает точечный рисунок в матрицу значений усредненного канала.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Прямоугольная матрица</returns>
        public static double[,] FromBitmap(this Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            double[,] rgb = FromBitmap(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Преобразовывает точечный рисунок в матрицу значений усредненного канала.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Прямоугольная матрица</returns>
        public unsafe static double[,] FromBitmap(this BitmapData bmData)
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
        /// Преобразовывает прямоугольную матрицу значений каналов в монохромный точечный рисунок.
        /// </summary>
        /// <param name="m">Прямоугольная матрица</param>
        /// <returns>Точечный рисунок</returns>
        public unsafe static Bitmap ToBitmap(this double[,] m)
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
        #endregion
    }
    #endregion

    #region Mathematics
    /// <summary>
    /// Используется для работы с яркостью, представленной в виде значения, принадлежащего интервалу [0, 1].
    /// </summary>
    public static class Intensity
    {
        #region Nonlinear methods components
        /// <summary>
        /// Реализует алгоритм коррекции яркости Single Scale Retinex.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="xlow">Яркость под воздействием фильтра</param>
        /// <param name="nbase">Основание логарифма</param>
        /// <param name="a">Коэффициент сжатия [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SingleScaleRetinex(double x, double xlow, double nbase, double a, double b)
        {
            // Singe scale retinex modified algorithm
            // by Asiryan Valeriy
            // 
            return Math.Exp(a * Math.Log(x / xlow, nbase) + b) - 0.5f;
        }
        /// <summary>
        /// Возвращает маску коррекции яркости Single Scale Retinex.
        /// </summary>
        /// <param name="nbase">Основание логарифма</param>
        /// <param name="a">Коэффициент сжатия (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Реализует алгоритм инверсии локального контраста.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="xlow">Яркость под воздействием НЧ-фильтра</param>
        /// <param name="a">Коэффициент сжатия (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double LocalContrastInversion(double x, double xlow, double a, double b)
        {
            return a * x / (xlow + b);
        }
        /// <summary>
        /// Возвращает маску инверсии локального контраста.
        /// </summary>
        /// <param name="a">Коэффициент сжатия (0, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Реализует алгоритм локального улучшения контраста.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="xlow">Яркость под воздействием НЧ-фильтра</param>
        /// <param name="a">Коэффициент сжатия [-1, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double LocalContrastEnhancement(double x, double xlow, double a)
        {
            return x + a * (x - xlow);

            //return (1.0 + a) * (x - xlow) + xlow;
        }
        /// <summary>
        /// Возвращает маску локального улучшения контраста.
        /// </summary>
        /// <param name="a">Коэффициент сжатия [-1, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Реализует алгоритм гомоморфной обработки.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="mu">Яркость под воздействием НЧ-фильтра</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double HomomorphicEnhancement(double x, double mu, double a, double b)
        {
            return Math.Exp(Math.Log(x) - a * Math.Log(mu + b));
        }
        /// <summary>
        /// Возвращает маску гомоморфной обработки.
        /// </summary>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение (0, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Реализует алгоритм ξ-функции повышения контраста.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="mu">Математическое ожидание</param>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение [-1, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double KsiContrastEnchancement(double x, double mu, double a, double b)
        {
            // Пусть:
            // x ∈ [0, 1], μ ∈ [0, 1] - математическое ожидание величины x.
            // σ - дисперсия, ξ - кси-коэффициент.
            //
            // σ = x - μ
            // ξ = σ / μ
            // Тогда итоговое значение яркости будет вычисляться по формуле:
            // x' = x + α * ξ + β, где α, β ∈ [-1, 1].

            double sigma = x - mu;
            double ksi = sigma / mu;
            return x + a * ksi + b;
        }
        /// <summary>
        /// Возвращает маску ξ-функции повышения контраста.
        /// </summary>
        /// <param name="a">Коэффициент сжатия контраста [-1, 1]</param>
        /// <param name="b">Смещение [-1, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Реализует SAUCE-алгоритм.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="mu">Математическое ожидание</param>
        /// <param name="d">Степень отличия [0, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Возвращает маску SAUCE-алгоритма.
        /// </summary>
        /// <param name="a">Коэффициент сжатия [-1, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        #endregion

        #region LogStretrch Histogram method components
        /// <summary>
        /// Натуральный логарифм от 0.5.
        /// </summary>
        public const double log05 = -0.693147180559945;
        /// <summary>
        /// Погрешность вычисления логарифма.
        /// </summary>
        public const double logEpsilon = 1e-9;
        /// <summary>
        /// Реализует логарифимическое сжатие случайной величины.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="mu">Математическое ожидание</param>
        /// <param name="s">Уровень теней</param>
        /// <param name="l">Уровень светов</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double LogStretch(double x, double mu, double s, double l)
        {
            return Intensity.LogPow(x, Maths.Range(Intensity.log05 / Math.Log(mu), s, l));
        }
        /// <summary>
        /// Возвращает маску логарифимического сжатия.
        /// </summary>
        /// <param name="s">Уровень теней</param>
        /// <param name="l">Уровень светов</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
        /// Возвращает число, возведенное в логарифмическую степень.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="power">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double LogPow(double a, double power)
        {
            return Math.Exp(Math.Log(a) * power);
        }
        #endregion

        #region Gamma method components
        /// <summary>
        /// Получает маску для гамма-коррекции.
        /// </summary>
        /// <param name="g">Гамма</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует гамму-коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="g">Гамма</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Gamma(double x, double g)
        {
            return Math.Pow(x, g);
        }
        #endregion

        #region Shift method components
        /// <summary>
        /// Получает маску для коррекции ко смещением случайной величины.
        /// </summary>
        /// <param name="b">Сдвиг (-0.5, 0.5)</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует коррекцию ко смещением случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="b">Сдвиг (-0.5, 0.5)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Shift(double x, double b)
        {
            double v = log05 / Math.Log(0.5 - b);
            return LogPow(x, v);
        }
        #endregion

        #region Binarization method components
        /// <summary>
        /// Получает маску для бинаризации.
        /// </summary>
        /// <param name="threshold">Пороговое значение [0, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует бинаризацию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="threshold">Пороговое значение [0, 1]</param>
        /// <returns>Целое число без знака</returns>
        public static double Bin(double x, double threshold)
        {
            return (x > threshold) ? 1 : 0;
        }
        #endregion

        #region Exposure method components
        /// <summary>
        /// Получает маску для экспозиционной коррекции.
        /// </summary>
        /// <param name="average">Среднее число</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует экспозиционную коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="average">Среднее число</param>
        /// <returns>Целое число без знака</returns>
        public static double Exposure(double x, double average)
        {
            double T = 255.0 / average;
            return 1 - Math.Exp(-T * x);
        }
        #endregion

        #region Sin method components
        /// <summary>
        /// Получает маску для синусоидальной коррекции.
        /// </summary>
        /// <param name="delta">Дельта</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует синусоидальную коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="delta">Дельта</param>
        /// <returns>Целое число без знака</returns>
        public static double Sin(double x, double delta)
        {
            return 0.5 * Math.Sin((3.14 * x) - (3.14 / 2)) + 0.5 + delta;
        }
        #endregion

        #region Cos method components
        /// <summary>
        /// Получает маску для косинусоидальной коррекции.
        /// </summary>
        /// <param name="delta">Дельта</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует косинусоидальную коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="delta">Дельта</param>
        /// <returns>Целое число без знака</returns>
        public static double Cos(double x, double delta)
        {
            return 0.5 * Math.Cos((3.14 * x) - 3.14) + 0.5 + delta;
        }
        #endregion

        #region Log method components
        /// <summary>
        /// Получает маску для логарифмической коррекции.
        /// </summary>
        /// <param name="a">Основание логарифма</param>
        /// <param name="delta">Дельта</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует логарифмическую коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="a">Основание логарифма</param>
        /// <param name="delta">Дельта</param>
        /// <returns>Целое число без знака</returns>
        public static double Log(double x, double a, double delta)
        {
            return Math.Log((1.0 + (x + delta) / 0.5), a);
        }
        #endregion

        #region Add method components
        /// <summary>
        /// Получает маску для коррекции по значению Y = (X + V).
        /// </summary>
        /// <param name="value">Значение</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Add(double value, int length)
        {
            double[] table = new double[length];
            for (int x = 0; x < length; x++)
            {
                table[x] = x / (double)length + value;
            }
            return table;
        }
        #endregion

        #region Contrast method components
        /// <summary>
        /// Получает маску для коррекции контрастности.
        /// </summary>
        /// <param name="value">Значение</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует коррекцию контрастности случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="value">Контрастность</param>
        /// <returns>Целое число без знака</returns>
        public static double Contrast(double x, double value)
        {
            value = (1 + value); // Вычисление общего значения.
            double xc = x;
            xc -= 0.5f;
            xc *= value;
            xc += 0.5f;
            return xc;
        }
        #endregion

        #region Log contrast method components
        /// <summary>
        /// Получает маску для коррекции контрастности.
        /// </summary>
        /// <param name="power">Значение</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует алгоритм глобальной коррекции контраста.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="power">Степень</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double LogContrast(double x, double power)
        {
            if (x <= 0.5)
            {
                return Intensity.LogPow(x * 2, power) * 0.5;
            }
            return 1.0 - Intensity.LogPow((1 - x) * 2, power) * 0.5;
        }
        #endregion

        #region Invert method components
        /// <summary>
        /// Получает маску инверсии.
        /// </summary>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Инвертирует случайную величину.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Целое число без знака</returns>
        public static double Invert(double x)
        {
            return 1.0 - x;
        }
        #endregion

        #region Equalization components
        /// <summary>
        /// Эквализирует случайную величину относительно {min, max} диапазона.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="max">Максимальное значение, которое принимает случайная величина</param>
        /// <param name="min">Минимальное значение, которое принимает случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Equalize(double x, double min, double max)
        {
            double a = max - min;
            double b = x - min;
            double c = (a != 0) ? b / a : x;
            return Maths.Double(c);
        }
        /// <summary>
        /// Получает маску для эквализации случайной величины относительно {min, max} диапазона.
        /// </summary>
        /// <param name="max">Максимальное значение, которое принимает случайная величина</param>
        /// <param name="min">Минимальное значение, которое принимает случайная величина</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Equalize(double min, double max, int length)
        {
            double[] table = new double[length];

            for (int x = 0; x < length; x++)
            {
                table[x] = Intensity.Equalize(x / (double)length, min, max);
            }
            return table;
        }
        #endregion

        #region Linear method components
        /// <summary>
        /// Получает маску для линейной коррекции.
        /// </summary>
        /// <param name="range">Пара чисел Max и Min</param>
        /// <param name="delta">Дельта</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Linear(RangeDouble range, double delta, int length)
        {
            return Intensity.Linear(range.Max, range.Min, delta, length);
        }
        /// <summary>
        /// Получает маску для линейной коррекции.
        /// </summary>
        /// <param name="xmax">Максимальное значение, принимаемое случайной величиной</param>
        /// <param name="xmin">Минимальное значение, принимаемое случайной величиной</param>
        /// <param name="delta">Дельта</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Реализует линейную коррекцию случайной величины.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <param name="xmax">Максимальное значение, принимаемое случайной величиной</param>
        /// <param name="xmin">Минимальное значение, принимаемое случайной величиной</param>
        /// <param name="delta">Дельта</param>
        /// <returns>Целое число без знака</returns>
        public static double Linear(double x, double xmax, double xmin, double delta)
        {
            return (x - xmin) / (xmax - xmin) + delta;
        }
        #endregion

        #region Levels method components
        /// <summary>
        /// Получает маску для коррекции уровней.
        /// </summary>
        /// <param name="xmin">Минимальное значение входного параметра</param>
        /// <param name="xmax">Максимальное значение входного параметра</param>
        /// <param name="ymin">Минимальное значение выходного параметра</param>
        /// <param name="ymax">Максимальное значение выходного параметра</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
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
        /// Получает маску для коррекции уровней.
        /// </summary>
        /// <param name="input">Входные значения</param>
        /// <param name="output">Выходные значения</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Levels(RangeDouble input, RangeDouble output, int length)
        {
            return Intensity.Levels(input.Min, input.Max, output.Min, output.Max, length);
        }
        #endregion

        #region Quantization method components
        /// <summary>
        /// Получает маску для квантования по уровню.
        /// </summary>
        /// <param name="levels">Количество уровней квантования</param>
        /// <param name="length">Количество уровней представления</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Quantize(int levels, int length)
        {
            if (levels > length)
                throw new Exception("Количество уровней квантования не может быть больше количества уровней представления");

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

        #region Bradley method components
        /// <summary>
        /// Реализует алгоритм коррекции яркости Single Scale Retinex.
        /// </summary>
        /// <param name="x">Яркость</param>
        /// <param name="xlow">Яркость под воздействием фильтра</param>
        /// <param name="difference">Предел разницы яркости между пикселем обработки и средним значением локальных пикселей [0, 1]</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Bradley(double x, double xlow, double difference = 0.15)
        {
            // Bradley local threshold void.
            // Derek Bradley, Gerhard Roth (2005). Adaptive Thresholding Using the Integral Image.
            // Retrieved from http://www.scs.carleton.ca/~roth/iit-publications-iti/docs/gerh-50002.pdf

            double z = 1.0 - difference;
            return (x < xlow * z) ? 0 : 1;
        }
        /// <summary>
        /// Возвращает маску коррекции яркости Single Scale Retinex.
        /// </summary>
        /// <param name="difference">Предел разницы яркости между пикселем обработки и средним значением локальных пикселей [0, 1]</param>
        /// <param name="length">Размерность массива</param>
        /// <returns>Матрица</returns>
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
    }
    /// <summary>
    /// Используется для смешивания слоев.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// http://www.pegtop.net/delphi/articles/blendmodes/index.htm
    /// </remarks>
    /// </summary>
    public static class BlendMode
    {
        #region Blend mode components
        /// <summary>
        /// Реализует функцию усреднения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Average(double a, double b)
        {
            return (a + b) / 2.0;
        }
        /// <summary>
        /// Реализует функцию экранирования.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Screen(double a, double b)
        {
            return 1 - (1 - a) * (1 - b);
        }
        /// <summary>
        /// Реализует функцию вычитания.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Difference(double a, double b)
        {
            return Math.Abs(a - b);
        }
        /// <summary>
        /// Реализует функцию обращения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Negation(double a, double b)
        {
            return 1 - Math.Abs(1 - a - b);
        }
        /// <summary>
        /// Реализует функцию исключения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Exclusion(double a, double b)
        {
            return a + b - 2 * a * b;
        }
        /// <summary>
        /// Реализует функцию наложения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Overlay(double a, double b)
        {
            if (a < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Реализует функцию "жесткий свет".
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double HardLight(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b;
            }
            return 1 - 2 * (1 - a) * (1 - b);
        }
        /// <summary>
        /// Реализует функцию поворота цвета.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Dodge(double a, double b)
        {
            return a / (1 - b);
        }
        /// <summary>
        /// Реализует функцию "умного" поворота цвета.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SoftDodge(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * a / (1 - b);
            }
            return 1 - 0.5 * (1 - b) / a;
        }
        /// <summary>
        /// Реализует функцию обратного экранирования.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Burn(double a, double b)
        {
            return 1 - (1 - a) / b;
        }
        /// <summary>
        /// Реализует функцию "умного" обратного экранирования.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double SoftBurn(double a, double b)
        {
            if (a + b < 1)
            {
                return 0.5 * b / (1 - a);
            }
            return 1 - 0.5 * (1 - a) / b;
        }
        /// <summary>
        /// Реализует функцию отражения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Reflect(double a, double b)
        {
            return a * a / (1 - b);
        }
        /// <summary>
        /// Реализует функцию свечения.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Glow(double a, double b)
        {
            return b * b / (1 - a);
        }
        /// <summary>
        /// Реализует функцию печати.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Stamp(double a, double b)
        {
            return a + 2 * b - 1;
        }
        /// <summary>
        /// ВРеализует функцию "freeze".
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Freeze(double a, double b)
        {
            double x = 1 - a;
            return 1 - x * x / b;
        }
        /// <summary>
        /// Реализует функцию "heat".
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Heat(double a, double b)
        {
            double x = 1 - b;
            return 1 - x * x / a;
        }
        /// <summary>
        /// Реализует функцию интерполяции.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Interpolation(double a, double b)
        {
            return 0.5 - 0.25 * Math.Cos(Math.PI * a) - 0.25 * Math.Cos(Math.PI * b);
        }
        /// <summary>
        /// Реализует функцию "умный свет" Adobe Photoshop.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Photoshop(double a, double b)
        {
            if (b < 0.5)
            {
                return 2 * a * b + a * a * (1 - 2 * a);
            }
            return 2 * a * (1 - b) + Math.Sqrt(a) * (2 * b - 1);
        }
        /// <summary>
        /// Реализует функцию "умный свет" Illusions.hu.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Illusions(double a, double b)
        {
            double x = 2 * (0.5 - b);
            double y = Math.Pow(2, x);
            return Math.Pow(a, y);
        }
        /// <summary>
        /// Реализует функцию "умный свет" Pegtop.
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Pegtop(double a, double b)
        {
            return (1 - 2 * b) * a * a + 2 * b * a;
        }
        /// <summary>
        /// Реализует функцию "Cairo".
        /// </summary>
        /// <param name="a">Первый слой</param>
        /// <param name="b">Второй слой</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Fw3c(double a, double b)
        {
            if (b <= 0.5)
            {
                return a - (1 - 2 * b) * a * (1 - a);
            }
            return a + (2 * b - 1) * (BlendMode.Gw3c(a) - a);
        }
        /// <summary>
        /// Реализует функцию "Cairo".
        /// </summary>
        /// <param name="a">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
    /// Используется для работы со статистическими характеристиками изображения.
    /// </summary>
    public static class Statistics
    {
        #region Histogram
        /// <summary>
        /// Получает гистограмму изображения.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Одномерный массив</returns>
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
        /// Получает гистограмму изображения.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        /// <returns>Одномерный массив</returns>
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
        /// Получает данные одновременно о красном, зеленом и синем каналах.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Одномерный массив</returns>
        public static int[] Histogram(Bitmap Data)
        {
            BitmapData bmData = BitmapConverter.Lock32bpp(Data);
            int[] rgb = Histogram(bmData);
            BitmapConverter.Unlock(Data, bmData);
            return rgb;
        }
        /// <summary>
        /// Получает данные одновременно о красном, зеленом и синем каналах.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <param name="channel">Цветовой канал модели RGBA</param>
        /// <returns>Одномерный массив</returns>
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
        /// Получает массив значений общей функции плотности яркости.
        /// </summary>
        /// <param name="H">Гистограмма</param>
        /// <returns>Одномерный массив</returns>
        public static int[] CDF(int[] H)
        {
            int length = H.Length;
            int[] cdf = new int[length];
            cdf[0] = H[0];

            // Рекурсивный метод расчета функции плотности яркости
            // взамен реккурентному:
            for (int i = 1; i < length; i++)
            {
                cdf[i] = H[i] + cdf[i - 1];
            }
            return cdf;
        }
        /// <summary>
        /// Получает массив значений эквализированной гистограммы путем пересчета функции плотности яркости.
        /// </summary>
        /// <param name="H">Гистограмма</param>
        /// <returns>Одномерный массив</returns>
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
        /// Вычисляет оптимальный порог методом Оцу для исходного точечного рисунка.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Целое число со знаком</returns>
        public static int OtsuThreshold(BitmapData bmData)
        {
            double[] v = new double[256];
            int[] h = Statistics.Histogram(bmData); // Получаем гистограмму исходного точечного рисунка.
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
        /// Вычисляет оптимальный порог методом Оцу для исходного точечного рисунка.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Вычисляет оптимальный порог для исходного точечного рисунка.
        /// </summary>
        /// <param name="bmData">Атрибуты точечного изображения</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Вычисляет оптимальный порог для исходного точечного рисунка.
        /// </summary>
        /// <param name="Data">Точечный рисунок</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Получает индекс максимального элемента массива.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Получает индекс максимального элемента массива.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Получает индекс минимального элемента массива.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Получает индекс минимального элемента массива.
        /// </summary>
        /// <param name="data">Одномерный массив</param>
        /// <returns>Целое число со знаком</returns>
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
        /// Получает значение вероятности класса.
        /// </summary>
        /// <param name="init">Начальный момент</param>
        /// <param name="end">Конечный момент</param>
        /// <param name="h">Гистограмма точечного рисунка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Получает значение среднего арифметического класса.
        /// </summary>
        /// <param name="init">Начальный момент</param>
        /// <param name="end">Конечный момент</param>
        /// <param name="h">Гистограмма точечного рисунка</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
    /// Используется для работы с матрицами точек.
    /// </summary>
    public static class PointMatrix
    {
        #region Filters
        /// <summary>
        /// Получает массив упорядоченных пар чисел X и Y для эффекта отображения по оси Y.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта отображения по оси X.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта сдвига по оси X.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Сдвиг</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта сдвига по оси Y.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Сдвиг</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта шума.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Плотность [0, 100]</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Плотность [0, 100]</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта мозаики.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Плотность [0, 100]</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
        /// Получает массив упорядоченных пар чисел X и Y для эффекта воды.
        /// </summary>
        /// <param name="width">Ширина изображения</param>
        /// <param name="height">Высота изображения</param>
        /// <param name="value">Плотность [0, 100]</param>
        /// <returns>Массив упорядоченных пар чисел X и Y</returns>
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
    /// Определяет шум Перлина.
    /// <remarks>
    /// Perlin noise (Шум Перлина, также иногда Классический шум Перлина) — математический алгоритм по генерированию процедурной текстуры 
    /// псевдо-случайным методом. Используется в компьютерной графике для увеличения реализма или графической сложности поверхности геометрических объектов. 
    /// Также может использоваться для генерации эффектов дыма, тумана и т.д.
    /// 
    /// Более подробную информацию можно найти на сайте:
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
        /// Инициализирует шум Перлина.
        /// </summary>
        /// <param name="octaves">Количество октав [1, 32]</param>
        /// <param name="persistence">Продолжительность</param>
        /// <param name="frequency">Частота</param>
        /// <param name="amplitude">Амплитуда</param>
        public PerlinNoise(int octaves = 4, double persistence = 0.65, double frequency = 1, double amplitude = 1)
        {
            Octaves = octaves; Persistence = persistence; Frequency = frequency; Amplitude = amplitude;
        }
        /// <summary>
        /// Получает или задает частоту.
        /// </summary>
        public double Frequency
        {
            get { return initFrequency; }
            set { initFrequency = value; }
        }
        /// <summary>
        /// Получает или задает амплитуду.
        /// </summary>
        public double Amplitude
        {
            get { return initAmplitude; }
            set { initAmplitude = value; }
        }
        /// <summary>
        /// Получает или задает продолжительность.
        /// </summary>
        public double Persistence
        {
            get { return persistence; }
            set { persistence = value; }
        }
        /// <summary>
        /// Получает или задает количество октав [1, 32].
        /// </summary>
        public int Octaves
        {
            get { return octaves; }
            set { octaves = System.Math.Max(1, System.Math.Min(32, value)); }
        }
        /// <summary>
        /// Одномерная функция шума Перлина.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Двумерная функция шума Перлина.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="y">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// Одномерная функция шума.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private static double Noise(int x)
        {
            int n = (x << 13) ^ x;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// Двумерная функция шума.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="y">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private static double Noise(int x, int y)
        {
            int n = x + y * 57;
            n = (n << 13) ^ n;

            return (1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
        }
        /// <summary>
        /// Одномерная функция сглаженного шума.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private static double SmoothedNoise(double x)
        {
            int xInt = (int)x;
            double xFrac = x - xInt;

            return PerlinNoise.CosineInterpolate(Noise(xInt), Noise(xInt + 1), xFrac);
        }
        /// <summary>
        /// Двумерная функция сглаженного шума.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="y">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
        /// <summary>
        /// Косинусная интерполяция.
        /// </summary>
        /// <param name="x1">Неизвестно</param>
        /// <param name="x2">Неизвестно</param>
        /// <param name="a">Параметр</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
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
    /// Определеяет цветовое пространство.
    /// </summary>
    public enum Space
    {
        /// <summary>
        /// Цветовое пространство RGB.
        /// </summary>
        RGB,
        /// <summary>
        /// Цветовое пространство HSB.
        /// </summary>
        HSB,
        /// <summary>
        /// Цветовое пространство HSB.
        /// </summary>
        HSL,
        /// <summary>
        /// Цветовое пространство YCbCr.
        /// </summary>
        YCbCr,
        /// <summary>
        /// Оттенки серого.
        /// </summary>
        Grayscale,
    }
    /// <summary>
    /// Определяет цветовой канал модели RGBA.
    /// </summary>
    public enum RGBA
    {
        /// <summary>
        /// Альфа-канал.
        /// </summary>
        Alpha = 3,
        /// <summary>
        /// Красный.
        /// </summary>
        Red = 2,
        /// <summary>
        /// Зеленый.
        /// </summary>
        Green = 1,
        /// <summary>
        /// Синий.
        /// </summary>
        Blue = 0,
    }
    /// <summary>
    /// Указывает направление вектора градиента.
    /// </summary>
    public enum Gradient
    {
        /// <summary>
        /// Северное направление вектора градиента.
        /// </summary>
        North = 0,
        /// <summary>
        /// Северо-Западное направление вектора градиента.
        /// </summary>
        NorthWest = 1,
        /// <summary>
        /// Западное направление вектора градиента.
        /// </summary>
        West = 2,
        /// <summary>
        /// Юго-Западное направление вектора градиента.
        /// </summary>
        SouthWest = 3,
        /// <summary>
        /// Южное направление вектора градиента.
        /// </summary>
        South = 4,
        /// <summary>
        /// Юго-Восточное направление вектора градиента.
        /// </summary>
        SouthEast = 5,
        /// <summary>
        /// Восточное направление вектора градиента.
        /// </summary>
        East = 6,
        /// <summary>
        /// Северо-Восточное направление вектора градиента.
        /// </summary>
        NorthEast = 7,
    }
    #endregion
}
