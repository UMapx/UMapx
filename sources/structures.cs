// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Runtime.Serialization;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Structures: Range, Point, Size
    /// <summary>
    /// Определяет пару дробных чисел, представляющих отрезок.
    /// </summary>
    public struct RangeDouble : ICloneable, ISerializable
    {
        #region Private data
        private double max;
        private double min;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих отрезок.
        /// </summary>
        /// <param name="min">Нижняя граница отрезка</param>
        /// <param name="max">Верхняя граница отрезка</param>
        public RangeDouble(double min, double max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Получает или задает нижнюю границу отрезка.
        /// </summary>
        public double Min
        {
            get
            {
                return this.min;
            }
            set
            {
                this.min = value;
            }
        }
        /// <summary>
        /// Получает или задает верхнюю границу отрезка.
        /// </summary>
        public double Max
        {
            get
            {
                return this.max;
            }
            set
            {
                this.max = value;
            }
        }
        /// <summary>
        /// Проверяет находится ли случайная величина в заданном интервале.
        /// </summary>
        /// <param name="x">Число</param>
        /// <returns>Логическое значение</returns>
        public bool IsOnRange(double x)
        {
            if ((x >= this.min) && (x <= this.max))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return min.GetHashCode() ^ max.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа RangeDouble.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is RangeDouble) ? (this == (RangeDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа RangeDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(RangeDouble a, RangeDouble b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа RangeDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(RangeDouble a, RangeDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new RangeDouble(min, max);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public RangeDouble Clone()
        {
            return new RangeDouble(min, max);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Min", this.Min);
            info.AddValue("Max", this.Max);
        }
        #endregion
    }
    /// <summary>
    /// Определяет пару целых чисел, представляющих отрезок.
    /// </summary>
    public struct RangeInt : ICloneable, ISerializable
    {
        #region Private data
        private int max;
        private int min;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих отрезок.
        /// </summary>
        /// <param name="min">Нижняя граница отрезка</param>
        /// <param name="max">Верхняя граница отрезка</param>
        public RangeInt(int min, int max)
        {
            this.min = min;
            this.max = max;
        }
        /// <summary>
        /// Получает или задает нижнюю границу отрезка.
        /// </summary>
        public int Min
        {
            get
            {
                return this.min;
            }
            set
            {
                this.min = value;
            }
        }
        /// <summary>
        /// Получает или задает верхнюю границу отрезка.
        /// </summary>
        public int Max
        {
            get
            {
                return this.max;
            }
            set
            {
                this.max = value;
            }
        }
        /// <summary>
        /// Проверяет находится ли случайная величина в заданном интервале.
        /// </summary>
        /// <param name="x">Число</param>
        /// <returns>Логическое значение</returns>
        public bool IsOnRange(int x)
        {
            if ((x >= this.min) && (x <= this.max))
            {
                return true;
            }
            return false;
        }
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return min.GetHashCode() ^ max.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", min, max);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа RangeInt.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is RangeInt) ? (this == (RangeInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа RangeInt между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(RangeInt a, RangeInt b)
        {
            return (a.Max == b.Max && a.Min == b.Min);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа RangeInt между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(RangeInt a, RangeInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new RangeInt(min, max);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public RangeInt Clone()
        {
            return new RangeInt(min, max);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Min", this.Min);
            info.AddValue("Max", this.Max);
        }
        #endregion
    }
    /// <summary>
    /// Определяет пару дробных чисел, представляющих упорядоченную пару координат X и Y.
    /// </summary>
    public struct PointDouble : ICloneable, ISerializable
    {
        #region Private data
        private double y;
        private double x;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих упорядоченную пару координат X и Y.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        public PointDouble(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
        /// <summary>
        /// Получает или задает координату X.
        /// </summary>
        public double X
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
        /// Получает или задает координату Y.
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
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return x.GetHashCode() ^ y.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", x, y);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа PointDouble.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is PointDouble) ? (this == (PointDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа PointDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(PointDouble a, PointDouble b)
        {
            return (a.X == b.X && a.Y == b.Y);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа PointDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(PointDouble a, PointDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new PointDouble(x, y);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public PointDouble Clone()
        {
            return new PointDouble(x, y);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("X", this.x);
            info.AddValue("Y", this.y);
        }
        #endregion
    }
    /// <summary>
    /// Определяет пару целых чисел, представляющих упорядоченную пару координат X и Y.
    /// </summary>
    public struct PointInt : ICloneable, ISerializable
    {
        #region Private data
        private int y;
        private int x;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих упорядоченную пару координат X и Y.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        public PointInt(int x, int y)
        {
            this.x = x;
            this.y = y;
        }
        /// <summary>
        /// Получает или задает координату X.
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
        /// Получает или задает координату Y.
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
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return x.GetHashCode() ^ y.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", x, y);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа PointInt.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is PointInt) ? (this == (PointInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа PointInt между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(PointInt a, PointInt b)
        {
            return (a.X == b.X && a.Y == b.Y);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа PointInt между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(PointInt a, PointInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new PointInt(x, y);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public PointInt Clone()
        {
            return new PointInt(x, y);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("X", this.x);
            info.AddValue("Y", this.y);
        }
        #endregion
    }
    /// <summary>
    /// Определяет пару дробных чисел, представляющих упорядоченную пару ширины и высоты.
    /// </summary>
    public struct SizeDouble : ICloneable, ISerializable
    {
        #region Private data
        private double width;
        private double height;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих упорядоченную пару ширины и высоты.
        /// </summary>
        /// <param name="width">Ширина</param>
        /// <param name="height">Высота</param>
        public SizeDouble(double width, double height)
        {
            this.height = height;
            this.width = width;
        }
        /// <summary>
        /// Получает или задает нижнюю границу отрезка.
        /// </summary>
        public double Height
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
        /// Получает или задает верхнюю границу отрезка.
        /// </summary>
        public double Width
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
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return width.GetHashCode() ^ height.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", width, height);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа SizeDouble.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is SizeDouble) ? (this == (SizeDouble)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа SizeDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(SizeDouble a, SizeDouble b)
        {
            return (a.Width == b.Width && a.Height == b.Height);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа SizeDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(SizeDouble a, SizeDouble b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new SizeDouble(width, height);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public SizeDouble Clone()
        {
            return new SizeDouble(width, height);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Width", this.width);
            info.AddValue("Height", this.height);
        }
        #endregion
    }
    /// <summary>
    /// Определяет пару целых чисел, представляющих упорядоченную пару ширины и высоты.
    /// </summary>
    public struct SizeInt : ICloneable, ISerializable
    {
        #region Private data
        private int width;
        private int height;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует пару целых чисел, представляющих упорядоченную пару ширины и высоты.
        /// </summary>
        /// <param name="width">Ширина</param>
        /// <param name="height">Высота</param>
        public SizeInt(int width, int height)
        {
            this.height = height;
            this.width = width;
        }
        /// <summary>
        /// Получает или задает нижнюю границу отрезка.
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
        /// Получает или задает верхнюю границу отрезка.
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
        #endregion

        #region Overrides
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return width.GetHashCode() ^ height.GetHashCode();
        }
        /// <summary>
        /// Преобразует логическое значение в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return string.Format("({0}, {1})", width, height);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземпляр заданному значению типа SizeInt.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is SizeInt) ? (this == (SizeInt)obj) : false;
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа SizeInt между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(SizeInt a, SizeInt b)
        {
            return (a.Width == b.Width && a.Height == b.Height);
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа SizeDouble между собой.
        /// </summary>
        /// <param name="a">Пара чисел</param>
        /// <param name="b">Пара чисел</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(SizeInt a, SizeInt b)
        {
            return !(a == b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        object ICloneable.Clone()
        {
            return new SizeInt(width, height);
        }
        /// <summary>
        /// Создает копию пары чисел.
        /// </summary>
        /// <returns>Пара чисел</returns>
        public SizeInt Clone()
        {
            return new SizeInt(width, height);
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context)
        {
            info.AddValue("Width", this.width);
            info.AddValue("Height", this.height);
        }
        #endregion
    }
    #endregion
}
