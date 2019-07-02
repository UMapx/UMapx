// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Runtime.Serialization;
using System.Text.RegularExpressions;

namespace UMapx.Core
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                 UMAPX.CORE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Complex
    /// <summary>
    /// Определяет комплексное число.
    /// </summary>
    public struct Complex : ICloneable, ISerializable
    {
        #region Private data
        /// <summary>
        /// Действительная часть комплексного числа.
        /// </summary>
        public double Re;
        /// <summary>
        /// Мнимая часть комплексного числа.
        /// </summary>
        public double Im;
        #endregion

        #region Structure components
        /// <summary>
        /// Инициализирует комплексное число.
        /// </summary>
        /// <param name="re">Действительная часть комплексного числа</param>
        /// <param name="im">Мнимая часть комплексного числа</param>
        public Complex(double re, double im)
        {
            this.Re = re;
            this.Im = im;
        }
        /// <summary>
        /// Получает значение модуля.
        /// </summary>
        public double Abs
        {
            get
            {
                return Math.Sqrt(Re * Re + Im * Im);
            }
        }
        /// <summary>
        /// Получает значение угла фазы (аргумента).
        /// </summary>
        public double Angle
        {
            get
            {
                return Math.Atan2(Im, Re);
            }
        }
        /// <summary>
        /// Возвращает мнимую единицу.
        /// </summary>
        public static Complex I
        {
            get
            {
                return new Complex(0, 1);
            }
        }
        /// <summary>
        /// Возвращает действительную единицу.
        /// </summary>
        public static Complex One
        {
            get
            {
                return new Complex(1, 0);
            }
        }
        /// <summary>
        /// Возвращает комплексный ноль.
        /// </summary>
        public static Complex Zero
        {
            get
            {
                return new Complex(0, 0);
            }
        }
        /// <summary>
        /// Возвращает комплексно-сопряженное число.
        /// </summary>
        public Complex Conjugate
        {
            get
            {
                return new Complex(this.Re, -this.Im);
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
            return this.Re.GetHashCode() ^ this.Im.GetHashCode();
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземляр заданному значению типа Complex.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is Complex) ? (this == (Complex)obj) : false;
        }
        /// <summary>
        /// Преобразует комплексное число в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return this.ToString("G3");
        }
        /// <summary>
        /// Преобразует комплексное число в соответствующее ему строковое представление.
        /// </summary>
        /// <param name="format">Строка числового формата</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public string ToString(string format)
        {
            return StringOptions.Disp(new double[] { this.Re, this.Im }, format, StringOptions.C);
        }
        #endregion

        #region Bools
        /// <summary>
        /// Проверяет равны ли два объекта типа Comlex между собой.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(Complex a, Complex b)
        {
            return ((a.Re == b.Re) && (a.Im == b.Im));
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа Complex между собой.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(Complex a, Complex b)
        {
            return !(a == b);
        }
        #endregion

        #region Operators
        /// <summary>
        /// Сумма двух комплексных чисел
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator +(Complex a, Complex b)
        {
            return new Complex(a.Re + b.Re, a.Im + b.Im);
        }
        /// <summary>
        /// Сумма комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator +(Complex a, double b)
        {
            return new Complex(a.Re + b, a.Im);
        }
        /// <summary>
        /// Сумма комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator +(double a, Complex b)
        {
            return new Complex(b.Re + a, b.Im);
        }


        /// <summary>
        /// Разность двух комплексных чисел
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator -(Complex a, Complex b)
        {
            return new Complex(a.Re - b.Re, a.Im - b.Im);
        }
        /// <summary>
        /// Разность комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator -(Complex a, double b)
        {
            return new Complex(a.Re - b, a.Im);
        }
        /// <summary>
        /// Разность комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator -(double a, Complex b)
        {
            return new Complex(a - b.Re, b.Im);
        }
        /// <summary>
        /// Инвертирует комплексное число.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator -(Complex a)
        {
            return new Complex(-a.Re, -a.Im);
        }


        /// <summary>
        /// Произведение двух комплексных чисел
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator *(Complex a, Complex b)
        {
            double aRe = a.Re, aIm = a.Im;
            double bRe = b.Re, bIm = b.Im;

            return new Complex(aRe * bRe - aIm * bIm, aRe * bIm + aIm * bRe);
        }
        /// <summary>
        /// Произведение комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator *(double a, Complex b)
        {
            return new Complex(b.Re * a, b.Im * a);
        }
        /// <summary>
        /// Произведение комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator *(Complex a, double b)
        {
            return new Complex(a.Re * b, a.Im * b);
        }


        /// <summary>
        /// Частное двух комплексных чисел
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator /(Complex a, Complex b)
        {
            double aRe = a.Re, aIm = a.Im;
            double bRe = b.Re, bIm = b.Im;
            double abs = bRe * bRe + bIm * bIm;
            double inv = 1 / abs;

            return new Complex((aRe * bRe + aIm * bIm) * inv, (aIm * bRe - aRe * bIm) * inv);
        }
        /// <summary>
        /// Частное комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Комплексное число</param>
        /// <param name="b">Число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator /(Complex a, double b)
        {
            return new Complex(a.Re / b, a.Im / b);
        }
        /// <summary>
        /// Частное комплексного числа и действительного числа.
        /// </summary>
        /// <param name="a">Числа</param>
        /// <param name="b">Комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static Complex operator /(double a, Complex b)
        {
            if (b.Im == 0)
            {
                return new Complex(a / b.Re, 0);
            }
            else if (b.Re == 0)
            {
                return new Complex(0, a / b.Im);
            }
            return new Complex(a / b.Re, a / b.Im);
        }
        #endregion

        #region Conversion operators
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(double value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(float value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(long value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(ulong value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(short value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(ushort value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(int value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(uint value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(byte value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(sbyte value)
        {
            return new Complex(value, 0);
        }
        /// <summary>
        /// Определяет явное преобразование числа в комплексное число.
        /// </summary>
        /// <param name="value">Значение, преобразуемое в комплексное число</param>
        /// <returns>Комплексное число</returns>
        public static implicit operator Complex(decimal value)
        {
            return new Complex((double)value, 0);
        }
        #endregion

        #region Parsing
        /// <summary>
        /// Переводит исходную строку в комплексное число.
        /// <remarks>
        /// Примеры входной строки: "1 + 2i", "0.321 + 11i", ".1i".
        /// </remarks>
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public static Complex Parse(string s)
        {
            return StringOptions.Compar(s);
        }
        /// <summary>
        /// Пробует перевести исходную строку в комплексное число.
        /// </summary>
        /// <param name="complex">Исходная строка</param>
        /// <param name="result">Комплексное число</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string complex, ref Complex result)
        {
            try
            {
                result = Complex.Parse(complex);
                return true;
            }
            catch (FormatException)
            {
                result = new Complex();
                return false;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию комплексного числа.
        /// </summary>
        /// <returns>Комплексное число</returns>
        object ICloneable.Clone()
        {
            return new Complex(this.Re, this.Im);
        }
        /// <summary>
        /// Создает копию комплексного числа.
        /// </summary>
        /// <returns>Комплексное число</returns>
        public Complex Clone()
        {
            return new Complex(this.Re, this.Im);
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
            info.AddValue("Real", this.Re);
            info.AddValue("Imaginary", this.Im);
        }
        #endregion
    }
    #endregion

    #region Quaternion
    /// <summary>
    /// Определяет кватернион.
    /// <remarks>
    /// Кватернион - это система гиперкомплексных чисел, образующая векторное пространство размерностью четыре над полем вещественных чисел.
    /// </remarks>
    /// </summary>
    public struct Quaternion : ICloneable, ISerializable
    {
        #region Public data
        /// <summary>
        /// Значение X координаты вектора кватерниона.
        /// </summary>
        public double X;
        /// <summary>
        /// Значение Y координаты вектора кватерниона.
        /// </summary>
        public double Y;
        /// <summary>
        /// Значение Z координаты вектора кватерниона.
        /// </summary>
        public double Z;
        /// <summary>
        /// Получает координату поворота кватерниона.
        /// </summary>
        public double W;
        #endregion

        #region Quaternion components
        /// <summary>
        /// Создает кватернион на основе заданных координат.
        /// </summary>
        /// <param name="x">Координата X</param>
        /// <param name="y">Координата Y</param>
        /// <param name="z">Координата Z</param>
        /// <param name="w">Координата поворота</param>
        public Quaternion(double x, double y, double z, double w)
        {
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.W = w;
        }
        /// <summary>
        /// Получает кватернион, который представляет отсутствие вращения.
        /// </summary>
        public static Quaternion Identity
        {

            get
            {
                return new Quaternion(0, 0, 0, 1);
            }
        }
        /// <summary>
        /// Получает значение, указывающее, является ли текущий экземпляр единичным кватернионом.
        /// </summary>
        public bool IsIdentity
        {
            get
            {
                return this.X == 0.0 && this.Y == 0.0 && this.Z == 0.0 && this.W == 1.0;
            }
        }
        /// <summary>
        /// Возвращает значение модуля кватерниона.
        /// </summary>
        public double Abs
        {
            get
            {
                return Math.Sqrt(this.X * this.X + this.Y * this.Y + this.Z * this.Z + this.W * this.W);
            }
        }
        /// <summary>
        /// Вычисляет модуля кватерниона в квадрате.
        /// </summary>
        public double SquaredAbs
        {
            get
            {
                return this.X * this.X + this.Y * this.Y + this.Z * this.Z + this.W * this.W;
            }
        }
        #endregion

        #region Operations
        /// <summary>
        /// Делит каждую координату указанного кватерниона на его длину.
        /// </summary>
        public Quaternion Normalize
        {
            get
            {
                double norm = 1.0 / this.Abs;
                return new Quaternion(this.X * norm, this.Y * norm, this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Возвращает сопряженный объект заданного кватерниона.
        /// </summary>
        public Quaternion Conjugate
        {
            get
            {
                return new Quaternion(-this.X, -this.Y, -this.Z, this.W);
            }
        }
        /// <summary>
        /// Возвращает инверсный объект кватерниона.
        /// </summary>
        public Quaternion Inverse
        {
            get
            {
                double norm = 1.0 / this.SquaredAbs;
                return new Quaternion(-this.X * norm, -this.Y * norm, -this.Z * norm, this.W * norm);
            }
        }
        /// <summary>
        /// Создает новый кватернион на основе заданного значения нутации, прецессии и собственного вращения.
        /// </summary>
        /// <param name="yaw">Угол нутации вокруг оси Y в радианах</param>
        /// <param name="pitch">Угол прецессии вокруг оси X в радианах</param>
        /// <param name="roll">Угол собственного вращения вокруг оси Z в радианах</param>
        /// <returns>Кватернион</returns>
        public static Quaternion FromYPR(double yaw, double pitch, double roll)
        {
            double a = roll * 0.5;
            double b = Math.Sin(a);
            double c = Math.Cos(a);
            double d = pitch * 0.5;
            double e = Math.Sin(d);
            double f = Math.Cos(d);
            double g = yaw * 0.5;
            double h = Math.Sin(g);
            double i = Math.Cos(g);

            return new Quaternion(
                i * e * c + h * f * b,
                h * f * c - i * e * b,
                i * f * b - h * e * c,
                i * f * c + h * e * b);
        }
        /// <summary>
        /// Вычисляет скалярное произведение двух кватернионов.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static double Dot(Quaternion a, Quaternion b)
        {
            return a.X * b.X + a.Y * b.Y + a.Z * b.Z + a.W * b.W;
        }
        /// <summary>
        /// Выполняет интерполяцию между двумя кватернионами, используя сферическую линейную интерполяцию.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <param name="amount">Относительный вес второго кватерниона в интерполяции</param>
        /// <returns>Кватернион</returns>
        public static Quaternion Slerp(Quaternion a, Quaternion b, double amount)
        {
            double d, e, dot = Quaternion.Dot(a, b);
            bool flag = false;

            if (dot < 0.0)
            {
                flag = true;
                dot = -dot;
            }
            if (dot > 0.999999)
            {
                d = 1.0 - amount;
                e = (flag ? (-amount) : amount);
            }
            else
            {
                double f = Math.Acos(dot);
                double g = (1.0 / Math.Sin(f));
                d = Math.Sin(((1.0 - amount) * f)) * g;
                e = (flag ? ((-Math.Sin((amount * f))) * g) : (Math.Sin((amount * f)) * g));
            }

            return new Quaternion(
                d * a.X + e * b.X,
                d * a.Y + e * b.Y,
                d * a.Z + e * b.Z,
                d * a.W + e * b.W);
        }
        /// <summary>
        /// Выполняет линейную интерполяцию между двумя кватернионами на основе значения, указывающего взвешивание второго кватерниона.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <param name="amount">Относительный вес второго кватерниона в интерполяции</param>
        /// <returns>Кватернион</returns>
        public static Quaternion Lerp(Quaternion a, Quaternion b, double amount)
        {
            double f = 1.0 - amount;
            Quaternion quaternion3 = default(Quaternion);
            double dot = Dot(a, b);

            if (dot >= 0.0)
            {
                quaternion3.X = f * a.X + amount * b.X;
                quaternion3.Y = f * a.Y + amount * b.Y;
                quaternion3.Z = f * a.Z + amount * b.Z;
                quaternion3.W = f * a.W + amount * b.W;
            }
            else
            {
                quaternion3.X = f * a.X - amount * b.X;
                quaternion3.Y = f * a.Y - amount * b.Y;
                quaternion3.Z = f * a.Z - amount * b.Z;
                quaternion3.W = f * a.W - amount * b.W;
            }

            double norm = 1.0 / quaternion3.Abs;
            quaternion3.X *= norm;
            quaternion3.Y *= norm;
            quaternion3.Z *= norm;
            quaternion3.W *= norm;
            return quaternion3;
        }
        /// <summary>
        /// Сцепляет два кватерниона.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion Concatenate(Quaternion a, Quaternion b)
        {
            double x = b.X, y = b.Y, z = b.Z, w = b.W;
            double x2 = a.X, y2 = a.Y, z2 = a.Z, w2 = a.W;

            double e = y * z2 - z * y2;
            double f = z * x2 - x * z2;
            double c = x * y2 - y * x2;
            double d = x * x2 + y * y2 + z * z2;

            return new Quaternion(
                x * w2 + x2 * w + e,
                y * w2 + y2 * w + f,
                z * w2 + z2 * w + c,
                w * w2 - d);
        }
        #endregion

        #region Operators
        /// <summary>
        /// Обращает знак каждой координаты кватерниона.
        /// </summary>
        /// <param name="q">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator -(Quaternion q)
        {
            return new Quaternion(-q.X, -q.Y, -q.Z, -q.W);
        }
        /// <summary>
        /// Складывает каждый элемент в одном кватернионе с соответствующим элементом во втором кватернионе.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator +(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.X + b.X,
                a.Y + b.Y,
                a.Z + b.Z,
                a.W + b.W);
        }
        /// <summary>
        /// Вычитает каждый элемент во втором кватернионе из соответствующего элемента в первом кватернионе.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator -(Quaternion a, Quaternion b)
        {
            return new Quaternion(
                a.X - b.X,
                a.Y - b.Y,
                a.Z - b.Z,
                a.W - b.W);
        }
        /// <summary>
        /// Возвращает кватернион, являющийся результатом перемножения двух кватернионов.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator *(Quaternion a, Quaternion b)
        {
            double x = a.X;
            double y = a.Y;
            double z = a.Z;
            double w = a.W;
            double x2 = b.X;
            double y2 = b.Y;
            double z2 = b.Z;
            double w2 = b.W;
            double d = y * z2 - z * y2;
            double e = z * x2 - x * z2;
            double g = x * y2 - y * x2;
            double h = x * x2 + y * y2 + z * z2;

            return new Quaternion(
                 x * w2 + x2 * w + d,
                 y * w2 + y2 * w + e,
                 z * w2 + z2 * w + g,
                 w * w2 - h
                );
        }
        /// <summary>
        /// Возвращает кватернион, получаемый в результате масштабирования всех координат заданного кватерниона на скалярный множитель.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Множитель</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator *(Quaternion a, double b)
        {
            return new Quaternion(
                a.X * b,
                a.Y * b,
                a.Z * b,
                a.W * b);
        }
        /// <summary>
        /// Делит один кватернион на второй кватернион.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator /(Quaternion a, Quaternion b)
        {
            double x = a.X;
            double y = a.Y;
            double z = a.Z;
            double w = a.W;
            double d = 1.0 / b.SquaredAbs;
            double e = -b.X * d;
            double f = -b.Y * d;
            double g = -b.Z * d;
            double i = b.W * d;
            double j = y * g - z * f;
            double k = z * e - x * g;
            double l = x * f - y * e;
            double m = x * e + y * f + z * g;

            return new Quaternion(
                x * i + e * w + j,
                y * i + f * w + k,
                z * i + g * w + l,
                w * i - m);
        }
        /// <summary>
        /// Возвращает кватернион, получаемый в результате масштабирования всех координат заданного кватерниона на скалярный множитель.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Делитель</param>
        /// <returns>Кватернион</returns>
        public static Quaternion operator /(Quaternion a, double b)
        {
            return new Quaternion(
                a.X / b,
                a.Y / b,
                a.Z / b,
                a.W / b);
        }
        #endregion

        #region Bools & overrides
        /// <summary>
        /// Проверяет равны ли два объекта типа Кватернион между собой.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Логическое значение</returns>
        public static bool operator ==(Quaternion a, Quaternion b)
        {
            return a.X == b.X && a.Y == b.Y && a.Z == b.Z && a.W == b.W;
        }
        /// <summary>
        /// Проверяет не равны ли два объекта типа Кватернион между собой.
        /// </summary>
        /// <param name="a">Кватернион</param>
        /// <param name="b">Кватернион</param>
        /// <returns>Логическое значение</returns>
        public static bool operator !=(Quaternion a, Quaternion b)
        {
            return !(a == b);
        }
        /// <summary>
        /// Возвращает значение, указывающее, равен ли данный экземляр заданному значению типа Кватернион.
        /// </summary>
        /// <param name="obj">Объект</param>
        /// <returns>Логическое значение</returns>
        public override bool Equals(object obj)
        {
            return (obj is Quaternion) ? (this == (Quaternion)obj) : false;
        }
        /// <summary>
        /// Преобразует комплексное число в соответствующее ему строковое представление.
        /// </summary>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public override string ToString()
        {
            return this.ToString("G3");
        }
        /// <summary>
        /// Преобразует комплексное число в соответствующее ему строковое представление.
        /// </summary>
        /// <param name="format">Строка числового формата</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public string ToString(string format)
        {
            return StringOptions.Disp(new double[] { this.X, this.Y, this.Z, this.W }, format, StringOptions.Q);
        }
        /// <summary>
        /// Возвращает хэш-код для данного объекта.
        /// </summary>
        /// <returns>Целое число со знаком</returns>
        public override int GetHashCode()
        {
            return this.X.GetHashCode() + this.Y.GetHashCode() + this.Z.GetHashCode() + this.W.GetHashCode();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию кватерниона.
        /// </summary>
        /// <returns>Кватернион</returns>
        object ICloneable.Clone()
        {
            return new Quaternion(this.X, this.Y, this.Z, this.W);
        }
        /// <summary>
        /// Создает копию кватерниона.
        /// </summary>
        /// <returns>Кватернионо</returns>
        public Quaternion Clone()
        {
            return new Quaternion(this.X, this.Y, this.Z, this.W);
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
            info.AddValue("X", this.X);
            info.AddValue("Y", this.Y);
            info.AddValue("Z", this.Z);
            info.AddValue("W", this.W);
        }
        #endregion

        #region Parsing
        /// <summary>
        /// Переводит исходную строку в кватернион.
        /// </summary>
        /// <remarks>
        /// Пример входной строки: "[1, -2; 3.2, -.13]";
        /// Размерность вектора должна быть равна 4.
        /// </remarks>
        /// <param name="s">Исходная строка</param>
        /// <returns>Кватернион</returns>
        public static Quaternion Parse(string s)
        {
            string[] cols = StringOptions.Matpar(s);
            string[] nums = cols[0].Split('|');

            if (cols.Length > 1 || nums.Length != 4)
                throw new Exception("Входная строка имела неверный формат");

            return new Quaternion(double.Parse(nums[0]),
                                  double.Parse(nums[1]),
                                  double.Parse(nums[2]),
                                  double.Parse(nums[3]));
        }
        /// <summary>
        /// Пробует перевести исходную строку в кватернион.
        /// </summary>
        /// <param name="quaternion">Исходная строка</param>
        /// <param name="result">Кватернион</param>
        /// <returns>Логическое значение</returns>
        public static bool TryParse(string quaternion, ref Quaternion result)
        {
            try
            {
                result = Quaternion.Parse(quaternion);
                return true;
            }
            catch (FormatException)
            {
                result = new Quaternion();
                return false;
            }
        }
        #endregion
    }
    #endregion

    #region Internal class
    /// <summary>
    /// Определяет класс операций со строками.
    /// </summary>
    internal class StringOptions
    {
        #region String voids
        /// <summary>
        /// Формат комплексного числа.
        /// </summary>
        public static string[] C
        {
            get
            {
                return new string[] { "", "i" };
            }
        }
        /// <summary>
        /// Формат кватерниона.
        /// </summary>
        public static string[] Q
        {
            get
            {
                return new string[] { "i", "j", "k", "" };
            }
        }
        /// <summary>
        /// Функция преобразования массива чисел в строку.
        /// </summary>
        /// <param name="v">Массив чисел</param>
        /// <param name="format">Строка числового формата</param>
        /// <param name="symbol">Массив символов</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public static string Disp(double[] v, string format, string[] symbol)
        {
            // Переменые:
            int length = v.Length, i;
            int start = -1;

            // Вычисление индекса первого
            // ненулевого элемента:
            for (i = 0; i < length; i++)
            {
                if (v[i] != 0)
                {
                    start = i;
                    break;
                }
            }

            // Если был найден ненулевой элемент:
            if (start != -1)
            {
                // Отображение первого ненуевого элемента:
                string result = Disp(v[start], format, true, symbol[start]);

                for (i = start + 1; i < length; i++)
                {
                    result += Disp(v[i], format, false, symbol[i]);
                }
                return result;
            }
            return "0";
        }
        /// <summary>
        /// Функция преобразования числового значения в строку.
        /// </summary>
        /// <param name="v">Число</param>
        /// <param name="format">Строка числового формата</param>
        /// <param name="s">Первое в ряду или нет</param>
        /// <param name="symbol">Символ</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public static string Disp(double v, string format, bool s, string symbol)
        {
            if (v == 0)
            {
                return "";
            }
            else if (v < 0)
            {
                return (s) ? "-" + (-v).ToString(format) + symbol : " - " + (-v).ToString(format) + symbol;
            }
            return (s) ? v.ToString(format) + symbol : " + " + v.ToString(format) + symbol;
        }
        /// <summary>
        /// Определяет общий метод приведения исходной строки к матричному виду.
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <returns>Массив строк</returns>
        public static string[] Matpar(string s)
        {
            // example: s = "[1,2,3,4]".
            // Regex options:
            Regex regex = new Regex(@"\[(?<matrice>.*)]", RegexOptions.None);
            Match match = regex.Match(s);

            // success?
            if (match.Success)
            {
                // get new string:
                return Regex.Split(match.Result("${matrice}").Replace(",", "|"), ";");
            }
            throw new Exception("Входная строка имела неверный формат");
        }
        /// <summary>
        /// Переводит исходную строку в комплексное число.
        /// <remarks>
        /// Примеры входной строки: "1 + 2i", "0.321 + 11i", ".1i".
        /// </remarks>
        /// </summary>
        /// <param name="s">Исходная строка</param>
        /// <returns>Текст как последовательность знаков Юникода</returns>
        public static Complex Compar(string s)
        {
            string u = s.Replace(" ", ""); // приведение строки к каноничному виду,
            int i, k = 0;                  // переменные,
            int length = u.Length;         // новая длина строки,
            string re = "" + u[0];         // действительная часть числа,
            string im = "";                // мнимая часть числа.

            // Случай, когда строка не содержит мнимой единицы,
            // значит комплексное число не содержит мнимой части:
            if (!u.Contains("i"))
            {
                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                return new Complex(double.Parse(re), 0);
            }
            // случай, когда строка содержит мнимую единицу,
            // значит комплексное число содержит мнимую часть:
            else
            {
                // случай, когда строка содержит только мнимую единицу:
                if (u == "i") return new Complex(0, 1);

                for (i = 1; i < length; i++)
                {
                    if (u[i] != '+' && u[i] != '-')
                    {
                        re += u[i]; k++;
                    }
                    else break;
                }

                // случай, когда комплексное число содержит действительную часть и
                // мнимую часть: c = re + im * i.
                if (k != length - 1)
                {
                    int k1 = k + 1, k2 = k + 2; im += u[k1];

                    // случай, когда мнимая единица обозначена, только как 'i'.
                    if (u[k2] == 'i') return new Complex(double.Parse(re), double.Parse(im + '1'));

                    for (i = k2; i < length; i++)
                    {
                        if (u[i] != 'i')
                        {
                            im += u[i];
                        }
                        else break;
                    }
                    return new Complex(double.Parse(re), double.Parse(im));
                }
                else // случай, когда комплексное число содержит только мнимую часть: c = re * i.
                {
                    return new Complex(0, double.Parse(re.Replace("i", "")));
                }
            }
        }
        #endregion
    }
    #endregion
}
