// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Runtime.Serialization;
using UMapx.Core;

namespace UMapx.Response
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                               UMAPX.RESPONSE
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Filter with response impulse reaction
    /// <summary>
    /// Определяет фильтр с бесконечной импульсной характеристикой.
    /// <remarks>
    /// Фильтр с бесконечной импульсной характеристикой (рекурсивный фильтр, БИХ-фильтр или IIR-фильтр) — линейный электронный фильтр, 
    /// использующий один или более своих выходов в качестве входа, то есть образует обратную связь. Основным свойством таких фильтров 
    /// является то, что их импульсная переходная характеристика имеет бесконечную длину во временной области, а передаточная функция 
    /// имеет дробно-рациональный вид.
    /// </remarks>
    /// </summary>
    public class IIR : IResponse, ICloneable, ISerializable
    {
        #region Private data
        private double[] a;
        private double[] b;
        #endregion

        #region IIR Components
        /// <summary>
        /// Инициализирует фильтр с бесконечной импульсной характеристикой.
        /// </summary>
        public IIR() { }
        /// <summary>
        /// Инициализирует фильтр с бесконечной импульсной характеристикой.
        /// </summary>
        /// <param name="a">Массив коэффициентов обратной связи</param>
        /// <param name="b">Массив коэффициентов сигнала</param>
        public IIR(double[] a, double[] b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Получает или задает массив коэффициентов обратной связи.
        /// </summary>
        public double[] A
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
        /// Получает или задает массив коэффициентов сигнала.
        /// </summary>
        public double[] B
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
        /// Возвращает массив значений реакции фильтра, при подаче дискретной функции.
        /// </summary>
        /// <param name="u">Одномерный массив</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Reaction(double[] u)
        {
            int length = u.Length;
            double[] y = new double[length];

            double input, output;
            int t, P = b.Length, Q = a.Length;
            int n, i, k;

            for (n = 0; n < length; n++)
            {
                input = 0; output = 0;

                for (i = 0; i < P; i++)
                {
                    t = n - i;
                    if (t < 0) continue;
                    input += b[i] * u[t];
                }

                for (k = 1; k < Q; k++)
                {
                    t = n - k;
                    if (t < 0) continue;
                    output += a[k] * y[t];
                }

                y[n] = input + output;

            }
            return y;
        }
        /// <summary>
        /// Возвращает амплитудно-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Amplitude(double[] w)
        {
            int i, j;
            Complex K1;
            Complex K2;
            double[] amplitude = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;
                K2 = Complex.One;

                for (i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                amplitude[j] = K1.Abs / K2.Abs;
            }
            return amplitude;
        }
        /// <summary>
        /// Возвращает фазо-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Phase(double[] w)
        {
            int i, j;
            Complex K1;
            Complex K2;
            double[] phase = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;
                K2 = Complex.One;

                for (i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i); }
                for (i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w[j] * i); }

                phase[j] = K1.Angle - K2.Angle;
            }
            return phase;
        }
        /// <summary>
        /// Возвращает значение амплитуды на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Amplitude(double w)
        {
            int i;
            Complex K1 = new Complex(0, 0);
            Complex K2 = new Complex(1, 0);

            for (i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Abs / K2.Abs;
        }
        /// <summary>
        /// Возвращает значение фазы на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Phase(double w)
        {
            int i;
            Complex K1 = new Complex(0, 0);
            Complex K2 = new Complex(1, 0);

            for (i = 0; i < b.Length; i++) { K1 += b[i] * Maths.Exp(-Maths.I * w * i); }
            for (i = 0; i < a.Length; i++) { K2 -= a[i] * Maths.Exp(-Maths.I * w * i); }

            return K1.Angle - K2.Angle;
        }
        /// <summary>
        /// Проверяет является ли заданный фильтр устойчивым.
        /// </summary>
        public bool Stability
        {
            get
            {
                Complex sum = Complex.Zero;

                for (int j = 0; j < a.Length; j++) { sum += a[j] * Maths.Exp(-Maths.I); }

                if (sum == Complex.Zero)
                {
                    return true;
                }
                return false;
            }
        }
        #endregion

        #region Sample filters
        /// <summary>
        /// Получает готовый фильтр нижних частот.
        /// </summary>
        public static IIR LowPass
        {
            get
            {
                IIR cv = new IIR();
                cv.A = new double[3] { 0, 0.5, 0.5 };
                cv.B = new double[3] { 1, 1, 0 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый фильтр верхних частот.
        /// </summary>
        public static IIR HighPass
        {
            get
            {
                IIR cv = new IIR();
                cv.A = new double[3] { 0, 0.5, 0.5 };
                cv.B = new double[3] { 1, -1, 0 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый полосовой фильтр.
        /// </summary>
        public static IIR BandPass
        {
            get
            {
                IIR cv = new IIR();
                cv.A = new double[3] { 0, 0.5, -0.5 };
                cv.B = new double[3] { 1, 1, -1 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый режекторный фильтр.
        /// </summary>
        public static IIR Notch
        {
            get
            {
                IIR cv = new IIR();
                cv.A = new double[3] { 0, 0.5, 0.5 };
                cv.B = new double[3] { 1, -1, 1 };
                return cv;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию фильтра.
        /// </summary>
        /// <returns>Фильтр</returns>
        object ICloneable.Clone()
        {
            return new IIR(this.A, this.B);
        }
        /// <summary>
        /// Создает копию фильтра.
        /// </summary>
        /// <returns>Фильтр</returns>
        public IIR Clone()
        {
            return new IIR(this.A, this.B);
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
            info.AddValue("A", this.A);
            info.AddValue("B", this.B);
        }
        #endregion
    }
    /// <summary>
    /// Определяет фильтр с конечной импульсной характеристикой.
    /// <remarks>
    /// Фильтр с конечной импульсной характеристикой (трансверсальный фильтр, КИХ-фильтр или FIR-фильтр) — один из видов линейных 
    /// цифровых фильтров, характерной особенностью которого является ограниченность по времени его импульсной характеристики 
    /// (с какого-то момента времени она становится точно равной нулю). Такой фильтр называют ещё нерекурсивным из-за отсутствия обратной связи. 
    /// Знаменатель передаточной функции такого фильтра — некая константа.
    /// </remarks>
    /// </summary>
    public class FIR : IResponse, ICloneable, ISerializable
    {
        #region Private data
        private double[] b;
        #endregion

        #region FIR Components
        /// <summary>
        /// Инициализирует фильтр с конечной импульсной характеристикой.
        /// </summary>
        public FIR() { }
        /// <summary>
        /// Инициализирует фильтр с конечной импульсной характеристикой.
        /// </summary>
        /// <param name="b">Массив коэффициентов сигнала</param>
        public FIR(double[] b)
        {
            B = b;
        }
        /// <summary>
        /// Получает или задает массив коэффициентов сигнала.
        /// </summary>
        public double[] B
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
        /// Возвращает массив значений реакции фильтра, при подаче дискретной функции.
        /// </summary>
        /// <param name="u">Одномерный массив</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Reaction(double[] u)
        {
            int length = u.Length;
            double[] y = new double[length];

            double input;
            int t, P = b.Length;
            int n, i;

            for (n = 0; n < length; n++)
            {
                input = 0;

                for (i = 0; i < P; i++)
                {
                    t = n - i;
                    if (t < 0) continue;
                    input += b[i] * u[t];
                }

                y[n] = input;

            }
            return y;
        }
        /// <summary>
        /// Возвращает амплитудно-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Amplitude(double[] w)
        {
            int i, j, length = b.Length;
            Complex K1;
            double[] amplitude = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;

                for (i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                amplitude[j] = K1.Abs;
            }
            return amplitude;
        }
        /// <summary>
        /// Возвращает фазо-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        public double[] Phase(double[] w)
        {
            int j, i, length = b.Length;
            Complex K1;
            double[] phase = new double[w.Length];

            for (j = 0; j < w.Length; j++)
            {
                K1 = Complex.Zero;

                for (i = 0; i < length; i++)
                {
                    K1 += b[i] * Maths.Exp(-Maths.I * w[j] * i);
                }
                phase[j] = K1.Angle;
            }
            return phase;
        }
        /// <summary>
        /// Возвращает значение амплитуды на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Amplitude(double w)
        {
            int i;
            int length = b.Length;
            Complex K1 = new Complex(0, 0);

            for (i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Abs;
        }
        /// <summary>
        /// Возвращает значение фазы на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Phase(double w)
        {
            int i;
            int length = b.Length;
            Complex K1 = new Complex(0, 0);

            for (i = 0; i < length; i++)
            {
                K1 += b[i] * Maths.Exp(-Maths.I * w * i);
            }
            return K1.Angle;
        }
        #endregion

        #region Sample filters
        /// <summary>
        /// Получает готовый фильтр нижних частот.
        /// </summary>
        public static FIR LowPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new double[3] { 1, 1, 0 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый фильтр верхних частот.
        /// </summary>
        public static FIR HighPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new double[3] { 1, -1, 0 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый полосовой фильтр.
        /// </summary>
        public static FIR BandPass
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new double[3] { 1, 1, -1 };
                return cv;
            }
        }
        /// <summary>
        /// Получает готовый режекторный фильтр.
        /// </summary>
        public static FIR Notch
        {
            get
            {
                FIR cv = new FIR();
                cv.B = new double[3] { 1, -1, 1 };
                return cv;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию фильтра.
        /// </summary>
        /// <returns>Фильтр</returns>
        object ICloneable.Clone()
        {
            return new FIR(this.B);
        }
        /// <summary>
        /// Создает копию фильтра.
        /// </summary>
        /// <returns>Фильтр</returns>
        public FIR Clone()
        {
            return new FIR(this.B);
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
            info.AddValue("B", this.B);
        }
        #endregion
    }
    #endregion

    #region Interfaces
    /// <summary>
    /// Определяет общий интерфейс фильтров отклика.
    /// </summary>
    public interface IResponse
    {
        #region Interface Components
        /// <summary>
        /// Возвращает массив значений реакции фильтра, при подаче дискретной функции.
        /// </summary>
        /// <param name="u">Одномерный массив</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        double[] Reaction(double[] u);
        /// <summary>
        /// Возвращает амплитудно-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        double[] Amplitude(double[] w);
        /// <summary>
        /// Возвращает фазо-частотную характеристику фильтра.
        /// </summary>
        /// <param name="w">Массив частот (рад/с)</param>
        /// <returns>Дискретная функция в декартовой системе координат</returns>
        double[] Phase(double[] w);
        /// <summary>
        /// Возвращает значение амплитуды на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Amplitude(double w);
        /// <summary>
        /// Возвращает значение фазы на заданной частоте.
        /// </summary>
        /// <param name="w">Частота (рад/с)</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Phase(double w);
        #endregion
    }
    #endregion
}
