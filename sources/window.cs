// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using UMapx.Core;

namespace UMapx.Window
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                                UMAPX.WINDOW
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Window functions
    /// <summary>
    /// Определяет оконную функцию Планка.
    /// </summary>
    public class Planck : WindowBase
    {
        #region Private data
        private double a = 0.15;
        #endregion

        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Планка.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="a">Параметр формы</param>
        public Planck(double frameSize, double a = 0.15)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Получает или задает значение параметра формы.
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
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Planck taper window:
            double n = this.frameSize - 1;
            double b = a * n;
            double c = (1 - a) * n;

            // Creating:
            if (x >= 0 && x < b)
            {
                return 1.0 / (Math.Exp(Z(x, true)) + 1);
            }
            else if (x >= b && x <= c)
            {
                return 1.0;
            }
            else if (x > c && x <= n)
            {
                return 1.0 / (Math.Exp(Z(x, false)) + 1);
            }
            return 0;
        }
        /// <summary>
        /// Функция Z+-(x, a).
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="p">Знак</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private double Z(double x, bool p)
        {
            // params:
            double t = p ? 1 : -1;
            double y = 2 * x / (this.frameSize - 1) - 1;

            // function:
            double u = 1.0 / (1 +         t * y);
            double v = 1.0 / (1 - 2 * a + t * y);
            return 2 * a * (u + v);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Тьюки.
    /// </summary>
    public class Tukey : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Тьюки.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="a">Параметр формы</param>
        public Tukey(double frameSize, double a = 1)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Получает или задает значение параметра формы.
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
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Tukey window:
            double n = this.frameSize - 1;
            double d = n * (1 - a / 2.0);
            double b = n / 2.0;
            double c = a * b;

            // Creating:
            if (x >= 0 && x < c)
            {
                return 0.5 * (1 + Math.Cos(Math.PI * (x / c - 1)));
            }
            else if (x >= c && x <= d)
            {
                return 1.0;
            }
            else if (x > d && x <= n)
            {
                return 0.5 * (1 + Math.Cos(Math.PI * (x / c - 2.0 / a + 1)));
            }
            return 0;
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет закрытую оконную функцию Гаусса.
    /// </summary>
    public class Confined : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Инициализирует закрытую оконную функцию Гаусса.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="sigma">Среднеквадратическое отклонение (0.14 * N)</param>
        public Confined(double frameSize, double sigma = 1)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Инициализирует закрытую оконную функцию Гаусса.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Confined(double frameSize)
        {
            this.Sigma = 0.14 * frameSize;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (>0).
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (sigma <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Вычисление функции:
            double a = G(-0.5) * (G(x + this.frameSize) + G(   x - this.frameSize));
            double b = G(-0.5         + this.frameSize) + G(-0.5 - this.frameSize);
            return G(x) - a / b;
        }
        /// <summary>
        /// Функция G(x).
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private double G(double x)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (2 * sigma);
            return Math.Exp(-t * t);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет обобщенную оконную нормальную функцию.
    /// </summary>
    public class Normal : WindowBase
    {
        #region Private data
        private double sigma = 1;
        private double p = 2;
        #endregion

        #region Window components
        /// <summary>
        /// Инициализирует обобщенную оконную нормальную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="sigma">Среднеквадратическое отклонение (>0)</param>
        /// <param name="pow">Степень<remarks>При p = 2 - окно Гаусса</remarks></param>
        public Normal(double frameSize, double sigma = 1, double pow = 2)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
            this.p = pow;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения (>0).
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (sigma <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Степень.
        /// </summary>
        public double Pow
        {
            get
            {
                return this.p;
            }
            set
            {
                this.p = value;
            }
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (sigma * a);
            return Math.Exp(-Math.Pow(t, p));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Кайзера.
    /// </summary>
    public class Kaiser : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion
        
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Кайзера.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="a">Параметр формы</param>
        public Kaiser(double frameSize, double a = 3)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Получает или задает значение параметра формы.
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
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Kaiser window:
            double u = 2 * x / (this.frameSize - 1);
            double v = Math.Sqrt(1 - u * u);
            double z = Math.PI * this.a;
            double q = Special.I0(z * v);
            return q / Special.I0(z    );
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1) / 2;
            double[] x = Matrice.Compute(-t, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию "Welch".
    /// </summary>
    public class Welch : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию "Welch".
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Welch(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Welch function:
            double t = (this.frameSize - 1) / 2;
            double a = (x - t) / t;
            return 1 - a * a;
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Ланцоша.
    /// </summary>
    public class Lanzcos : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Ланцоша.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Lanzcos(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Lanczos function:
            return Special.Sinc(2 * x / (this.frameSize - 1) - 1, Math.PI);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Парзена.
    /// </summary>
    public class Parzen : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Парзена.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Parzen(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // coefficients:
            double y = Math.Abs(x);
            double c = frameSize / 2.0;
            double a = x / c;
            double b = y / c;

            // props:
            if ((y >= 0) &&
                (y <= frameSize / 4))
            {
                return 1.0 - 6.0 * a * a * (1.0 - b);
            }
            else if (y >= frameSize / 4 &&
                y <= frameSize / 2)
            {
                return 2 * Math.Pow(1 - b, 3);
            }
            return 0.0;
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию "Flat-Top".
    /// </summary>
    public class FlatTop : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию "Flat-Top".
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public FlatTop(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 1 - 1.93 * Cosine.Function(2 * x, frameSize) + 1.29 * Cosine.Function(4 * x, frameSize) - 0.388 * Cosine.Function(6 * x, frameSize) + 0.028 * Cosine.Function(8 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Ньюттолла.
    /// </summary>
    public class Nuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Ньюттолла.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Nuttall(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 0.355768 - 0.487396 * Cosine.Function(2 * x, frameSize) + 0.144232 * Cosine.Function(4 * x, frameSize) - 0.012604 * Cosine.Function(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Блэкмана-Ньюттолла.
    /// </summary>
    public class BlackmanNuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Блэкмана-Ньюттолла.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public BlackmanNuttall(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 0.3635819 - 0.4891775 * Cosine.Function(2 * x, frameSize) + 0.1365995 * Cosine.Function(4 * x, frameSize) - 0.0106411 * Cosine.Function(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Блэкмана-Харриса.
    /// </summary>
    public class BlackmanHarris : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Блэкмана-Харриса.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public BlackmanHarris(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 0.35875 - 0.48829 * Cosine.Function(2 * x, frameSize) + 0.14128 * Cosine.Function(4 * x, frameSize) - 0.01168 * Cosine.Function(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Блэкмана.
    /// </summary>
    public class Blackman : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Блэкмана.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Blackman(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 0.42 - 0.5 * Cosine.Function(2 * x, frameSize) + 0.08 * Cosine.Function(4 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Барлетта-Ханна.
    /// </summary>
    public class BarlettHann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Барлетта-Ханна.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public BarlettHann(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            // Berlett-Hann function:
            double a = Math.Abs(Math.Abs(x / (this.frameSize - 1)) - 0.5);
            return 0.62 - 0.48 * a - 0.38 * Cosine.Function(2 * x, this.frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Ханна (Хеннинга).
    /// </summary>
    public class Hann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Ханна (Хеннинга).
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Hann(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return Math.Pow(Sine.Function(x, frameSize), 2);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Хэмминга.
    /// </summary>
    public class Hamming : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Хэмминга.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Hamming(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return 0.53836 - 0.46164 * Cosine.Function(2 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion
    }
    /// <summary>
    /// Определяет косинусную оконную функцию.
    /// </summary>
    public class Cosine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует косинусную оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Cosine(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return Math.Cos(Math.PI * x / (this.frameSize - 1));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            double t = (this.frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Множитель окна</returns>
        public static double Function(double x, double frameSize)
        {
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
    /// <summary>
    /// Определяет синусную оконную функцию.
    /// </summary>
    public class Sine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Инициализирует синусную оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Sine(double frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x)
        {
            return Math.Sin(Math.PI * x / (this.frameSize - 1));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow()
        {
            // window function on a discrete time:
            double t = this.frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return Matrice.Compute(x, Function);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public static double Function(double x, double frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
    /// <summary>
    /// Определяет общий класс для оконных функций.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Window_function
    /// </remarks>
    /// </summary>
    public abstract class WindowBase : IWindow
    {
        #region Private data
        /// <summary>
        /// Размер окна.
        /// </summary>
        protected double frameSize;
        #endregion

        #region Window components
        /// <summary>
        /// Получает или задает размер окна.
        /// </summary>
        public double FrameSize
        {
            get
            {
                return this.frameSize;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Размер окна должен быть строго больше 0");

                this.frameSize = value;
            }
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public abstract double Function(double x);
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public abstract double[] GetWindow();
        #endregion
    }
    /// <summary>
    /// Определяет общий вид оконных функций.
    /// </summary>
    public interface IWindow
    {
        #region Interface
        /// <summary>
        /// Получает или задает размер окна.
        /// </summary>
        double FrameSize { get; set; }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Function(double x);
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        double[] GetWindow();
        #endregion
    }
    #endregion
}
