// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2019
// Asiryan Valeriy
// Moscow, Russia
// Version 4.0.0

using System;
using System.Threading.Tasks;
using UMapx.Core;
using UMapx.Transform;

namespace UMapx.Window
{
    // **************************************************************************
    //                              WINDOW TOOLBOX
    //                            UMAPX.NET FRAMEWORK
    // **************************************************************************
    // Window Toolbox provides a wide functionality for the study discrete 
    // and continuous wavelets. It includes algorithms for discrete one- and two-
    // dimensional wavelet transforms of real and complex signals.
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
        /// <param name="a">Параметр формы [0, 0.5]</param>
        public Planck(int frameSize, double a = 0.15)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Получает или задает значение параметра формы [0, 0.5].
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = Maths.Range(value, 0, 0.5);
            }
        }
        /// <summary>
        /// Функция Z+-(x, a).
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="p">Знак</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private double Z(double x, bool p, int frameSize)
        {
            // params:
            double t = p ? 1 : -1;
            double y = 2 * x / (frameSize - 1) - 1;

            // function:
            double u = 1.0 / (1 + t * y);
            double v = 1.0 / (1 - 2 * a + t * y);
            return 2 * a * (u + v);
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Planck taper window:
            double n = frameSize - 1;
            double b = a * n;
            double c = (1 - a) * n;

            // Creating:
            if (x >= 0 && x < b)
            {
                return 1.0 / (Math.Exp(Z(x, true, frameSize)) + 1);
            }
            else if (x >= b && x <= c)
            {
                return 1.0;
            }
            else if (x > c && x <= n)
            {
                return 1.0 / (Math.Exp(Z(x, false, frameSize)) + 1);
            }
            return 0;
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        /// <param name="a">Параметр формы [0, 1]</param>
        public Tukey(int frameSize, double a = 1)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Получает или задает значение параметра формы [0, 1].
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = Maths.Double(value);
            }
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Tukey window:
            double n = frameSize - 1;
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
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Confined(int frameSize, double sigma = 1)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Инициализирует закрытую оконную функцию Гаусса.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        public Confined(int frameSize)
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
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Вычисление функции:
            double a = G(-0.5) * (G(x + frameSize) + G(   x - frameSize));
            double b = G(-0.5         + frameSize) + G(-0.5 - frameSize);
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
        public override double[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Normal(int frameSize, double sigma = 1, double pow = 2)
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
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (sigma * a);
            return Math.Exp(-Math.Pow(t, p));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Kaiser(int frameSize, double a = 3)
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
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Kaiser window:
            double u = 2 * x / (frameSize - 1);
            double r = 1 - u * u;
            double v = r >= 0 ? Math.Sqrt(1 - u * u) : 0;
            double z = Math.PI * this.a;
            double q = Special.I0(z * v);
            return q / Special.I0(z    );
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
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
        public Welch(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Welch function:
            double t = (frameSize - 1) / 2.0;
            double a = (x - t) / t;
            return 1 - a * a;
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Lanzcos(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Lanczos function:
            return Special.Sinc(2 * x / (frameSize - 1) - 1, Math.PI);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Parzen(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
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
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
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
        public FlatTop(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 1 - 1.93 * Cosine.cosinefunc(2 * x, frameSize) + 1.29 * Cosine.cosinefunc(4 * x, frameSize) - 0.388 * Cosine.cosinefunc(6 * x, frameSize) + 0.028 * Cosine.cosinefunc(8 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Nuttall(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.355768 - 0.487396 * Cosine.cosinefunc(2 * x, frameSize) + 0.144232 * Cosine.cosinefunc(4 * x, frameSize) - 0.012604 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public BlackmanNuttall(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.3635819 - 0.4891775 * Cosine.cosinefunc(2 * x, frameSize) + 0.1365995 * Cosine.cosinefunc(4 * x, frameSize) - 0.0106411 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            // window function on a discrete time:
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public BlackmanHarris(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.35875 - 0.48829 * Cosine.cosinefunc(2 * x, frameSize) + 0.14128 * Cosine.cosinefunc(4 * x, frameSize) - 0.01168 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Blackman(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.42 - 0.5 * Cosine.cosinefunc(2 * x, frameSize) + 0.08 * Cosine.cosinefunc(4 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public BarlettHann(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Berlett-Hann function:
            double a = Math.Abs(Math.Abs(x / (frameSize - 1)) - 0.5);
            return 0.62 - 0.48 * a - 0.38 * Cosine.cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Hann(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Pow(Sine.sinefunc(x, frameSize), 2);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Hamming(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.53836 - 0.46164 * Cosine.cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
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
        public Cosine(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Множитель окна</returns>
        internal static double cosinefunc(double x, int frameSize)
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
        public Sine(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        internal static double sinefunc(double x, int frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконную функцию Габора.
    /// </summary>
    public class Gabor : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Инициализирует оконную функцию Габора.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <param name="sigma">Параметр сжатия</param>
        public Gabor(int frameSize, double sigma = 1)
        {
            this.FrameSize = frameSize;
            this.Sigma = sigma;
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
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public override double Function(double x, int frameSize)
        {
            // Gabor window function
            double y = x / frameSize;
            double z = Math.Pow(2 * Math.PI * y / sigma, 2);
            return Math.Exp(-z);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
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
        protected int frameSize;
        #endregion

        #region Window components
        /// <summary>
        /// Получает или задает размер окна.
        /// </summary>
        public int FrameSize
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
        public double Function(double x)
        {
            return this.Function(x, this.frameSize);
        }
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <returns>Одномерный массив</returns>
        public double[] GetWindow()
        {
            return this.GetWindow(this.frameSize);
        }
        /// <summary>
        /// Возвращает массив значений оконной функции.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public double[] Function(double[] x, int frameSize)
        {
            int length = x.Length;
            double[] H = new double[length];

            for (int i = 0; i < length; i++)
            {
                H[i] = Function(x[i], frameSize);
            }

            return H;
        }
        /// <summary>
        /// Возвращает массив значений оконной функции.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Function(double[] x)
        {
            return this.Function(x, this.frameSize);
        }
        /// <summary>
        /// Возвращает значение оконной функции.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public abstract double Function(double x, int frameSize);
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        public abstract double[] GetWindow(int frameSize);
        #endregion
    }
    #endregion

    #region Short-time Fourier analysis
    /// <summary>
    /// Определяет быстрое оконное преобразование Фурье.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    public class FastShortTimeFourierTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FastFourierTransform FFT;
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Коэффициенты оконной функции.
        /// </summary>
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое оконное преобразование Фурье.
        /// </summary>
        /// <param name="function">Оконная функция</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public FastShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.FFT = new FastFourierTransform(normalized, direction);
            Direction = direction;
            Window = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-64);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-time Fourier transform
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            // params
            int N = A.Length, i, j;
            Complex[] B = new Complex[N];
            int frame = coefs.Length;

            // Short-Time Fourier Transform
            for (i = 0; i < N; i += frame)
            {
                Complex[] data = new Complex[frame];

                for (j = 0; j < frame; j++)
                    data[j] = A[i + j] * coefs[Maths.Mod(i - frame / 2, frame)];

                data = FFT.Forward(data);

                for (j = 0; j < frame; j++)
                    B[i + j] = data[j];
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length, i, j;
            Complex[] A = new Complex[N];
            int frame = coefs.Length;

            for (i = 0; i < N; i += frame)
            {
                Complex[] data = new Complex[frame];

                for (j = 0; j < frame; j++)
                {
                    data[j] = B[i + j];
                }

                data = FFT.Backward(data);

                for (j = 0; j < frame; j++)
                {
                    A[i + j] = data[j] / coefs[Maths.Mod(i - frame / 2, frame)];
                }
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            // Fourier transform:
            Complex[,] B = (Complex[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Определяет оконное преобразование Фурье.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    public class ShortTimeFourierTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private FourierTransform FFT;
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        /// <summary>
        /// Коэффициенты оконной функции.
        /// </summary>
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует оконное преобразование Фурье.
        /// </summary>
        /// <param name="function">Оконная функция</param>
        /// <param name="normalized">Нормализированное преобразование или нет</param>
        /// <param name="direction">Направление обработки</param>
        public ShortTimeFourierTransform(IWindow function, bool normalized = true, Direction direction = Direction.Vertical)
        {
            // fourier transform initialization:
            this.FFT = new FourierTransform(normalized, direction);
            Direction = direction;
            Window = function;

            // sampling window function:
            this.coefs = function.GetWindow().Add(1e-64);
        }
        /// <summary>
        /// Нормализированное преобразование или нет.
        /// </summary>
        public bool Normalized
        {
            get
            {
                return this.FFT.Normalized;
            }
            set
            {
                this.FFT.Normalized = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-time Fourier transform
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            // params
            int N = A.Length, i, j;
            Complex[] B = new Complex[N];
            int frame = coefs.Length;

            // Short-Time Fourier Transform
            for (i = 0; i < N; i += frame)
            {
                Complex[] data = new Complex[frame];

                for (j = 0; j < frame; j++)
                    data[j] = A[i + j] * coefs[Maths.Mod(i - frame / 2, frame)];

                data = FFT.Forward(data);

                for (j = 0; j < frame; j++)
                    B[i + j] = data[j];
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length, i, j;
            Complex[] A = new Complex[N];
            int frame = coefs.Length;

            for (i = 0; i < N; i += frame)
            {
                Complex[] data = new Complex[frame];

                for (j = 0; j < frame; j++)
                {
                    data[j] = B[i + j];
                }

                data = FFT.Backward(data);

                for (j = 0; j < frame; j++)
                {
                    A[i + j] = data[j] / coefs[Maths.Mod(i - frame / 2, frame)];
                }
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            // Fourier transform:
            Complex[,] B = (Complex[,])A.Clone();
            int N = A.GetLength(0);
            int M = A.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });

            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = Forward(col);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = Forward(row);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное оконное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            if (direction == Direction.Both)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }
                    
                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                // 2-d vertical short-time fourier transform:
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = Backward(col);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                // 2-d horizontal short-time fourier transform:
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }
                    
                    row = Backward(row);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    #endregion

    #region Gabor analysis
    /// <summary>
    /// Определяет группу ортогональных базисов и дискретных преобразований Вейля-Гейзенберга.
    /// <remarks>
    /// Базисы Вейля-Гейзенберга используются для получения частотно-временных характеристик сигнала.
    /// Более подробную информацию можно найти на сайте:
    /// https://elibrary.ru/item.asp?id=29767333
    /// </remarks>
    /// </summary>
    public class WeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Количество сдвигов по частоте.
        /// </summary>
        private int m;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует группу ортогональных базисов и преобразований Вейля-Гейзенберга.
        /// </summary>
        /// <param name="window">Оконная функция</param>
        /// <param name="m">Количество сдвигов по частоте [4, N]</param>
        /// <param name="direction">Направление обработки</param>
        public WeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Получает или задает количество сдвигов по частоте [4, N].
        /// <remarks>
        /// Четное число.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Неверное значение аргумента");

                this.m = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg static components
        /// <summary>
        /// Возвращает комплексную матрицу базиса Вейля-Гейзенберга.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="window">Оконная функция</param>
        /// <param name="N">Количество отсчетов, вмещающих функцию</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <param name="orthogonalize">Ортогонализированная матрица или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] WeylHeisenberg(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.GetPacket(window, N), M, orthogonalize);
        }
        /// <summary>
        /// Возвращает комплексную матрицу базиса Вейля-Гейзенберга.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Формирующая функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <param name="orthogonalize">Ортогонализированная матрица или нет</param>
        /// <returns>Матрица</returns>
        public static Complex[,] WeylHeisenberg(double[] g0, int M, bool orthogonalize = true)
        {
            if (orthogonalize)
            {
                return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.Zak(g0, M), M);
            }

            return WeylHeisenbergTransform.WeylHeisenberg(g0, M);
        }
        /// <summary>
        /// Возвращает вектор значений оконной функции.
        /// </summary>
        /// <param name="window">Оконная функция</param>
        /// <param name="length">Количество отсчетов, вмещаемых функцию</param>
        /// <returns>Одномерный массив</returns>
        public static double[] GetPacket(IWindow window, int length)
        {
            // exeption by length
            if (window.FrameSize > length)
                return WeylHeisenbergTransform.nSymmetry(window, length);

            // params for approximation
            double[] w = WeylHeisenbergTransform.nSymmetry(window, (int)window.FrameSize);
            int n = w.Length;
            double min = Math.Min(w[0], w[n - 1]);
            double[] g = new double[length];
            int i, j = (length - n) / 2;
            int k = Math.Min(length - 2 * j, n);
            int z = j + k;

            // do job for intervals
            for (i = 0; i < j; i++)
                g[i] = min;

            for (i = j; i < z; i++)
                g[i] = w[i - j];

            for (i = z; i < length; i++)
                g[i] = min;

            return g;
        }
        /// <summary>
        /// Возвращает вектор значений оконной функции, удовлетворяющей условию N-1 симметрии.
        /// </summary>
        /// <param name="window">Оконная функция</param>
        /// <param name="length">Количество отсчетов функции</param>
        /// <returns>Одномерный массив</returns>
        private static double[] nSymmetry(IWindow window, int length)
        {
            // creaing window function
            double[] g = window.GetWindow(length + 1);
            double[] w = new double[length];

            // N-1 symmetric
            for (int i = 0; i < length; i++)
                w[i] = g[i];

            return w;
        }
        /// <summary>
        /// Возвращает комплексную матрицу базиса Вейля-Гейзенберга.
        /// <remarks>
        /// Размерность матрицы [N, N], где N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Формирующая функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Матрица</returns>
        private static Complex[,] WeylHeisenberg(double[] g0, int M)
        {
            int N = g0.Length, L = N / M;                               // Определение параметров,

            if (L <= 0)
                throw new Exception("Количество сдвигов по частоте определено неверно");

            Complex[,] G = new Complex[N, N];                           // комплексный прямоугольный сигнальный базис,
            Complex c = 2 * Maths.Pi * Maths.I;                         // комплексный коэффициент преобразования,
            double a = M / 2.0;                                         // вычисление оптимального фазового параметра.

            Parallel.For(0, N, n =>                                     // Использование параллельных циклов.
            {
                double phase = n - a / 2.0;                             // Фазовый сдвиг.
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);                 // Экспоненциальный коэффициент.

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;                                  // Элемент матрицы,
                        i = Maths.Mod(n - l * M, N);                    // сдвиг по времени,
                        j = Maths.Mod(n + M / 2 - l * M, N);            // сдвиг по частоте,

                        psi = new Complex(
                            (g0[i] * exp).Re,                           // Комлпексный сигнальный базис,
                            (Maths.I * g0[j] * exp).Re);                // < Ψ Re, Ψ Im >

                        G[n, u] = psi;
                    }
                }
            });

            return G;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] B = Matrice.Dot(A, U.Hermitian());
            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] A = Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            int N = A.GetLength(0), M = A.GetLength(1);
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(this.window, M, this.m, true);
            Complex[,] B;

            if (direction == Direction.Both)
            {
                B = U.Hermitian().Dot(A).Dot(V);
            }
            else if (direction == Direction.Vertical)
            {
                B = U.Hermitian().Dot(A);
            }
            else
            {
                B = A.Dot(V);
            }
            return B;
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            int N = B.GetLength(0), M = B.GetLength(1);
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[,] V = WeylHeisenbergTransform.WeylHeisenberg(this.window, M, this.m, true);
            Complex[,] A;

            if (direction == Direction.Both)
            {
                A = U.Dot(B).Dot(V.Hermitian());
            }
            else if (direction == Direction.Vertical)
            {
                A = U.Dot(B);
            }
            else
            {
                A = B.Dot(V.Hermitian());
            }
            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Zak components
        /// <summary>
        /// Преобразование Фурье.
        /// </summary>
        private static FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        /// <summary>
        /// Быстрое преобразование Фурье.
        /// </summary>
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Vertical);
        /// <summary>
        /// Реализует Zak-ортогонализацию вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одномерный массив</returns>
        public static double[] Zak(double[] v, int M)
        {
            // Быстрый алгоритм ортогонализации формирующей 
            // WH-функции с использованием дискретного Zak-преобразования.
            // В.П. Волчков, Д.А. Петров и В.М. Асирян.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            double[] vort = new double[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            // Дискретное Zak-преобразование:
            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            // Быстрое или обычное 
            // преобразование Фурье по столбцам:
            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            // Параметры:
            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            // Вычисление формирующей WH-функции:
            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            // Обратное дискретное Zak-преобразование:
            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = (sum / L2).Re;
            }

            return vort;
        }
        /// <summary>
        /// Реализует Zak-ортогонализацию вектора.
        /// </summary>
        /// <param name="v">Одномерный массив</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одномерный массив</returns>
        public static Complex[] Zak(Complex[] v, int M)
        {
            // Быстрый алгоритм ортогонализации формирующей 
            // WH-функции с использованием дискретного Zak-преобразования.
            // В.П. Волчков, Д.А. Петров и В.М. Асирян.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            Complex[] vort = new Complex[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            // Дискретное Zak-преобразование:
            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            // Быстрое или обычное 
            // преобразование Фурье по столбцам:
            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            // Параметры:
            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

            // Вычисление формирующей WH-функции:
            for (i = 0; i < L; i++)
            {
                for (j = 0; j < N; j++)
                {
                    z1 = Z[i, j];
                    z2 = Z[L + i, j];

                    even = Math.Pow(z1.Abs, 2);
                    odd = Math.Pow(z2.Abs, 2);
                    phi = w / Math.Sqrt(even + odd);

                    Z[i, j] = z1 * phi;
                    Z[L + i, j] = z2 * phi;
                }
            }

            // Обратное дискретное Zak-преобразование:
            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = sum / L2;
            }

            return vort;
        }
        #endregion
    }
    /// <summary>
    /// Определяет класс быстрого преобразования Вейля-Гейзенберга.
    /// <remarks>
    /// Класс представляет вычислительно эффективную реализацию одномерного и двумерного дикретных ортогональных
    /// преобразований Вейля-Гейзенберга.
    /// Более подробную информацию можно найти на сайте:
    /// https://elibrary.ru/title_about.asp?id=58245
    /// </remarks>
    /// </summary>
    public class FastWeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Преобразование Фурье.
        /// <remarks>
        /// Используется как вспомогательный компонент.
        /// </remarks>
        /// </summary>
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Horizontal);
        /// <summary>
        /// Оконная функция.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Количество сдвигов по частоте.
        /// </summary>
        private int m;
        /// <summary>
        /// Направление обработки.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Инициализирует быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="window">Оконная функция</param>
        /// <param name="m">Количество сдвигов по частоте [4, N]</param>
        /// <param name="direction">Направление обработки</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Получает или задает количество сдвигов по частоте [4, N].
        /// <remarks>
        /// Четное число.
        /// </remarks>
        /// </summary>
        public int M
        {
            get
            {
                return this.m;
            }
            set
            {
                if (value <= 2 || !Maths.IsEven(value))
                    throw new Exception("Неверное значение аргумента");

                this.m = value;
            }
        }
        /// <summary>
        /// Получает или задает направление обработки.
        /// </summary>
        public Direction Direction
        {
            get
            {
                return this.direction;
            }
            set
            {
                this.direction = value;
            }
        }
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Прямое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Forward(Complex[] A)
        {
            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, A.Length), this.m);
            return FastWeylHeisenbergTransform.WHT(A, g0, m);
        }
        /// <summary>
        /// Обратное дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public Complex[] Backward(Complex[] B)
        {
            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, B.Length), this.m);
            return FastWeylHeisenbergTransform.IWHT(B, g0, m);
        }
        /// <summary>
        /// Прямое дискретное быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Forward(Complex[,] A)
        {
            Complex[,] B = (Complex[,])A.Clone();
            int N = B.GetLength(0), M = B.GetLength(1);

            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, N), this.m);
            double[] g1 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, M), this.m);

            if (direction == Direction.Both)
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });

                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = B[i, j];
                    }

                    col = WHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        B[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = B[i, j];
                    }

                    row = WHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        B[i, j] = row[j];
                    }
                });
            }

            return B;
        }
        /// <summary>
        /// Обратное дискретное быстрое преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public Complex[,] Backward(Complex[,] B)
        {
            Complex[,] A = (Complex[,])B.Clone();
            int N = B.GetLength(0);
            int M = B.GetLength(1);

            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, N), this.m);
            double[] g1 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, M), this.m);

            if (direction == Direction.Both)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });

                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }
            else if (direction == Direction.Vertical)
            {
                Parallel.For(0, M, j =>
                {
                    Complex[] col = new Complex[N];
                    int i;

                    for (i = 0; i < N; i++)
                    {
                        col[i] = A[i, j];
                    }

                    col = IWHT(col, g0, m);

                    for (i = 0; i < N; i++)
                    {
                        A[i, j] = col[i];
                    }
                });
            }
            else
            {
                Parallel.For(0, N, i =>
                {
                    Complex[] row = new Complex[M];
                    int j;

                    for (j = 0; j < M; j++)
                    {
                        row[j] = A[i, j];
                    }

                    row = IWHT(row, g1, m);

                    for (j = 0; j < M; j++)
                    {
                        A[i, j] = row[j];
                    }
                });
            }

            return A;
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Прямое дискретное преобразование Фурье.
        /// </summary>
        /// <param name="A">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Обратное дискретное преобразование Фурье.
        /// </summary>
        /// <param name="B">Двумерный массив</param>
        /// <returns>Двумерный массив</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Прямое быстрое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="g0">Формирующая WH-функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одноменый массив</returns>
        public static Complex[] WHT(Complex[] input, double[] g0, int M)
        {
            // Функция реализует быстрый алгоритм прямого преобразования Вейля-Гейзенберга, 
            // изложенный в следующих статьях:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // В.М. Асирян, В.П. Волчков, "ЭФФЕКТИВНАЯ РЕАЛИЗАЦИЯ ПРЯМОГО ПРЕОБРАЗОВАНИЯ ВЕЙЛЯ-ГЕЙЗЕНБЕРГА" [2].
            // Алгоритм является вычислительно эффективным при больших значениях M.

            // Параметры преобразования и выходной сигнал:
            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);

            // Вспомогательные матрицы и переменные:
            Complex[,] s0 = new Complex[M, L];
            Complex[,] a0 = new Complex[M, L];
            Complex[,] b0 = new Complex[M, L];
            Complex[,] A0 = new Complex[L, M];
            Complex[,] B0 = new Complex[L, M];
            Complex[,] A1 = new Complex[L, M2];
            Complex[,] B1 = new Complex[L, M2];
            Complex c1re, c2re;
            Complex c1im, c2im;
            int k, i, j, u, n, m, l;

            // 1. Формирование матриц перестановок:
            for (m = 0; m < M; m++)
            {
                for (n = 0; n < L; n++)
                {
                    u = n * M;
                    i = Maths.Mod(m + M4 + u, N);
                    j = Maths.Mod(m - M4 + u, N);
                    k = Maths.Mod(-m - M4 - u, N);

                    s0[m, n] = input[k];
                    a0[m, n] = g0[i];
                    b0[m, n] = g0[j];
                }
            }

            // 2. Вычисление циклической свертки для матриц
            // полученных в пункте 1:
            for (l = 0; l < L; l++)
            {
                for (n = 0; n < L; n++)
                {
                    k = Maths.Mod(n - l, L);

                    for (m = 0; m < M; m++)
                    {
                        A0[l, m] += a0[m, n] * s0[m, k];
                        B0[l, m] += b0[m, n] * s0[m, k];
                    }
                }
            }

            // Переменные замены:
            Complex x, y, z, w;

            // 3. Вычисление новых последовательностей:
            for (l = 0; l < L; l++)
            {
                for (m = 0; m < M2; m++)
                {
                    // Замена переменных для действительной части:
                    x = A0[l, m];
                    y = A0[l, m + M2];
                    z = A0[l, M2 - m].Conjugate;
                    w = A0[l, Maths.Mod(M - m, M)].Conjugate;

                    // Вычисление выражений:
                    c1re = x + y + z + w;
                    c2re = x - y - z + w;

                    // Аналогичная замена переменных для мнимой части:
                    x = B0[l, m];
                    y = B0[l, m + M2];
                    z = B0[l, M2 - m].Conjugate;
                    w = B0[l, Maths.Mod(M - m, M)].Conjugate;

                    // Вычисление выражений:
                    c1im = x + y - z - w;
                    c2im = x - y + z - w;

                    // Вычисление элементов матриц:
                    A1[l, m] = 1.0 / (2) * (c1re + Maths.I * c2re * exp[m]);
                    B1[l, m] = 1.0 / (2 * Maths.I) * (c1im + Maths.I * c2im * exp[m]);
                }
            }

            // 4. Быстрое обратное M/2- точечное преобразование Фурье по строкам:
            A1 = FFT.Backward(Matrice.Conjugate(A1));
            B1 = FFT.Backward(Matrice.Conjugate(B1));

            // 5. Формирование комплексного вектора сигнала:
            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Четные и нечетные индексы фильтра:
                    i = l * M + 2 * k;
                    j = l * M + 2 * k + 1;

                    // Замена переменных:
                    x = A1[l, k];
                    y = B1[l, k];

                    // Вычисление элементов вектора:
                    output[i] = x.Re + Maths.I * y.Re;
                    output[j] = x.Im + Maths.I * y.Im;
                }
            }

            // Результат прямого WH-преобразования.
            return output;
        }
        /// <summary>
        /// Обратное быстрое дискретное преобразование Вейля-Гейзенберга.
        /// </summary>
        /// <param name="input">Одномерный массив</param>
        /// <param name="g0">Формирующая WH-функция</param>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одноменый массив</returns>
        public static Complex[] IWHT(Complex[] input, double[] g0, int M)
        {
            // Функция реализует быстрый алгоритм обратного преобразования Вейля-Гейзенберга, 
            // изложенный в следующих статьях:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // В.М. Асирян, В.П. Волчков, "ЭФФЕКТИВНАЯ РЕАЛИЗАЦИЯ ПРЯМОГО ПРЕОБРАЗОВАНИЯ ВЕЙЛЯ-ГЕЙЗЕНБЕРГА" [2].
            // Алгоритм является вычислительно эффективным при больших значениях M.

            // Вспомогательные матрицы и переменные:
            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[,] A1 = new Complex[L, M];
            Complex[,] B1 = new Complex[L, M];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);
            Complex s;
            int n, k, l;

            // 1. Формирование вещественных матриц сдвигов:
            for (k = 0; k < M; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Замена переменных:
                    s = input[k + l * M];

                    // Вычисление элементов матриц:
                    A1[l, k] = s.Re;
                    B1[l, k] = s.Im;
                }
            }

            // 2. Вычисление матриц:
            Complex[,] Za = new Complex[L, M2];
            Complex[,] Zb = new Complex[L, M2];

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    Za[l, k] = A1[l, k * 2] + Maths.I * A1[l, k * 2 + 1];
                    Zb[l, k] = B1[l, k * 2] + Maths.I * B1[l, k * 2 + 1];
                }
            }

            // 3. Быстрое обратное M/2- точечное преобразование Фурье.
            Za = Matrice.Conjugate(FFT.Backward(Za));
            Zb = Matrice.Conjugate(FFT.Backward(Zb));

            // Переменные замены:
            Complex a0, a1, b0, b1;
            Complex x, y, u, v;

            // 4. Формирование новых матриц:
            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    // Для матрицы A(l, k)
                    // Замена переменных:
                    a0 = Za[l, k]; a1 = Za[l, Maths.Mod(M - k, M2)].Conjugate;

                    // Вычисление выражений:
                    x = 1.0 / (2) * (a0 + a1);
                    y = 1.0 / (2 * Maths.I) * (a0 - a1);
                    y *= exp[k];

                    // Формирование матрицы:
                    A1[l, k] = x + y;
                    A1[l, k + M2] = x - y;


                    // Для матрицы B(l, k)
                    // Замена переменных:
                    b0 = Zb[l, k]; b1 = Zb[l, Maths.Mod(M - k, M2)].Conjugate;

                    // Вычисление выражений:
                    u = 1.0 / (2) * (b0 + b1);
                    v = 1.0 / (2 * Maths.I) * (b0 - b1);
                    v *= exp[k];

                    // Формирование матрицы:
                    B1[l, k] = u + v;
                    B1[l, k + M2] = u - v;
                }
            }

            // 5. Формирование выходного сигнала:
            for (l = 0; l < L; l++)
            {
                for (n = 0; n < N; n++)
                {
                    // Вычисление выражения:
                    output[n] += A1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M, N)] - Maths.I
                               * B1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M + M2, N)];
                }
            }

            // Результат обратного WH-преобразования.
            return output;
        }
        /// <summary>
        /// Возвращает одномерный массив фазовых поворотов.
        /// </summary>
        /// <param name="M">Количество сдвигов по частоте</param>
        /// <returns>Одномерный массив</returns>
        private static Complex[] GetRotation(int M)
        {
            int M2 = M / 2;
            Complex[] phase = new Complex[M2];
            for (int k = 0; k < M2; k++)
            {
                phase[k] = Maths.Exp(Maths.I * 2 * Math.PI / M * k);
            }
            return phase;
        }
        #endregion
    }
    #endregion

    #region Window interfaces
    /// <summary>
    /// Определяет общий вид оконных функций.
    /// </summary>
    public interface IWindow
    {
        #region Interface
        /// <summary>
        /// Получает или задает размер окна.
        /// </summary>
        int FrameSize { get; set; }
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
        /// <summary>
        /// Возвращает оконную функцию.
        /// </summary>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        double[] GetWindow(int frameSize);
        /// <summary>
        /// Возвращает массив значений оконной функции.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <param name="frameSize">Размер окна</param>
        /// <returns>Одномерный массив</returns>
        double[] Function(double[] x, int frameSize);
        /// <summary>
        /// Возвращает массив значений оконной функции.
        /// </summary>
        /// <param name="x">Одномерный массив</param>
        /// <returns>Одномерный массив</returns>
        double[] Function(double[] x);
        #endregion
    }
    /// <summary>
    /// Определяет общий интерфейс окнонных преобразований.
    /// </summary>
    public interface IWindowTransform
    {
        #region Interface
        /// <summary>
        /// Получает или задает оконную функцию.
        /// </summary>
        IWindow Window { get; set; }
        #endregion
    }
    #endregion
}
