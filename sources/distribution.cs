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

namespace UMapx.Distribution
{
    // **************************************************************************
    //                            UMAPX.NET FRAMEWORK
    //                             UMAPX.DISTRIBUTION
    // **************************************************************************
    // Designed by Asiryan Valeriy (c), 2015-2019
    // Moscow, Russia.
    // **************************************************************************

    #region Distributions
    /// <summary>
    /// Определяет распределение Гаусса.
    /// <remarks>
    /// Нормальное распределение, также называемое распределением Гаусса или Гаусса — Лапласа. 
    /// Носитель x ∈ (-inf, +inf), параметры: μ - математическое ожидание, σ > 0 - среднеквадратическое отклонение.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Normal_distribution
    /// </remarks>
    /// </summary>
    public class Gaussian : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double sigma = 1;
        private double mu = 0;
        #endregion;

        #region Gaussian components
        /// <summary>
        /// Инициализирует распределение Гаусса.
        /// </summary>
        public Gaussian() { }
        /// <summary>
        /// Инициализирует распределение Гаусса.
        /// </summary>
        /// <param name="sigma">Среднеквадратическое отклонение</param>
        /// <param name="mu">Коэффициент сдвига</param>
        public Gaussian(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента сдвига.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return sigma * sigma;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return this.mu;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return Maths.Exp(Maths.Pow((x - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            return 0.5 + 0.5 * Special.Erf((x - mu) / Maths.Sqrt(2.0 * sigma * sigma));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(sigma * Maths.Sqrt(2 * Maths.Pi * Maths.E), Maths.E);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Gaussian(sigma, mu);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Gaussian Clone()
        {
            return new Gaussian(sigma, mu);
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
            info.AddValue("Sigma", sigma);
            info.AddValue("Mu", mu);
        }
        #endregion
    }
    /// <summary>
    /// Определяет логарифмическое распределение Гаусса.
    /// <remarks>
    /// Двухпараметрическое семейство абсолютно непрерывных распределений. 
    /// Если случайная величина имеет логнормальное распределение, то её логарифм имеет нормальное распределение.
    /// Носитель x ∈ (0, +inf), параметры: μ - математическое ожидание, σ >= 0 - среднеквадратическое отклонение.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Log-normal_distribution
    /// </remarks>
    /// </summary>
    public class LogGaussian : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double sigma = 1;
        private double mu = 0;
        #endregion;

        #region GaussianLog components
        /// <summary>
        /// Инициализирует логарифмическое распределение Гаусса.
        /// </summary>
        public LogGaussian() { }
        /// <summary>
        /// Инициализирует логарифмическое распределение Гаусса.
        /// </summary>
        /// <param name="sigma">Среднеквадратическое отклонение</param>
        /// <param name="mu">Коэффициент сдвига</param>
        public LogGaussian(double sigma, double mu)
        {
            Sigma = sigma;
            Mu = mu;
        }
        /// <summary>
        /// Получает или задает значение среднеквадратического отклонения.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента сдивга.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return Maths.Exp(mu + sigma * sigma / 2);
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return (Maths.Exp(sigma * sigma) - 1) * Maths.Exp(2 * mu + sigma * sigma);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Maths.Exp(mu);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return Maths.Exp(mu - sigma * sigma);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Maths.Exp(sigma * sigma) + 2.0) * Maths.Sqrt(Maths.Exp(sigma * sigma) - 1.0);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return Maths.Exp(4 * sigma * sigma) + 2.0 * Maths.Exp(3 * sigma * sigma) + 3.0 * Maths.Exp(3 * sigma * sigma) - 6.0;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Exp(Maths.Pow((Maths.Log(x) - mu), 2) / (-2.0 * sigma * sigma)) / (Maths.Sqrt(2.0 * Maths.Pi) * sigma * x);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 0.5 + 0.5 * Special.Erf((Maths.Log(x) - mu) / Maths.Sqrt(sigma * 1.414));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return 0.5 + 0.5 * Maths.Log(2 * Maths.Pi * sigma * sigma) + mu;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new LogGaussian(sigma, mu);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public LogGaussian Clone()
        {
            return new LogGaussian(sigma, mu);
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
            info.AddValue("Sigma", sigma);
            info.AddValue("Mu", mu);
        }
        #endregion
    }
    /// <summary>
    /// Определяет полукруговое распределение Винера.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    /// </remarks>
    /// </summary>
    public class Wigner : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double r;
        #endregion;

        #region Wigner components
        /// <summary>
        /// Инициализирует полукруговое распределение Винера.
        /// </summary>
        /// <param name="r">Радиус</param>
        public Wigner(double r)
        {
            R = r;
        }
        /// <summary>
        /// Получает или задает значение радиуса.
        /// </summary>
        public double R
        {
            get
            {
                return this.r;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.r = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(-r, r);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return r * r / 4.0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return -1;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = 2 / (Math.PI * r2);
            return b * a;
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (Math.Abs(x) > r)
            {
                return double.NaN;
            }

            double r2 = r * r, x2 = x * x;
            double a = Math.Sqrt(r2 - x2);
            double b = x / (Math.PI * r2);
            double c = Math.Asin(x / r) / Maths.Pi;
            return 0.5 + b * a + c;
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(Maths.Pi * r) - 1.0 / 2;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Wigner(this.r);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Wigner Clone()
        {
            return new Wigner(this.r);
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
            info.AddValue("R", this.r);
        }
        #endregion
    }
    /// <summary>
    /// Определяет логарифмическое распределение Рэлея.
    /// <remarks>
    /// Распределение вероятностей случайной величины, введенное впервые в 1880 г. Джоном Уильямом Стреттом (лордом Рэлеем) 
    /// в связи с задачей сложения гармонических колебаний со случайными фазами.
    /// Носитель x ∈ [0, +inf), параметры: σ > 0 - параметр масштаба.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Rayleigh_distribution
    /// </remarks>
    /// </summary>
    public class Rayleigh : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double sigma = 1;
        #endregion;

        #region Rayleigh components
        /// <summary>
        /// Инициализирует логарифмическое распределение Рэлея.
        /// </summary>
        public Rayleigh() { }
        /// <summary>
        /// Инициализирует логарифмическое распределение Рэлея.
        /// </summary>
        /// <param name="sigma">Параметр мастаба</param>
        public Rayleigh(double sigma)
        {
            Sigma = sigma;
        }
        /// <summary>
        /// Получает или задает значение параметра мастаба.
        /// </summary>
        public double Sigma
        {
            get
            {
                return this.sigma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return Maths.Sqrt(Maths.Pi / 2.0) * sigma;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return (2.0 - Maths.Pi / 2.0) * Maths.Pow(sigma);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return sigma * Maths.Sqrt(Maths.Log(4));
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return sigma;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 * Maths.Sqrt(Maths.Pi) * (Maths.Pi - 3) / Maths.Pow(4 - Maths.Pi, 1.5);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return -(6 * Maths.Pi - 24 * Maths.Pi + 16) / Maths.Pow(4 - Maths.Pi);
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return x / sigma / sigma * Maths.Exp(-(x * x) / (2 * sigma * sigma));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-(x * x) / (2 * sigma * sigma));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return 1 + Maths.Log(sigma / Maths.Log(2)) + Maths.Gamma / 2;
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Rayleigh(sigma);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Rayleigh Clone()
        {
            return new Rayleigh(sigma);
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
            info.AddValue("Sigma", sigma);
        }
        #endregion
    }
    /// <summary>
    /// Определяет экспоненциальное распределение.
    /// <remarks>
    /// Экспоненциальное распределение — абсолютно непрерывное распределение, моделирующее время между двумя последовательными 
    /// свершениями одного и того же события.
    /// Носитель x ∈ [0, +inf), параметры: λ > 0 - интенсивность.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Exponential_distribution
    /// </remarks>
    /// </summary>
    public class Exponential : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double l = 1;
        #endregion

        #region Exp components
        /// <summary>
        /// Инициализирует экспоненциальное распределение.
        /// </summary>
        public Exponential() { }
        /// <summary>
        /// Инициализирует экспоненциальное распределение.
        /// </summary>
        /// <param name="lambda">Параметр интенсивности (0, +inf)</param>
        public Exponential(double lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Получает или задает значение параметра интенсивности (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.l = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return Maths.Pow(l, -1);
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return Maths.Pow(l, -2);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0.0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Maths.Log(2) / l;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2.0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return l * Maths.Exp(-l * x);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-l * x);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return 1 - Maths.Log(l);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Exponential(this.l);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Exponential Clone()
        {
            return new Exponential(this.l);
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
            info.AddValue("Lambda", this.l);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Коши.
    /// <remarks>
    /// Класс абсолютно непрерывных распределений. 
    /// Случайная величина, имеющая распределение Коши, является стандартным примером величины, не имеющей математического ожидания и дисперсии.
    /// Носитель x ∈ (-inf, +inf), параметры: γ > 0 - коэффициент масштаба, x0 - коэффициент сдвига.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Cauchy_distribution
    /// </remarks>
    /// </summary>
    public class Cauchy : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double g = 0.5;
        private double x0 = 0;
        #endregion

        #region Caushi components
        /// <summary>
        /// Инициализирует распределение Коши.
        /// </summary>
        public Cauchy() { }
        /// <summary>
        /// Инициализирует распределение Коши.
        /// </summary>
        /// <param name="gamma">Коэффициент масштаба (0, +inf)</param>
        /// <param name="x0">Коэффициент сдвига</param>
        public Cauchy(double gamma, double x0)
        {
            Gamma = gamma;
            X0 = x0;
        }
        /// <summary>
        /// Получает или задает значение коэффициента масштаба (0, +inf).
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.g;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.g = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента сдвига.
        /// </summary>
        public double X0
        {
            get
            {
                return this.x0;
            }
            set
            {
                this.x0 = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return double.PositiveInfinity;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return x0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return x0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return 1.0 / (Maths.Pi * g * (1.0 + Maths.Pow((x - x0) / g)));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            return 1.0 / Maths.Pi * Maths.Atg((x - x0) / g) + 0.5;
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * g);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Cauchy(g, x0);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Cauchy Clone()
        {
            return new Cauchy(g, x0);
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
            info.AddValue("Gamma", g);
            info.AddValue("X0", x0);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Вейбулла.
    /// <remarks>
    /// Распределение Вейбулла в теории вероятностей — двухпараметрическое семейство абсолютно непрерывных распределений. 
    /// Названо в честь Валодди Вейбулла, детально охарактеризовавшего его в 1951, хотя впервые его определил Фреше в 1927, 
    /// а применено оно было ещё в 1933 для описания распределения размеров частиц.
    /// Носитель x ∈ [0, +inf), параметры: λ > 0 - коэффициент масштаба, k - коэффициент формы.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Weibull_distribution
    /// </remarks>
    /// </summary>
    public class Weibull : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double l = 1;
        private double k = 1;
        #endregion

        #region Weibull components
        /// <summary>
        /// Инициализирует распределение Вейбулла.
        /// </summary>
        public Weibull() { }
        /// <summary>
        /// Инициализирует распределение Вейбулла.
        /// </summary>
        /// <param name="lambda">Коэффициент масштаба (0, +inf)</param>
        /// <param name="k">Коэффициент формы (0, +inf)</param>
        public Weibull(double lambda, double k)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Получает или задает значение коэффициента масштаба (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.l = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента формы (0, +inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                this.k = Maths.Max(0.0000001, value);
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return l * Special.Gamma(1.0 + 1.0 / k);
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return l * l * (Special.Gamma(1.0 + 2.0 / k) - Special.Gamma(1.0 + 1.0 / k));
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return l * Maths.Pow(k - 1, 1.0 / k) / Maths.Pow(k, 1.0 / k);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return l * Maths.Pow(Maths.Log(2), 1.0 / k);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (Special.Gamma(1.0 + 3.0 / k) * Maths.Pow(l, 3) - 3 * Mean * Special.Gamma(1 + 2.0 / k) * Maths.Pow(l, 2) + 2 * Maths.Pow(Mean, 3)) / (Variance * Maths.Sqrt(Variance));
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return (Special.Gamma(1.0 + 4.0 / k) * Maths.Pow(l, 4) - 4 * Mean * Special.Gamma(1 + 3.0 / k) * Maths.Pow(l, 3) + 6 * Maths.Pow(Mean, 2) * Math.Pow(l, 2) * Special.Gamma(1.0 + 2.0 / k) - 3 * Maths.Pow(Mean, 4)) / (Variance * Variance);
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return (k / l) * Maths.Pow(x / l, k - 1) * Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Exp(-Maths.Pow(x / l, k));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Maths.Gamma * (1.0 - 1.0 / k) + Maths.Pow(l / k, k) + Maths.Log(l / k);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Weibull(l, k);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Weibull Clone()
        {
            return new Weibull(l, k);
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
            info.AddValue("Lambda", l);
            info.AddValue("K", k);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Лапласа.
    /// <remarks>
    /// В теории вероятностей и статистике распределение Лапласа является непрерывным распределением вероятностей, 
    /// названным по имени Пьера-Симона Лапласа. Его также иногда называют двойным экспоненциальным распределением, 
    /// потому что его можно рассматривать как два экспоненциальных распределения (с дополнительным параметром местоположения), 
    /// соединенных вместе друг с другом, хотя термин «двойное экспоненциальное распределение» также иногда используется для обозначения 
    /// Распределение Гумбеля.
    /// Носитель x ∈ (-inf, +inf), параметры: a > 0 - коэффициент масштаба, b - коэффициент сдвига.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Laplace_distribution
    /// </remarks>
    /// </summary>
    public class Laplace : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a = 1;
        private double b = 0;
        #endregion

        #region Laplace components
        /// <summary>
        /// Инициализирует распределение Лапласа.
        /// </summary>
        public Laplace() { }
        /// <summary>
        /// Инициализирует распределение Лапласа.
        /// </summary>
        /// <param name="alfa">Коэффициент масштаба (0, +inf)</param>
        /// <param name="beta">Коэффициент сдвига</param>
        public Laplace(double alfa, double beta)
        {
            Alfa = alfa;
            Beta = beta;
        }
        /// <summary>
        /// Получает или задает значение коэффициента масштаба (0, +inf).
        /// </summary>
        public double Alfa
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.a = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента сдвига.
        /// </summary>
        public double Beta
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
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return 2.0 / a / a;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return b;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 3;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return a / 2.0 * Maths.Exp(-a * Maths.Abs(x - b));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= b)
            {
                return 0.5 * Maths.Exp(a * (x - b));
            }
            return 1 - 0.5 * Maths.Exp(-a * (x - b));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Maths.Log(4 * Maths.Pi * a);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Laplace(a, b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Laplace Clone()
        {
            return new Laplace(a, b);
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
            info.AddValue("Alfa", a);
            info.AddValue("Beta", b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет геометрическое распределение.
    /// <remarks>
    /// Распределение дискретной случайной величины равной количеству испытаний случайного эксперимента до наблюдения первого «успеха».
    /// Носитель x ∈ [0, +inf), параметры: p ∈ [0, 1] - вероятность успеха.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Geometric_distribution
    /// </remarks>
    /// </summary>
    public class Geometric : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double p = 0.2;
        private double q = 0.8;
        #endregion

        #region Geometric components
        /// <summary>
        /// Инициализирует геометрическое распределение.
        /// </summary>
        public Geometric() { }
        /// <summary>
        /// Инициализирует геометрическое распределение.
        /// </summary>
        /// <param name="p">Вероятность "успеха" [0, 1]</param>
        public Geometric(double p)
        {
            P = p;
        }
        /// <summary>
        /// Получает или задает значение вероятности "успеха" [0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value < 0 || value > 1)
                    throw new ArgumentException("Неверное значение аргумента");

                this.p = value;
                this.q = 1 - this.p;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return q / p;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return q / p / p;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (2 - p) / Maths.Sqrt(1 - p);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6 + p * p / q;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Maths.Pow(q, x) * p;
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return 1 - Maths.Pow(q, x);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return -Maths.Log2(p) - q / p * Maths.Log2(q);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Geometric(this.p);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Geometric Clone()
        {
            return new Geometric(this.p);
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
            info.AddValue("P", p);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Пуассона.
    /// <remarks>
    /// Вероятностное распределение дискретного типа, моделирует случайную величину, представляющую собой число событий, произошедших за фиксированное
    /// время, при условии, что данные события происходят с некоторой фиксированной средней интенсивностью и независимо друг от друга.
    /// Распределение Пуассона играет ключевую роль в теории массового обслуживания.
    /// Носитель x ∈ [0, +inf), параметры: λ > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Poisson_distribution
    /// </remarks>
    /// </summary>
    public class Poisson : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double l = 1;
        #endregion

        #region Poisson components
        /// <summary>
        /// Инициализирует распределение Пуассона.
        /// </summary>
        public Poisson() { }
        /// <summary>
        /// Инициализирует распределение Пуассона.
        /// </summary>
        /// <param name="lambda">Параметр λ (0, +inf)</param>
        public Poisson(double lambda)
        {
            Lambda = lambda;
        }
        /// <summary>
        /// Получает или задает значение параметра λ (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.l;
            }
            set
            {
                if (value < 0)
                    throw new ArgumentException("Неверное значение аргумента");

                this.l = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return l;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Floor(l);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(l + 0.333 - 0.02 / l);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Pow(l, -0.5);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return Math.Pow(l, -1.0);
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Exp(-l) * Math.Pow(l, x) / Special.Factorial(x);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaQ(x + 1, l);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return l * (1 - Math.Log(l)) + Math.Exp(-l) * Row(l);
            }
        }
        /// <summary>
        /// Возвращает значение суммы ряда энтропии.
        /// </summary>
        /// <param name="l">Лямба</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private double Row(double l)
        {
            double sum = 0;
            int k, n = 20;
            double fac;

            for (k = 0; k < n; k++)
            {
                fac = Special.Factorial(k);
                sum += Math.Pow(l, k) * Math.Log(fac) / fac;
            }

            return sum;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Poisson(this.l);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Poisson Clone()
        {
            return new Poisson(this.l);
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
            info.AddValue("Lambda", this.l);
        }
        #endregion
    }
    /// <summary>
    /// Определяет равномерное распределение.
    /// <remarks>
    /// В теории вероятностей и статистике непрерывное равномерное распределение или прямоугольное распределение представляет собой семейство симметричных 
    /// вероятностных распределений, так что для каждого члена семейства все интервалы одинаковой длины на носителе распределения равновероятны. Опора 
    /// определяется двумя параметрами: a и b, которые являются его минимальными и максимальными значениями. Распределение часто сокращается U (a, b). 
    /// Это максимальное распределение вероятности энтропии для случайной переменной X при отсутствии каких-либо ограничений, кроме того, что она содержится 
    /// в поддержке распределения.
    /// Носитель x ∈ (-inf, +inf), параметры: a, b ∈ (-inf, inf).
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
    /// </remarks>
    /// </summary>
    public class Uniform : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Uniform components
        /// <summary>
        /// Инициализирует равномерное распределение.
        /// </summary>
        public Uniform() { }
        /// <summary>
        /// Инициализирует равномерное распределение.
        /// </summary>
        /// <param name="a">Параметр сдвига a</param>
        /// <param name="b">Параметр сдвига b</param>
        public Uniform(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Получает или задает параметр сдвига a.
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
        /// Получает или задает параметр сдвига b.
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
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return (a + b) / 2.0;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return Math.Pow(b - a, 2)  / 12.0;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return (a - b) / 2.0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Mean;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return -1.2;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 0;
            }
            return 1.0 / (b - a);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < a)
            {
                return 0;
            }
            else if (x > b)
            {
                return 1;
            }
            return (x - a) / (b - a);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return Math.Log(b - a);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Uniform(this.a, this.b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Uniform Clone()
        {
            return new Uniform(this.a, this.b);
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
            info.AddValue("A", this.a);
            info.AddValue("B", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет бета-распределение.
    /// <remarks>
    /// В теории вероятностей и статистике бета-распределение представляет собой семейство непрерывных вероятностных распределений, определенных на 
    /// интервале [0, 1], параметризованных двумя положительными параметрами формы, обозначенными как α и β, которые фигурируют как показатели случайной 
    /// величины и управляют формой распределения.
    /// Носитель x ∈ [0, 1], параметры: a > 0, b > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Beta_distribution
    /// </remarks>
    /// </summary>
    public class Beta : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region Beta components
        /// <summary>
        /// Инициализирует бета-распределение.
        /// </summary>
        public Beta() { }
        /// <summary>
        /// Инициализирует бета-распределение.
        /// </summary>
        /// <param name="a">Параметр a</param>
        /// <param name="b">Параметр b</param>
        public Beta(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Получает или задает параметр a.
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
        /// Получает или задает параметр b.
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
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, 1);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return a / (a + b);
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return (a * b) / Maths.Pow(a + b) / (a + b + 1);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (a > 1 && b > 1)
                {
                    return (a - 1) / (a + b - 2);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 * (b - a) * Math.Sqrt(a + b + 1) / (a + b + 2) / Math.Sqrt(a * b);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                double a2 = a * a, b2 = b * b, a3 = a2 * a;
                return 6 * (a3 - a2 * (2 * b - 1) + b2 * (b + 1) - 2 * a * b * (b + 2)) / (a * b * (a + b + 2) * (a + b + 3));
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, a - 1) * Math.Pow(1 - x, b - 1) / Special.Beta(a, b);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x > 1)
            {
                return 0;
            }
            else if (x < 0)
            {
                return 0;
            }
            return Special.Beta(a, b, x);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Beta(this.a, this.b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Beta Clone()
        {
            return new Beta(this.a, this.b);
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
            info.AddValue("A", this.a);
            info.AddValue("B", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет Г-распределение.
    /// <remarks>
    /// Г-распределение является двухпараметрическим семейством непрерывных вероятностных распределений, определенных на интервале (0, +inf), параметризованных 
    /// двумя положительными параметрами формы, обозначенными как θ и k, которые фигурируют как показатели случайной величины и управляют формой распределения.
    /// Носитель x ∈ [0, +inf), параметры: k > 0, θ > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gamma_distribution
    /// </remarks>
    /// </summary>
    public class Gamma : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double thetta = 1;
        private double k = 1;
        #endregion

        #region Gamma components
        /// <summary>
        /// Инициализирует Г-распределение.
        /// </summary>
        public Gamma() { }
        /// <summary>
        /// Инициализирует Г-распределение.
        /// </summary>
        /// <param name="thetta">Параметр θ (0, +inf)</param>
        /// <param name="k">Параметр k (0, +inf)</param>
        public Gamma(double thetta, double k)
        {
            Thetta = thetta; K = k;
        }
        /// <summary>
        /// Получает или задает параметр θ (0, +inf).
        /// </summary>
        public double Thetta
        {
            get
            {
                return this.thetta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.thetta = value;
            }
        }
        /// <summary>
        /// Получает или задает параметр k (0, +inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return thetta * k;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return k * thetta * thetta;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 1)
                {
                    return (k - 1) * thetta;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Math.Pow(x, k - 1) * Math.Exp(-x / thetta) / (Special.Gamma(k) * Math.Pow(thetta, k));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            return Special.GammaP(k, x / thetta);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return k * thetta + (1 - k) * Math.Log(thetta) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Gamma(this.thetta, this.k);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Gamma Clone()
        {
            return new Gamma(this.thetta, this.k);
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
            info.AddValue("Thetta", this.thetta);
            info.AddValue("K", this.k);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Парето.
    /// <remarks>
    /// Двухпараметрическое семейство абсолютно непрерывных распределений, являющихся степенными. Называется по имени Вилфредо Парето. Встречается при исследовании 
    /// различных явлений, в частности, социальных, экономических, физических и других. Вне области экономики иногда называется также распределением Брэдфорда.
    /// Носитель x ∈ [0, +inf), параметры: Xm > 0 - коэффицимент масштаба, k > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Pareto_distribution
    /// </remarks>
    /// </summary>
    public class Pareto : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double xm = 1;
        private double k = 1;
        #endregion

        #region Pareto components
        /// <summary>
        /// Инициализирует распределение Парето.
        /// </summary>
        public Pareto() { }
        /// <summary>
        /// Инициализирует распределение Парето.
        /// </summary>
        /// <param name="xm">Коэффициент масштаба θ (0, +inf)</param>
        /// <param name="k">Параметр k (0, +inf)</param>
        public Pareto(double xm, double k)
        {
            Xm = xm; K = k;
        }
        /// <summary>
        /// Получает или задает коэффициент масштаба Xm (0, +inf).
        /// </summary>
        public double Xm
        {
            get
            {
                return this.xm;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.xm = value;
            }
        }
        /// <summary>
        /// Получает или задает параметр k (0, +inf).
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return (k * xm) / (k - 1);
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                if (k > 2)
                {
                    return Maths.Pow(xm / k - 1) * (k / (k - 2));
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return xm;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return xm * Maths.Sqrt(2, k);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (k > 3)
                {
                    return 2 * (1 + k) / (k - 3) * Math.Sqrt((k - 2) / k);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                double k2 = k * k;
                double k3 = k2 * k;
                return 6 * (k3 + k2 + 6 * k - 2) / (k * (k - 3) * (k - 4));
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return k * Math.Pow(xm, k) / Math.Pow(x, k + 1);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < xm)
            {
                return 0;
            }
            return 1.0 - Math.Pow(xm / x, k);
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return k * xm + (1 - k) * Math.Log(xm) + Special.GammaLog(k); // + (1 - k) * Special.Ksi(k);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Pareto(this.xm, this.k);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Pareto Clone()
        {
            return new Pareto(this.xm, this.k);
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
            info.AddValue("Xm", this.xm);
            info.AddValue("K", this.k);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Бернулли.
    /// <remarks>
    /// Дискретное распределение вероятностей, моделирующее случайный эксперимент произвольной природы, при заранее известной вероятности успеха или неудачи.
    /// Носитель x ∈ {0, 1}, параметры: p ∈ [0, 1].
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Bernoulli_distribution
    /// </remarks>
    /// </summary>
    public class Bernoulli : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double p = 0.2;
        private double q = 0.8;
        #endregion

        #region Bernoully components
        /// <summary>
        /// Инициализирует распределение Бернулли.
        /// </summary>
        public Bernoulli() { }
        /// <summary>
        /// Инициализирует распределение Бернулли.
        /// </summary>
        /// <param name="p">Вероятность "успеха" [0, 1]</param>
        public Bernoulli(double p)
        {
            P = p;
        }
        /// <summary>
        /// Получает или задает значение вероятности "успеха" [0, 1].
        /// </summary>
        public double P
        {
            get
            {
                return this.p;
            }
            set
            {
                if (value < 0 || value > 1)
                    throw new ArgumentException("Неверное значение аргумента");

                this.p = value;
                this.q = 1 - this.p;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, 1);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return p;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return p * q;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                if (q > p)
                {
                    return 0;
                }
                else if (q < p)
                {
                    return 1;
                }
                return 0.5;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (p - q) / Maths.Sqrt(p * q);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return (6 * p * p - 6 * p + 1) / p * (1 - p);
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x == 0)
            {
                return q;
            }
            return p;
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }
            else if (x >= 1)
            {
                return 1;
            }
            return q;
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return -q * Maths.Log(q) - p * Maths.Log(p);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Bernoulli(this.P);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Bernoulli Clone()
        {
            return new Bernoulli(this.P);
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
            info.AddValue("P", this.P);
        }
        #endregion
    }
    /// <summary>
    /// Определяет лог-логистическое распределение.
    /// <remarks>
    /// Распределение вероятности случайной величины, логарифм которой имеет логистическое распределение. По форме оно похоже на логарифмически нормальное распределение, 
    /// но имеет более тяжелые хвосты. В отличие от логарифмически нормального распределения, его кумулятивная функция распределения может быть записана в закрытом виде.
    /// Носитель x ∈ [0, +inf], параметры: a > 0, b > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Log-logistic_distribution
    /// </remarks>
    /// </summary>
    public class LogLogistic : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a = 1;
        private double b = 1;
        #endregion

        #region LogLogistic components
        /// <summary>
        /// Инициализирует лог-логистическое распределение.
        /// </summary>
        public LogLogistic() { }
        /// <summary>
        /// Инициализирует лог-логистическое распределение.
        /// </summary>
        /// <param name="a">Параметр a</param>
        /// <param name="b">Параметр b</param>
        public LogLogistic(double a, double b)
        {
            A = a; B = b;
        }
        /// <summary>
        /// Получает или задает значение параметра a.
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.a = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра b.
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
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (b > 1)
                {
                    return a * Math.Pow((b - 1) / (b + 1), 1 / b);
                }
                return 0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return (b / a) * Math.Pow(x / a, b - 1) / (1.0 + Math.Pow(Math.Pow(x / a, b), 2));
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            return 1.0 / (1 + Math.Pow(x / a, -b));
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new LogLogistic(this.a, this.b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public LogLogistic Clone()
        {
            return new LogLogistic(this.a, this.b);
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
            info.AddValue("A", this.a);
            info.AddValue("B", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет биномиальное распределение.
    /// <remarks>
    /// Распределение количества «успехов» в последовательности из n независимых случайных экспериментов, таких, что вероятность «успеха» в каждом из них постоянна и равна p.
    /// Носитель x ∈ [0, n], параметры: n >= 0, p ∈ [0, 1].
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Binomial_distribution
    /// </remarks>
    /// </summary>
    public class Binomial : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double n = 20;
        private double p = 0.5;
        private double q = 0.5;
        #endregion

        #region Binomial components
        /// <summary>
        /// Инициализирует биномиальное распределение.
        /// </summary>
        public Binomial() { }
        /// <summary>
        /// Инициализирует биномиальное распределение.
        /// </summary>
        /// <param name="n">Число испытаний (>0)</param>
        /// <param name="p">Вероятность успеха [0, 1]</param>
        public Binomial(double n, double p)
        {
            N = n; P = p;
        }
        /// <summary>
        /// Число испытаний.
        /// </summary>
        public double N
        {
            get 
            { 
                return n; 
            }
            set
            {
                this.n = value;
            }
        }
        /// <summary>
        /// Вероятность успеха [0, 1].
        /// </summary>
        public double P
        {
            get 
            { 
                return p; 
            }
            set
            {
                if (value > 1 || value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.p = value;
                this.q = 1.0 - p;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, this.n);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return n * p; 
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return n * p * q;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                double test = (n + 1) * p;

                if (test <= 0 || (int)test != test)
                    return Math.Floor(test);

                if (test <= n)
                    return test;

                if (test == n + 1)
                    return n;

                return Double.NaN;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Floor(n * p);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return (q - p) / Math.Sqrt(n * p * q);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return (1 - 6 * p * q) / Variance;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0 || x > n)
            {
                return 0;
            }

            double a = Special.LogBinomial(n, x);
            double b = x == 0 ? 0 : x * Math.Log(p);
            double c = (n - x);
            double d = Math.Log(1 - p);
            double log = a + b + c * d;

            return Math.Exp(log);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
                return 0;
            if (x >= n)
                return 1;

            double a = n - x;
            double b = x + 1;
            return Special.Beta(a, b, q);
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Binomial(this.n, this.p);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Binomial Clone()
        {
            return new Binomial(this.n, this.p);
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
            info.AddValue("N", this.n);
            info.AddValue("P", this.p);
        }
        #endregion
    }
    /// <summary>
    /// Определяет гипергеометрическое распределение.
    /// <remarks>
    /// Гипергеометрическое распределение в теории вероятностей моделирует количество удачных выборок без возвращения из конечной совокупности.
    /// Носитель x ∈ [0, K], параметры: N ∈ [0, +inf], D ∈ [0, N], K ∈ [0, N].
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Hypergeometric_distribution
    /// </remarks>
    /// </summary>
    public class Hypergeometric : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double n = 30; // эквивалент N
        private double k = 20; // эквивалент k
        private double d = 20; // эквивалент D
        #endregion

        #region Hypergeometric components
        /// <summary>
        /// Инициализирует гипергеометрическое распределение.
        /// </summary>
        public Hypergeometric() { }
        /// <summary>
        /// Инициализирует гипергеометрическое распределение.
        /// </summary>
        /// <param name="n">Параметр N [0, +inf]</param>
        /// <param name="k">Параметр D [0, N]</param>
        /// <param name="d">Параметр K [0, N]</param>
        public Hypergeometric(double n, double k, double d) 
        {
            N = n; K = k; D = d;
        }
        /// <summary>
        /// Получает или задает значение параметра N [0, +inf].
        /// </summary>
        public double N
        {
            get
            {
                return n;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.n = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра D [0, N].
        /// </summary>
        public double D
        {
            get
            {
                return d;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Неверное значение аргумента");

                this.d = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра k [0, N].
        /// </summary>
        public double K
        {
            get
            {
                return k;
            }
            set
            {
                if (value < 0 || value > N)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, k);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return k * d / n; 
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return k * (d / n) * ((n - d) / n) * ((n - k) / (n - 1.0)); }
        }
        /// <summary>
        /// Получает значени моды.
        /// </summary>
        public double Mode
        {
            get
            {
                double num = (k + 1) * (d + 1);
                double den = n + 2;
                return Math.Floor(num / den);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (n - 2 * d) * Math.Pow(n - 1, 0.5) * (n - 2 * k);
                double k2 = (n * d * (n - d) * Math.Pow(n - k, 0.5) * (n - 2));
                return k1 / k2;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                double n2 = n * n;
                double k1 = (n2 * (n - 1)) / (k * (n - 2) * (n - 3) * (n - k));
                double k2 = (n * (n + 1) - 6 * n * (n - k)) / (d * (n - d));
                double k3 = (3 * n * (n - k) * (n + 6)) / n2 - 6;
                return k1 * (k2 + k3);
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < Math.Max(0, k + d - n) || x > Math.Min(d, k))
            {
                return 0;
            }

            double a = Special.Binomial(d, x);
            double b = Special.Binomial(n - d, k - x);
            double c = Special.Binomial(n, k);
            return (a * b) / c;
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double sum = 0;
            int k = (int)x;
            for (int i = 0; i <= k; i++)
            {
                sum += Function(i);
            }

            return sum;
        }
        /// <summary>
        /// Получает значение энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Hypergeometric(n, k, d);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Hypergeometric Clone()
        {
            return new Hypergeometric(n, k, d);
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
            info.AddValue("N", this.n);
            info.AddValue("K", this.k);
            info.AddValue("D", this.d);
        }
        #endregion
    }
    /// <summary>
    /// Определяет логистическое распределение.
    /// <remarks>
    /// Один из видов абсолютно непрерывных распределений. Формой напоминает нормальное распределение, но имеет более «тяжёлые» концы и больший коэффициент 
    /// эксцесса.
    /// Носитель x ∈ (-inf, +inf), параметры: μ ∈ (-inf, +inf), s ∈ (0, +inf).
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Logistic_distribution
    /// </remarks>
    /// </summary>
    public class Logistic : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double mu = 5;
        private double s = 2;
        #endregion

        #region Logistic components
        /// <summary>
        /// Инициализирует логистическое распределение.
        /// </summary>
        /// <param name="mu">Параметр μ</param>
        /// <param name="s">Параметра s (0, +inf]</param>
        public Logistic(double mu, double s)
        {
            Mu = mu; S = s;
        }
        /// <summary>
        /// Инициализирует логистическое распределение.
        /// </summary>
        public Logistic() { }
        /// <summary>
        /// Получает или задает значение параметра μ.
        /// </summary>
        public double Mu
        { 
            get 
            { 
                return mu; 
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра s (0, +inf].
        /// </summary>
        public double S
        {
            get 
            {
                return s; 
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.s = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {   
            get 
            { 
                return mu; 
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return (s * s * Math.PI * Math.PI) / 3.0; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return mu; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 1.2;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get
            {
                return Math.Log(s) + 2;
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / s;
            return 1.0 / (1 + Math.Exp(-z));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            double z = (x - mu) / s;
            double num = Math.Exp(-z);
            double a = (1 + num);
            double den = s * a * a;

            return num / den;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Logistic(mu, s);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Logistic Clone()
        {
            return new Logistic(mu, s);
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
            info.AddValue("Mu", this.mu);
            info.AddValue("S", this.s);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Радемахера.
    /// <remarks>
    /// В теории вероятности и статистике распределение Радемахера (названное в честь Ганса Радемахера) представляет собой дискретное распределение 
    /// вероятности, при котором случайная величина x имеет 50%-ный шанс быть либо +1, либо -1.
    /// Носитель x ∈ {-1, +1}.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Rademacher_distribution
    /// </remarks>
    /// </summary>
    public class Rademacher : IDistribution, ICloneable
    {
        #region Rademacher components
        /// <summary>
        /// Инициализирует распределение Радемахера.
        /// </summary>
        public Rademacher() { }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(-1, 1);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return 1; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return -2;
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < -1)
            {
                return 0;
            }
            if (x >= 1)
            {
                return 1;
            }
            return 0.5;
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x == -1)
            {
                return 0.5;
            }
            if (x == 1)
            {
                return 0.5;
            }
            return 0.0;
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get
            {
                return Math.Log(2);
            }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Rademacher();
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Rademacher Clone()
        {
            return new Rademacher();
        }
        #endregion
    }
    /// <summary>
    /// Определяет треугольное распределение.
    /// <remarks>
    /// Треугольное распределение представляет собой непрерывное распределение вероятности с нижним пределом a, верхним пределом b и модой c.
    /// Носитель x ∈ [a, b], параметры: a, b, c.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Triangular_distribution
    /// </remarks>
    /// </summary>
    public class Triangular : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a;
        private double b;
        private double c;
        #endregion

        #region Triangular components
        /// <summary>
        /// Инициализирует треугольное распределение.
        /// </summary>
        public Triangular() { }
        /// <summary>
        /// Инициализирует треугольное распределение.
        /// </summary>
        /// <param name="a">Параметр a ∈ (-inf, +inf)</param>
        /// <param name="b">Параметр b ∈ (-inf, +inf)</param>
        /// <param name="c">Параметр c ∈ (-inf, +inf)</param>
        public Triangular(double a, double b, double c)
        {
            A = a; B = b; C = c;
        }
        /// <summary>
        /// Получает или задает значение параметра a ∈ (-inf, +inf).
        /// </summary>
        public double A
        { 
            get 
            { 
                return a; 
            }
            set
            {
                this.a = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра b ∈ (-inf, +inf).
        /// </summary>
        public double B 
        { 
            get 
            { 
                return b; 
            }
            set
            {
                this.b = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра c ∈ (-inf, +inf).
        /// </summary>
        public double C
        {
            get
            {
                return c;
            }
            set
            {
                this.c = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(a, b);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return (a + b + c) / 3.0; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return (a * a + b * b + c * c - a * b - a * c - b * c) / 18; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                double median;
                if (c >= (a + b) / 2.0)
                {
                    median = a + Math.Sqrt((b - a) * (c - a)) / 1.4142135623730950488016887242097;
                }
                else
                {
                    median = b - Math.Sqrt((b - a) * (b - c)) / 1.4142135623730950488016887242097;
                }

                return median;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return c; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                double k1 = (a + b - 2 * c) * (2 * a - b - c) * (a - 2 * b + c);
                double k2 = 5 * (a * a + b * b + c * c - a * b - a * c - b * c);
                return 1.4142135623730950488016887242097 * k1 / Math.Pow(k2, 1.5);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return -0.6;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get
            {
                return 0.5 + Math.Log((b - a) / 2);
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return ((x - a) * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return 1 - ((b - x) * (b - x)) / ((b - a) * (b - c));

            return 1;
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < a)
                return 0;

            if (x >= a && x <= c)
                return (2 * (x - a)) / ((b - a) * (c - a));

            if (x > c && x <= b)
                return (2 * (b - x)) / ((b - a) * (b - c));

            return 0;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Triangular(a, b, c);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Triangular Clone()
        {
            return new Triangular(a, b, c);
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
            info.AddValue("A", this.a);
            info.AddValue("B", this.b);
            info.AddValue("C", this.c);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Накагами.
    /// <remarks>
    /// Распределение Накагами или распределение Накагами-м является распределением вероятности, связанным с гамма-распределением. Он имеет два 
    /// параметра: параметр формы μ ≥ 0.5 и второй параметр, управляющий разбросом Ω > 0.
    /// Носитель x ∈ [0, +inf), параметры: μ ≥ 0.5, Ω > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Nakagami_distribution
    /// </remarks>
    /// </summary>
    public class Nakagami : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double mu;
        private double omega;
        private double constant; // 2 * μ ^ μ / (Γ(μ) * ω ^ μ))
        private double nratio;   // -μ / ω
        private double twoMu1;   // 2 * μ - 1.0
        #endregion

        #region Nakagami components
        /// <summary>
        /// Инициализирует распределение Накагами.
        /// </summary>
        public Nakagami()
        {
            Initialize(0.5, 1);
        }
        /// <summary>
        /// Инициализирует распределение Накагами.
        /// </summary>
        /// <param name="mu">Коэффициент формы</param>
        /// <param name="omega">Коэффициент распространения</param>
        public Nakagami(double mu, double omega)
        {
            Initialize(mu, omega);
        }
        /// <summary>
        /// Инициализация распределения.
        /// </summary>
        /// <param name="mu">Коэффициент формы</param>
        /// <param name="omega">Коэффициент распространения</param>
        private void Initialize(double mu, double omega)
        {
            Mu = mu;
            Omega = omega;

            double twoMuMu = 2.0 * Math.Pow(mu, mu);
            double gammaMu = Special.Gamma(mu);
            double spreadMu = Math.Pow(omega, mu);
            nratio = -mu / omega;
            twoMu1 = 2.0 * mu - 1.0;

            constant = twoMuMu / (gammaMu * spreadMu);
        }
        /// <summary>
        /// Получает или задает значение коэффициента формы.
        /// </summary>
        public double Mu
        {
            get 
            { 
                return mu; 
            }
            set
            {
                if (value <= 0.5)
                    throw new Exception("Неверное значение аргумента");

                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает значение коэффициента распространения.
        /// </summary>
        public double Omega
        {
            get 
            { 
                return omega; 
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.omega = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return (Special.Gamma(mu + 0.5) / Special.Gamma(mu)) * Math.Sqrt(omega / mu);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                double a = Math.Sqrt(2) / 2;
                double b = ((2 * mu - 1) * omega) / mu;
                return a * Math.Sqrt(b);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                double a = Special.Gamma(mu + 0.5) / Special.Gamma(mu);
                return omega * (1.0 - (1.0 / mu) * (a * a));
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        ///  Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return Special.GammaP(mu, (mu / omega) * (x * x));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }

            return constant * Math.Pow(x, twoMu1) * Math.Exp(nratio * x * x);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Nakagami(mu, omega);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Nakagami Clone()
        {
            return new Nakagami(mu, omega);
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
            info.AddValue("Mu", this.mu);
            info.AddValue("Omega", this.omega);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Леви.
    /// <remarks>
    /// Распределение Леви, названное в честь Пола Леви, является непрерывным распределением вероятностей для неотрицательной случайной величины. 
    /// В спектроскопии это распределение с частотой как зависимая переменная известно как профиль Ван-дер-Ваальса. Это частный случай обратного 
    /// гамма-распределения.
    /// Носитель x ∈ [μ, +inf), параметры: μ ∈ (-inf, +inf), c > 0 - коэффициент масштаба.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
    /// </remarks>
    /// </summary>
    public class Levy : IDistribution, ICloneable, ISerializable
    {
        #region Prviate data
        private double mu = 0;
        private double scale = 1;
        #endregion

        #region Levy components
        /// <summary>
        /// Инициализирует распределение Леви.
        /// </summary>
        /// <param name="mu">Коэффициент сдвига μ</param>
        /// <param name="c">Коэффициент масштаба (>0)</param>
        public Levy(double mu, double c)
        {
            Mu = mu; C = c;
        }
        /// <summary>
        /// Получает или задает коэффициент сдвига.
        /// </summary>
        public double Mu
        {
            get 
            { 
                return mu; 
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает коэффициент масштаба (>0).
        /// </summary>
        public double C
        {
            get 
            { 
                return scale; 
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.scale = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return Double.PositiveInfinity; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return Double.PositiveInfinity; }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (mu == 0)
                {
                    return scale / 3.0;
                }
                return Double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get
            {
                return (1.0 + 3.0 * Maths.Gamma + Math.Log(16 * Math.PI * scale * scale)) / 2.0;
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < mu)
            {
                return 0;
            }

            return Special.Erfc(Math.Sqrt(scale / (2 * (x - mu))));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < mu)
            {
                return 0;
            }
            double z = x - mu;
            double a = Math.Sqrt(scale / (2.0 * Math.PI));
            double b = Math.Exp(-(scale / (2 * z)));
            double c = Math.Pow(z, 3.0 / 2.0);

            return a * b / c;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Levy(mu, scale);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Levy Clone()
        {
            return new Levy(mu, scale);
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
            info.AddValue("Mu", this.mu);
            info.AddValue("C", this.scale);
        }
        #endregion
    }
    /// <summary>
    /// Определяет логарифмическое распределение.
    /// <remarks>
    /// Класс дискретных распределений. Логарифмическое распределение используется в различных приложениях, включая математическую генетику и физику.
    /// Носитель x ∈ [1, +inf), параметры: p ∈ (0, 1].
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Logarithmic_distribution
    /// </remarks>
    /// </summary>
    public class Logarithmic : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double p = 0.66;
        #endregion

        #region Logarithmic components
        /// <summary>
        /// Инициализирует логарифмическое распределение.
        /// </summary>
        /// <param name="p">Параметр</param>
        public Logarithmic(double p)
        {
            P = p;
        }
        /// <summary>
        /// Получает или задает значение параметра p ∈ (0, 1].
        /// </summary>
        public double P
        {
            get 
            { 
                return p; 
            }
            set 
            {
                if (value <= 0 || value > 1)
                    throw new Exception("Неверное значение аргумента");

                this.p = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(1, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return (-1 / Math.Log(1 - p)) * p / (1 - p); 
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get 
            {
                double k1 = p + Math.Log(1 - p);
                double k2 = Math.Pow(1 - p, 2) * Math.Pow(Math.Log(1 - p), 2);
                return -p * k1 / k2;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return 1; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            if (x > 1)
            {
                return 0;
            }
            return 1 + Special.Beta(x + 1, 0) / Math.Log(1 - p);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            if (x > 1)
            {
                return 0;
            }
            return -1 / Math.Log(1 - p) * Math.Pow(p, x) / x;
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Logarithmic(p);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Logarithmic Clone()
        {
            return new Logarithmic(p);
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
            info.AddValue("P", this.p);
        }
        #endregion
    }
    /// <summary>
    /// Определяет бета распределение второго рода.
    /// <remarks>
    /// В теории вероятностей и статистике бета-распределение (также известное как инвертированное бета-распределение или бета-распределение второго рода)
    /// представляет собой абсолютно непрерывное распределение вероятностей, определенное для x ∈ (0, +inf) и двух коэффициентов формы: α > 0, β > 0.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Beta_prime_distribution
    /// </remarks>
    /// </summary>
    public class BetaPrime : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double alpha = 1; // shape (α)
        private double beta = 1;  // shape (β)
        #endregion

        #region Beta-prime components
        /// <summary>
        /// Инициализирует бета распределение второго рода.
        /// </summary>
        /// <param name="alpha">Параметр α (0, +inf)</param>
        /// <param name="beta">Параметр β (0, +inf)</param>
        public BetaPrime(double alpha, double beta)
        {
            Alpha = alpha; Beta = beta;
        }
        /// <summary>
        /// Получает или задает значение параметра α ∈ (0, +inf).
        /// </summary>
        public double Alpha
        {
            get 
            { 
                return alpha; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.alpha = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get 
            { 
                return beta; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.beta = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                if (beta > 1)
                {
                    return alpha / (beta - 1);
                }

                return Double.PositiveInfinity;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (alpha >= 1)
                {
                    return (alpha - 1) / (beta + 1);
                }

                return 0.0;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                if (beta > 2.0)
                {
                    double num = alpha * (alpha + beta - 1);
                    double den = (beta - 2) * Math.Pow(beta - 1, 2);
                    return num / den;
                }
                else if (beta > 1.0)
                {
                    return Double.PositiveInfinity;
                }

                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            return Special.BetaIncomplete(alpha, beta, x / (1 + x));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 0;
            }

            double num = Math.Pow(x, alpha - 1) * Math.Pow(1 + x, -alpha - beta);
            double den = Special.Beta(alpha, beta);
            return num / den;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new BetaPrime(alpha, beta);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public BetaPrime Clone()
        {
            return new BetaPrime(alpha, beta);
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
            info.AddValue("α", this.alpha);
            info.AddValue("β", this.beta);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Барнибаума-Сондерса.
    /// <remarks>
    /// Распределение Барнибаума-Сондерса представляет собой распределение вероятности, широко используемое в приложениях надежности для моделирования 
    /// времен сбоев. В литературе имеется несколько альтернативных формулировок этого распределения. Оно названо в честь З. У. Бирнбаума и С. С. Сондерса.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Birnbaum–Saunders_distribution
    /// </remarks>
    /// </summary>
    public class BirnbaumSaunders : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        private double gamma = 1;
        #endregion

        #region Birnbaum-Saunders components
        /// <summary>
        /// Инициализирует распределение Барнибаума-Сондерса.
        /// </summary>
        /// <param name="mu">Коэффициент сдвига μ ∈ (0, +inf)</param>
        /// <param name="beta">Коэффициент масштаба β ∈ (0, +inf).</param>
        /// <param name="gamma">Коэффициент формы γ ∈ (0, +inf)</param>
        public BirnbaumSaunders(double mu, double beta, double gamma)
        {
            Mu = mu; Beta = beta; Gamma = gamma;
        }
        /// <summary>
        /// Получает или задает коэффициент сдвига μ ∈ (0, +inf).
        /// </summary>
        public double Mu
        {
            get 
            { 
                return mu; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает коэффициент масштаба β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get 
            { 
                return beta; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.beta = value;
            }
        }
        /// <summary>
        /// Получает или задает коэффициент формы γ ∈ (0, +inf).
        /// </summary>
        public double Gamma
        {
            get 
            { 
                return gamma; 
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support 
        {
            get
            {
                return new RangeDouble(this.mu, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get 
            { 
                return 1 + 0.5 * gamma * gamma; 
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return gamma * gamma * (1 + (5 * gamma * gamma) / 4); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double a = Math.Sqrt(x);
            double b = Math.Sqrt(1.0 / x);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            double c = x - mu;

            double a = Math.Sqrt(c / beta);
            double b = Math.Sqrt(beta / c);

            double alpha = (a + b) / (2 * gamma * c);
            double z = (a - b) / gamma;

            // Normal cumulative distribution function
            return alpha * Special.Erfc(-z / 1.4142135623731) * 0.5;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new BirnbaumSaunders(mu, beta, gamma);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public BirnbaumSaunders Clone()
        {
            return new BirnbaumSaunders(mu, beta, gamma);
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
            info.AddValue("μ", this.mu);
            info.AddValue("β", this.beta);
            info.AddValue("γ", this.gamma);
        }
        #endregion
    }
    /// <summary>
    /// Определяет кси-квадрат распределение.
    /// <remarks>
    /// Распределение xи-квадрат - это распределение суммы квадратов k независимых стандартных нормальных случайных величин.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Chi-squared_distribution
    /// </remarks>
    /// </summary>
    public class ChiSquare : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private int k = 1;
        #endregion

        #region Chi-square components
        /// <summary>
        /// Инициализирует кси-квадрат распределение.
        /// </summary>
        /// <param name="k">Число степеней свободы (0, +inf)</param>
        public ChiSquare(int k)
        {
            K = k;
        }
        /// <summary>
        /// Получает или задает число степеней свободы (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return k;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return k - 2.0 / 3;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (k >= 2)
                {
                    return k - 2;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return 2 * k;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return Math.Sqrt(8.0 / k);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / k;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get
            {
                double k2 = k / 2.0;
                double s1 = Math.Log(2.0 * Special.Gamma(k2));
                double s2 = (1.0 - k2) * Special.DiGamma(k2);
                return k2 + s1 + s2;
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return 0;
            }
            return Special.GammaP(k / 2.0, x / 2.0);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return 1;
            }
            return Special.GammaQ(k / 2.0, x / 2.0);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new ChiSquare(k);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public ChiSquare Clone()
        {
            return new ChiSquare(k);
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
            info.AddValue("K", this.k);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение распределение Гамбела.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gumbel_distribution
    /// </remarks>
    /// </summary>
    public class Gumbel : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double mu = 0;
        private double beta = 1;
        #endregion

        #region Gumbel components
        /// <summary>
        /// Инициализирует распределение Гамбела.
        /// </summary>
        /// <param name="mu">Коэффициент сдвига μ ∈ (-inf, +inf)</param>
        /// <param name="beta">Коэффициент масштаба β ∈ (0, +inf).</param>
        public Gumbel(double mu, double beta)
        {
            Mu = mu; Beta = beta;
        }
        /// <summary>
        /// Получает или задает коэффициент сдвига μ ∈ (-inf, +inf).
        /// </summary>
        public double Mu
        {
            get
            {
                return mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает коэффициент масштаба β ∈ (0, +inf).
        /// </summary>
        public double Beta
        {
            get
            {
                return beta;
            }
            set
            {
                if (value < 0)
                    throw new Exception("Неверное значение аргумента");

                this.beta = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return mu + beta * Maths.Gamma;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return mu - beta * Math.Log(Math.Log(2));
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return ((Math.PI * Math.PI) / 6.0) * beta * beta; }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return mu;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 1.14;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 12.0 / 5;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(beta) + Maths.Gamma + 1; }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double z = (x - mu) / beta;
            return Math.Exp(-Math.Exp(-z));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            double z = (x - mu) / beta;
            return (1 / beta) * Math.Exp(-(z + Math.Exp(-z)));
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Gumbel(mu, beta);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Gumbel Clone()
        {
            return new Gumbel(mu, beta);
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
            info.AddValue("μ", this.mu);
            info.AddValue("β", this.beta);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Стьюдента.
    /// <remarks>
    /// Распределение Стьюдента (t-распределение) в теории вероятностей — это однопараметрическое семейство абсолютно непрерывных распределений. 
    /// Названо в честь Уильяма Сили Госсета, который первым опубликовал работы, посвящённые этому распределению, под псевдонимом «Стьюдент».
    /// Распределение Стьюдента играет важную роль в некоторых широко используемых системах статистического анализа.
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Student%27s_t-distribution
    /// </remarks>
    /// </summary>
    public class Student : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double lambda;
        private double degrees;
        #endregion

        #region Student components
        /// <summary>
        /// Инициализирует распределение Стьюдента.
        /// </summary>
        /// <param name="n">Число степеней свободы n ∈ (0, +inf)</param>
        public Student(double n)
        {
            this.N = n;
            double num = Special.GammaLog((n + 1) / 2.0);
            double den = 0.5 * Math.Log(n * Maths.Pi) + Special.GammaLog(n / 2.0);
            this.lambda = num - den;
        }
        /// <summary>
        /// Получает или задает число степеней свободы n ∈ (0, +inf).
        /// </summary>
        public double N
        {
            get
            {
                return this.degrees;
            }
            set
            {
                if (value < 1)
                    throw new Exception("Неверное значение аргумента");

                this.degrees = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                if (degrees > 2)
                    return degrees / (degrees - 2);
                else if (degrees > 1)
                    return Double.PositiveInfinity;
                return Double.NaN;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (degrees > 3)
                {
                    return 0;
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                if (degrees > 4)
                {
                    return 6.0 / (degrees - 4);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get 
            {
                double a = Special.DiGamma((1 + degrees) / 2.0);
                double b = Special.DiGamma(degrees / 2.0);
                double c = (degrees + 1) / 2.0 * (a - b);
                double d = Math.Sqrt(degrees) * Special.Beta(degrees / 2.0, 0.5);

                return c + Math.Log(d);
            }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double v = degrees;
            double sqrt = Math.Sqrt(x * x + v);
            double u = (x + sqrt) / (2 * sqrt);
            return Special.BetaIncomplete(v / 2.0, v / 2.0, u);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return Math.Exp(LogFunction(x));
        }
        /// <summary>
        /// Логарифмическая плотность вероятности.
        /// </summary>
        /// <param name="x">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        private double LogFunction(double x)
        {
            return lambda - ((degrees + 1) / 2.0) * Math.Log((x * x) / degrees);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Student(degrees);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Student Clone()
        {
            return new Student(degrees);
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
            info.AddValue("n", this.degrees);
        }
        #endregion
    }
    /// <summary>
    /// Определяет U-квадратическое распределение.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/U-quadratic_distribution
    /// </remarks>
    /// </summary>
    public class UQuadratic : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        double a;
        double b;
        double alpha;
        double beta;
        #endregion

        #region UQuadratic components
        /// <summary>
        /// Инициализирует U-квадратическое распределение.
        /// </summary>
        /// <param name="a">Параметр a ∈ (0, +inf)</param>
        /// <param name="b">Параметр b ∈ (a, +inf)</param>
        public UQuadratic(double a, double b)
        {
            A = a; B = b;
            this.alpha = 12 / Math.Pow(b - a, 3);
            this.beta = (b + a) / 2;
        }
        /// <summary>
        /// Получает или задает значение параметра a ∈ (0, +inf).
        /// </summary>
        public double A
        {
            get
            {
                return this.a;
            }
            set
            {
                if (a < 0)
                    throw new Exception("Неверное значение аргумента");

                this.a = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра b ∈ (a, +inf).
        /// </summary>
        public double B
        {
            get
            {
                return this.b;
            }
            set
            {
                if (b < a)
                    throw new Exception("Значение параметра b должно быть либо больше, либо равно a");

                this.b = value;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return (int)a & (int)b;
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { return (a + b) / 2.0d; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return (3.0d / 20.0) * Math.Pow(b - a, 2.0); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 3.0 / 112.0 * Math.Pow(b - a, 4.0);
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(a, b); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 1;

            return (alpha / 3) * (Math.Pow(x - beta, 3) + Math.Pow(beta - a, 3));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < a)
                return 0;

            if (x > b)
                return 0;

            return alpha * Math.Pow(x - beta, 2);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new UQuadratic(a, b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public UQuadratic Clone()
        {
            return new UQuadratic(a, b);
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
            info.AddValue("a", this.a);
            info.AddValue("b", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Фишера.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/F-distribution
    /// </remarks>
    /// </summary>
    public class FisherSnedecor : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        // Distribution parameters
        private int d1;
        private int d2;

        // Derived values
        private double b;
        private double? mean;
        private double? variance;
        #endregion

        #region Fisher-Snedecor components
        /// <summary>
        /// Инициализирует распределение Фишера.
        /// </summary>
        /// <param name="d1">Первая степень свободы</param>
        /// <param name="d2">Вторая степень свободы</param>
        public FisherSnedecor(int d1 = 1, int d2 = 1)
        {
            if (d1 <= 0)
                throw new ArgumentOutOfRangeException("d1", "Значение должно быть больше нуля.");
            if (d2 <= 0)
                throw new ArgumentOutOfRangeException("d2", "Значение должно быть больше нуля.");

            this.d1 = d1;
            this.d2 = d2;

            this.b = Special.Beta(d1 * 0.5, d2 * 0.5);
        }
        /// <summary>
        /// Получает значение первой степени свободы.
        /// </summary>
        public int D1
        {
            get { return d1; }
        }
        /// <summary>
        /// Получает значение второй степени свободы.
        /// </summary>
        public int D2
        {
            get { return d2; }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                if (!mean.HasValue)
                {
                    if (d2 <= 2)
                    {
                        mean = Double.NaN;
                    }
                    else
                    {
                        mean = d2 / (d2 - 2.0);
                    }
                }

                return mean.Value;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                if (!variance.HasValue)
                {
                    if (d2 <= 4)
                    {
                        variance = Double.NaN;
                    }
                    else
                    {
                        variance = (2.0 * d2 * d2 * (d1 + d2 - 2)) /
                            (d1 * (d2 - 2) * (d2 - 2) * (d2 - 4));
                    }
                }

                return variance.Value;
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (d1 > 2)
                {
                    double a = (d1 - 2.0) / d1;
                    double b = d2 / (d2 + 2.0);
                    return a * b;
                }

                return Double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                if (d2 > 6)
                {
                    double v1 = 2 * d1 + d2 - 2;
                    double v2 = Math.Sqrt(8 * (d2 - 4));
                    double v3 = Math.Sqrt(d1 * (d1 + d2 - 2));
                    double v4 = d2 - 6;
                    return v1 * v2 / (v3 * v4);
                }
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, double.PositiveInfinity); }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
                return 0;

            double u = (d1 * x) / (d1 * x + d2);
            return Special.BetaIncomplete(d1 * 0.5, d2 * 0.5, u);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
                return 0;

            double u = Math.Pow(d1 * x, d1) * Math.Pow(d2, d2) /
                Math.Pow(d1 * x + d2, d1 + d2);
            return Math.Sqrt(u) / (x * b);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new FisherSnedecor(d1, d2);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public FisherSnedecor Clone()
        {
            return new FisherSnedecor(d1, d2);
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
            info.AddValue("d1", this.d1);
            info.AddValue("d2", this.d2);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Эрланга.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Erlang_distribution
    /// </remarks>
    /// </summary>
    public class Erlang : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private int k = 1;
        private double lambda = 0.5;
        #endregion

        #region Erlang distribution
        /// <summary>
        /// Инициализирует распределение Эрланга.
        /// </summary>
        /// <param name="k">Параметр формы k ∈ (0, +inf)</param>
        /// <param name="lambda">λ-параметр λ ∈ (0, +inf)</param>
        public Erlang(int k, double lambda)
        {
            K = k; Lambda = lambda;
        }
        /// <summary>
        /// Получает или задает значение параметра k ∈ (0, +inf).
        /// </summary>
        public int K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра λ ∈ (0, +inf).
        /// </summary>
        public double Lambda
        {
            get
            {
                return this.lambda;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.lambda = value;
            }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(0, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                return this.k / this.lambda;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return double.NaN;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                return this.k / Math.Pow(this.lambda, 2.0);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return (this.k - 1) / this.lambda;
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get
            {
                return 2.0 / Math.Sqrt(k);
            }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get
            {
                return 6.0 / k;
            }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Entropy
        {
            get
            {
                return (1 - k) * Special.DiGamma(k) + Math.Log(Special.Gamma(k) / lambda) + k;
            }
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Math.Pow(lambda, k) * Math.Pow(x, k - 1) * Math.Exp(-lambda * x) / Special.Factorial(k - 1);
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return double.NaN;
            }
            return Special.GammaIncomplete(k, lambda * x) / Special.Factorial(k - 1);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Erlang(this.k, this.lambda);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Erlang Clone()
        {
            return new Erlang(this.k, this.lambda);
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
            info.AddValue("k", this.k);
            info.AddValue("λ", this.lambda);
        }
        #endregion
    }
    /// <summary>
    /// Определяет компактное распределение Коши.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Wrapped_Cauchy_distribution
    /// </remarks>
    /// </summary>
    public class WrappedCauchy : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double mu;
        private double gamma;
        #endregion

        #region Wrapped distribution
        /// <summary>
        /// Инициалазирует компактное распределение Коши.
        /// </summary>
        /// <param name="mu">Параметр μ</param>
        /// <param name="gamma">Параметр γ > 0</param>
        public WrappedCauchy(double mu, double gamma)
        {
            this.mu = mu;
            this.gamma = gamma;
        }
        /// <summary>
        /// Получает или задает значение параметра μ.
        /// </summary>
        public double Mu
        {
            get
            {
                return this.mu;
            }
            set
            {
                this.mu = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра γ > 0.
        /// </summary>
        public double Gamma
        {
            get
            {
                return this.gamma;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.gamma = value;
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return mu; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return 1 - Math.Exp(-gamma); }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(-Math.PI, Math.PI); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(2 * Math.PI * (1 - Math.Exp(-2 * gamma))); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            double constant = (1.0 / (2 * Math.PI));
            return constant * Math.Sinh(gamma) / (Math.Cosh(gamma) - Math.Cos(x - mu));
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new WrappedCauchy(mu, gamma);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public WrappedCauchy Clone()
        {
            return new WrappedCauchy(mu, gamma);
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
            info.AddValue("μ", this.mu);
            info.AddValue("γ", this.gamma);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Кумарасвы.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Kumaraswamy_distribution
    /// </remarks>
    /// </summary>
    public class Kumaraswamy : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a;
        private double b;
        #endregion

        #region Distribution
        /// <summary>
        /// Инициализирует распределение Кумарасвы.
        /// </summary>
        /// <param name="a">Параметр формы распределения a > 0</param>
        /// <param name="b">Параметр формы распределения b > 0</param>
        public Kumaraswamy(double a, double b)
        {
            this.A = a;
            this.B = b;
        }
        /// <summary>
        /// Получает или задает параметр формы распределения a > 0.
        /// </summary>
        public double A
        {
            get
            {
                return a;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.a = value;
            }
        }
        /// <summary>
        /// Получает или задает параметр формы распределения b > 0.
        /// </summary>
        public double B
        {
            get
            {
                return b;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.b = value;
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get
            {
                double num = b * Special.Gamma(1 + (1 / a)) * Special.Gamma(b);
                double den = Special.Gamma(1 + (1 / a) + b);

                return num / den;
            }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get
            {
                double alpha = momentGeneratingFunction(2, a, b);
                double beta = Math.Pow(momentGeneratingFunction(1, a, b), 2);
                return alpha - beta;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(1 - Math.Pow(2, -1 / b), 1 / a);
            }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {

                if ((a >= 1) && (b >= 1) && (a != 1 && b != 1))
                {
                    double num = a - 1;
                    double den = a * b - 1;
                    return Math.Pow(num / den, 1 / a);
                }

                return Double.NaN;
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x > 1)
                return 1;

            if (x < 0)
                return 0;

            double xa = Math.Pow(x, a);
            return 1 - Math.Pow(1 - xa, b);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x > 1)
                return 0;

            if (x < 0)
                return 0;

            return a * b * Math.Pow(x, a - 1) * Math.Pow(1 - Math.Pow(x, a), b - 1);
        }
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        protected static double momentGeneratingFunction(int n, double a, double b)
        {
            return (b * Special.Beta(1.0 + ((double)n) / a, b));
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Kumaraswamy(a, b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Kumaraswamy Clone()
        {
            return new Kumaraswamy(a, b);
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
            info.AddValue("a", this.a);
            info.AddValue("b", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Гомперца.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Gompertz_distribution
    /// </remarks>
    /// </summary>
    public class Gompertz : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double eta;
        private double b;
        #endregion

        #region Gompertz distribution
        /// <summary>
        /// Инициализирует распределение Гомперца.
        /// </summary>
        /// <param name="eta">Параметр формы η > 0</param>
        /// <param name="b">Параметр масштаба b > 0</param>
        public Gompertz(double eta, double b)
        {
            Eta = eta; B = b;
        }
        /// <summary>
        /// Получает или задает значение параметра формы η > 0.
        /// </summary>
        public double Eta
        {
            get
            {
                return this.eta;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.eta = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра масштаба b > 0.
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
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                if (eta >= 1)
                    return 0;

                return (1 / b) * Math.Log(1 / eta);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return (1.0 / b) * Math.Log((-1 / eta) * Math.Log(0.5) + 1);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x < 0)
            {
                return 0;
            }

            double ebx = Math.Exp(b * x);
            return 1.0 - Math.Exp(-eta * (ebx - 1.0));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 0)
            {
                return 0;
            }

            double a1 = b * eta * Math.Exp(eta);
            double a2 = Math.Exp(b * x);
            double a3 = Math.Exp(-eta * Math.Exp(b * x));
            return a1 * a2 * a3;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Gompertz(this.eta, this.b);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Gompertz Clone()
        {
            return new Gompertz(this.eta, this.b);
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
            info.AddValue("η", this.eta);
            info.AddValue("b", this.b);
        }
        #endregion
    }
    /// <summary>
    /// Определяет гиперболическое секансное распределение.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Hyperbolic_secant_distribution
    /// </remarks>
    /// </summary>
    public class HyperbolicSecant : IDistribution, ICloneable, ISerializable
    {
        #region Distribution
        /// <summary>
        /// Ицинициализирует гиперболическое секансное распределение.
        /// </summary>
        public HyperbolicSecant() { }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return 1; }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { return 2; }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { return (4.0 / Math.PI) * Maths.G; }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            double angle = Math.Atan(Math.Exp(x * Math.PI / 2.0));
            return 2 * angle / Math.PI;
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            return 0.5 * Maths.Sch(x * (Math.PI / 2.0));
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new HyperbolicSecant();
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public HyperbolicSecant Clone()
        {
            return new HyperbolicSecant();
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context) { }
        #endregion
    }
    /// <summary>
    /// Определяет арксинусное распределение.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Arcsine_distribution
    /// </remarks>
    /// </summary>
    public class Arcsine : IDistribution, ICloneable, ISerializable
    {
        #region Distribution
        /// <summary>
        /// Ицинициализирует арксинусное распределение.
        /// </summary>
        public Arcsine() { }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { return 0.5; }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { return 1.0 / 8; }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { return double.NaN; }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { return 0; }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { return -1.5; }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(0, 1); }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { return Math.Log(Math.PI / 4); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            return 2.0 / Math.PI * Math.Asin(Math.Sqrt(x));
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x < 1 && x > 0)
            {
                return 1.0 / (Math.PI * Math.Sqrt(x * (1 - x)));
            }
            return double.NaN;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Arcsine();
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Arcsine Clone()
        {
            return new Arcsine();
        }
        #endregion

        #region Serialization members
        /// <summary>
        /// Получает информацию об объекте.
        /// </summary>
        /// <param name="info">Данные, необходимые для сериализации и диссериализации объекта</param>
        /// <param name="context">Источник и назначение заданного потока</param>
        public void GetObjectData(SerializationInfo info, StreamingContext context) { }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Бюрра.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Burr_distribution
    /// </remarks>
    /// </summary>
    public class Burr : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double c;
        private double k;
        #endregion

        #region Burr distribution
        /// <summary>
        /// Инициализирует распределение Бюрра.
        /// </summary>
        /// <param name="c">Параметр формы c > 0</param>
        /// <param name="k">Параметр масштаба k > 0</param>
        public Burr(double c, double k)
        {
            C = c; K = k;
        }
        /// <summary>
        /// Получает или задает значение параметра формы c > 0.
        /// </summary>
        public double C
        {
            get
            {
                return this.c;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.c = value;
            }
        }
        /// <summary>
        /// Получает или задает значение параметра масштаба k > 0.
        /// </summary>
        public double K
        {
            get
            {
                return this.k;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.k = value;
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { return k * Special.Beta(k - 1.0 / c, 1.0 + 1.0 / c); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return Math.Pow((c - 1) / (k * c + 1), 1.0 / c);
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get
            {
                return Math.Pow(Math.Pow(2, 1.0 / k) - 1.0, 1.0 / c);
            }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(double.Epsilon, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            if (x <= 0)
            {
                return double.NaN;
            }

            return 1.0 - Math.Pow(1.0 + Math.Pow(x, c), -k);
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            if (x <= 0)
            {
                return double.NaN;
            }

            double a = c * k;
            double b = Math.Pow(x, c - 1);
            double d = 1 + Math.Pow(x, c);
            return a * b / Math.Pow(d, k + 1);
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new Burr(this.c, this.k);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public Burr Clone()
        {
            return new Burr(this.c, this.k);
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
            info.AddValue("c", this.c);
            info.AddValue("k", this.k);
        }
        #endregion
    }
    /// <summary>
    /// Определяет Z-распределение Фишера.
    /// <remarks>
    /// Более подробную информацию можно найти на сайте:
    /// https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
    /// </remarks>
    /// </summary>
    public class FisherZ : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double d1;
        private double d2;
        #endregion

        #region FisherZ distribution
        /// <summary>
        /// Инициализирует Z-распределение Фишера.
        /// </summary>
        /// <param name="d1">Степень свободы d1 > 0</param>
        /// <param name="d2">Степень свободы d2 > 0</param>
        public FisherZ(double d1, double d2)
        {
            D1 = d1; D2 = d2;
        }
        /// <summary>
        /// Получает или задает значение степени свободы d1 > 0.
        /// </summary>
        public double D1
        {
            get
            {
                return this.d1;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.d1 = value;
            }
        }
        /// <summary>
        /// Получает или задает значение степени свободы d2 > 0.
        /// </summary>
        public double D2
        {
            get
            {
                return this.d2;
            }
            set
            {
                if (value <= 0)
                    throw new Exception("Неверное значение аргумента");

                this.d2 = value;
            }
        }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get
            {
                return 0;
            }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get { return new RangeDouble(Double.NegativeInfinity, Double.PositiveInfinity); }
        }
        /// <summary>
        /// Возвращает значение функции распределения вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double x)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Возвращает значение функции плотности вероятности.
        /// </summary>
        /// <param name="x">Случайная величина</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double x)
        {
            // helpers:

            double d12 = d1 / 2.0;
            double d22 = d2 / 2.0;

            // first equation:

            double a = Math.Pow(d1, d12);
            double b = Math.Pow(d2, d22);
            double c = 2 * a * b;
            double d = Special.Beta(d12, d22);
            double e = c / d;

            // second equation:

            double f = Math.Exp(d1 * x);
            double g = d1 * Math.Exp(2 * x) + d2;
            double h = d12 + d22;
            double j = f / Math.Pow(g, h);

            // result of F(x, d1, d2):
            return e * j;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new FisherZ(this.d1, this.d2);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public FisherZ Clone()
        {
            return new FisherZ(this.d1, this.d2);
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
            info.AddValue("d1", this.d1);
            info.AddValue("d2", this.d2);
        }
        #endregion
    }
    #endregion

    #region Unsorted distributions
    /// <summary>
    /// Определяет распределение конической формы.
    /// <remarks>
    /// Более подробную информацию можно получить на сайте:
    /// https://en.wikipedia.org/wiki/Cone-shape_distribution_function
    /// </remarks>
    /// </summary>
    public class ConeShape : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a;
        #endregion

        #region Cone-Shape components
        /// <summary>
        /// Инициализирует распределение конической формы.
        /// </summary>
        /// <param name="a">Коэффициент</param>
        public ConeShape(double a = 0.001)
        {
            A = a;
        }
        /// <summary>
        /// Получает или задает значение коэффициента.
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
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции плотности ядра.
        /// </summary>
        /// <param name="eta">Носитель</param>
        /// <param name="tau">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double eta, double tau)
        {
            double ksi = Maths.Pi * eta * tau;
            double psi = Math.Exp(-2 * Maths.Pi * a * tau * tau);
            return Math.Sin(ksi) / ksi * psi;
        }
        /// <summary>
        /// Возвращает значение функции распределения ядра.
        /// </summary>
        /// <param name="t">Носитель</param>
        /// <param name="tau">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double t, double tau)
        {
            if (Math.Abs(tau) >= 2 * Math.Abs(t))
            {
                return 1.0 / tau * Math.Exp(-2 * Maths.Pi * a * tau * tau);
            }
            return 0;
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new ConeShape(this.a);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public ConeShape Clone()
        {
            return new ConeShape(this.a);
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
            info.AddValue("α", this.a);
        }
        #endregion
    }
    /// <summary>
    /// Определяет распределение Чой-Вильямса.
    /// <remarks>
    /// Более подробную информацию можно получить на сайте:
    /// https://en.wikipedia.org/wiki/Choi%E2%80%93Williams_distribution_function
    /// </remarks>
    /// </summary>
    public class ChoiWilliams : IDistribution, ICloneable, ISerializable
    {
        #region Private data
        private double a;
        #endregion

        #region Choi-Williams components
        /// <summary>
        /// Инициализирует распределение Чой-Вильямса.
        /// </summary>
        /// <param name="a">Коэффициент</param>
        public ChoiWilliams(double a = 0.001)
        {
            A = a;
        }
        /// <summary>
        /// Получает или задает значение коэффициента.
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
        /// Получает значение математического ожидания.
        /// </summary>
        public double Mean
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        public double Variance
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        public double Median
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        public double Mode
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        public double Skewness
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        public double Excess
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        public RangeDouble Support
        {
            get
            {
                return new RangeDouble(double.NegativeInfinity, double.PositiveInfinity);
            }
        }
        /// <summary>
        /// Получает значение дифференциальной энтропии.
        /// </summary>
        public double Entropy
        {
            get { throw new NotSupportedException(); }
        }
        /// <summary>
        /// Возвращает значение функции плотности ядра.
        /// </summary>
        /// <param name="eta">Носитель</param>
        /// <param name="tau">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Function(double eta, double tau)
        {
            double ksi = eta * tau;
            return Math.Exp(-a * ksi * ksi);
        }
        /// <summary>
        /// Возвращает значение функции распределения ядра.
        /// </summary>
        /// <param name="t">Носитель</param>
        /// <param name="tau">Носитель</param>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        public double Distribution(double t, double tau)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Clone members
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        object ICloneable.Clone()
        {
            return new ChoiWilliams(this.a);
        }
        /// <summary>
        /// Создает копию распределения.
        /// </summary>
        /// <returns>Распределение</returns>
        public ChoiWilliams Clone()
        {
            return new ChoiWilliams(this.a);
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
            info.AddValue("α", this.a);
        }
        #endregion
    }
    #endregion

    #region Probabilities
    /// <summary>
    /// Определяет класс вероятностей Байеса.
    /// </summary>
    public class Bayes
    {
        #region Private data
        private double[] Pp; // массив апостериорных вероятностей,
        private double Pa;   // знчение полной вероятности.
        private int N;       // длинна массивов.
        #endregion

        #region Bayes components
        /// <summary>
        /// Инициализирует класс вероятностей Байеса.
        /// </summary>
        /// <param name="stat">Массив статистических вероятностей</param>
        /// <param name="prior">Массив априорных вероятностей (до опыта)</param>
        public Bayes(double[] stat, double[] prior)
        {
            if (stat.Length != prior.Length)
                throw new Exception("Массивы должны быть одинаковых размерностей.");

            // Запись данных:
            this.N = prior.Length;
            this.Pp = new double[N];
            this.Pa = 0;
            int i;

            // Вычисление полной вероятности:
            for (i = 0; i < N; i++)
            {
                Pa += stat[i] * prior[i];
            }
            // Вычисление апостериорных вероятностей:
            for (i = 0; i < N; i++)
            {
                Pp[i] = stat[i] * prior[i] / Pa;
            }
        }
        /// <summary>
        /// Возвращает значение полной вероятности.
        /// </summary>
        public double General
        {
            get
            {
                return Pa;
            }
        }
        /// <summary>
        /// Возвращает массив значений апостериорных вероятностей (после опыта).
        /// </summary>
        public double[] Probabilities
        {
            get
            {
                return Pp;
            }
        }
        #endregion
    }
    #endregion

    #region Interface
    /// <summary>
    /// Определяет общий интерфейс распределений.
    /// </summary>
    public interface IDistribution
    {
        #region Components
        /// <summary>
        /// Получает интервал носителя функции. 
        /// </summary>
        RangeDouble Support { get; }
        /// <summary>
        /// Получает значение математического ожидания.
        /// </summary>
        double Mean { get; }
        /// <summary>
        /// Получает значение дисперсии.
        /// </summary>
        double Variance { get; }
        /// <summary>
        /// Получает значение медианы.
        /// </summary>
        double Median { get; }
        /// <summary>
        /// Получает значение моды.
        /// </summary>
        double Mode { get; }
        /// <summary>
        /// Получает значение коэффициента ассиметрии.
        /// </summary>
        double Skewness { get; }
        /// <summary>
        /// Получает значение коэффициента эксцесса.
        /// </summary>
        double Excess { get; }
        /// <summary>
        /// Возвращает значение дифференциальной энтропии.
        /// </summary>
        /// <returns>Число двойной точности с плавающей запятой</returns>
        double Entropy { get; }
        #endregion
    }
    #endregion
}
