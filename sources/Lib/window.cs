// UMapx.NET framework
// Digital Signal Processing Library.
// 
// Copyright © UMapx.NET, 2015-2020
// Valery Asiryan
// Moscow, Russia

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
    // Window Toolbox includes a set of tools for synthesizing and 
    // orthogonalizing window functions. It implements discrete short-time 
    // Fourier transforms and Weyl-Heisenberg transforms (Gabor analysis) for 
    // real and complex signals.
    // **************************************************************************
    // Designed by Valery Asiryan (c), 2015-2020
    // Moscow, Russia.
    // **************************************************************************

    #region Window functions
    /// <summary>
    /// Defines the window function of Planck.
    /// </summary>
    [Serializable]
    public class Planck : WindowBase
    {
        #region Private data
        private double a = 0.15;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Planck window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter [0, 0.5]</param>
        public Planck(int frameSize, double a = 0.15)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter [0, 0.5].
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
        /// Function Z+-(x, a).
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="p">Sign</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the window function of Tukey.
    /// </summary>
    [Serializable]
    public class Tukey : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Tukey window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter [0, 1]</param>
        public Tukey(int frameSize, double a = 1)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter [0, 1].
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
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the closed Gaussian window.
    /// </summary>
    [Serializable]
    public class Confined : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the closed Gaussian window.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (0.14 * N)</param>
        public Confined(int frameSize, double sigma = 1)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Initializes a Gaussian window function closed.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Confined(int frameSize)
        {
            this.Sigma = 0.14 * frameSize;
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
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
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
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
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        private double G(double x)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (2 * sigma);
            return Math.Exp(-t * t);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
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
    /// Defines a generalized window normal function.
    /// </summary>
    [Serializable]
    public class Normal : WindowBase
    {
        #region Private data
        private double sigma = 1;
        private double p = 2;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes a generalized window normal function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Standard deviation (>0)</param>
        /// <param name="pow">Power<remarks>For p = 2 - Gaussian window</remarks></param>
        public Normal(int frameSize, double sigma = 1, double pow = 2)
        {
            this.Sigma = sigma;
            this.FrameSize = frameSize;
            this.p = pow;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
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
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Power.
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
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            double a = (frameSize - 1) / 2;
            double t = (x - a) / (sigma * a);
            return Math.Exp(-Math.Pow(t, p));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
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
    /// Defines the Kaiser window function.
    /// </summary>
    [Serializable]
    public class Kaiser : WindowBase
    {
        #region Private data
        private double a = 3;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Kaiser window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="a">Form parameter</param>
        public Kaiser(int frameSize, double a = 3)
        {
            this.FrameSize = frameSize;
            this.A = a;
        }
        /// <summary>
        /// Gets or sets the value of the form parameter.
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
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Kaiser window:
            double u = 2 * x / (frameSize - 1);
            double r = 1 - u * u;
            double v = r >= 0 ? Math.Sqrt(1 - u * u) : 0;
            double z = Math.PI * this.a;
            double q = Special.I(z * v, 0);
            return q / Special.I(z    , 0);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Welch window function.
    /// </summary>
    [Serializable]
    public class Welch : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Welch window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Welch(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Welch function:
            double t = (frameSize - 1) / 2.0;
            double a = (x - t) / t;
            return 1 - a * a;
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Lanczos window function.
    /// </summary>
    [Serializable]
    public class Lanzcos : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Lanczos window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Lanzcos(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Lanczos function:
            return Special.Sinc(2 * x / (frameSize - 1) - 1, Math.PI);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1);
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Parzen window function.
    /// </summary>
    [Serializable]
    public class Parzen : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Parzen window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Parzen(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
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
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the "Flat-Top" window function.
    /// </summary>
    [Serializable]
    public class FlatTop : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the "Flat-Top" window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public FlatTop(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 1 - 1.93 * Cosine.cosinefunc(2 * x, frameSize) + 1.29 * Cosine.cosinefunc(4 * x, frameSize) - 0.388 * Cosine.cosinefunc(6 * x, frameSize) + 0.028 * Cosine.cosinefunc(8 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
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
    /// Defines the Nuttall window function.
    /// </summary>
    [Serializable]
    public class Nuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Nuttall window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Nuttall(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.355768 - 0.487396 * Cosine.cosinefunc(2 * x, frameSize) + 0.144232 * Cosine.cosinefunc(4 * x, frameSize) - 0.012604 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Blackman-Nuttall window function.
    /// </summary>
    [Serializable]
    public class BlackmanNuttall : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Blackman-Nuttall window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BlackmanNuttall(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.3635819 - 0.4891775 * Cosine.cosinefunc(2 * x, frameSize) + 0.1365995 * Cosine.cosinefunc(4 * x, frameSize) - 0.0106411 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
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
    /// Defines the Blackman-Harris window function.
    /// </summary>
    [Serializable]
    public class BlackmanHarris : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Blackman-Harris window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BlackmanHarris(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.35875 - 0.48829 * Cosine.cosinefunc(2 * x, frameSize) + 0.14128 * Cosine.cosinefunc(4 * x, frameSize) - 0.01168 * Cosine.cosinefunc(6 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Blackman window function.
    /// </summary>
    [Serializable]
    public class Blackman : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Blackman window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Blackman(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.42 - 0.5 * Cosine.cosinefunc(2 * x, frameSize) + 0.08 * Cosine.cosinefunc(4 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Barlett-Hann window function.
    /// </summary>
    [Serializable]
    public class BarlettHann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Barlett-Hann window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public BarlettHann(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Berlett-Hann function:
            double a = Math.Abs(Math.Abs(x / (frameSize - 1)) - 0.5);
            return 0.62 - 0.48 * a - 0.38 * Cosine.cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Hann window function (Hanning).
    /// </summary>
    [Serializable]
    public class Hann : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Hann window function (Hanning).
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Hann(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Pow(Sine.sinefunc(x, frameSize), 2);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the Hamming window function.
    /// </summary>
    [Serializable]
    public class Hamming : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the Hamming window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Hamming(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return 0.53836 - 0.46164 * Cosine.cosinefunc(2 * x, frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the cosine window function.
    /// </summary>
    [Serializable]
    public class Cosine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the cosine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Cosine(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Factor</returns>
        internal static double cosinefunc(double x, int frameSize)
        {
            return Math.Cos(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
    /// <summary>
    /// Defines the sine window function.
    /// </summary>
    [Serializable]
    public class Sine : WindowBase
    {
        #region Window components
        /// <summary>
        /// Initializes the sine window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        public Sine(int frameSize)
        {
            this.FrameSize = frameSize;
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = frameSize - 1;
            double[] x = Matrice.Compute(0, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion

        #region Static components
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        internal static double sinefunc(double x, int frameSize)
        {
            return Math.Sin(Math.PI * x / (frameSize - 1));
        }
        #endregion
    }
    /// <summary>
    /// Defines the Gabor window function.
    /// </summary>
    [Serializable]
    public class Gabor : WindowBase
    {
        #region Private data
        private double sigma = 1;
        #endregion

        #region Window components
        /// <summary>
        /// Initializes the Gabor window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <param name="sigma">Scale parameter</param>
        public Gabor(int frameSize, double sigma = 1)
        {
            this.FrameSize = frameSize;
            this.Sigma = sigma;
        }
        /// <summary>
        /// Gets or sets the standard deviation (>0).
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
                    throw new Exception("Invalid argument value");

                this.sigma = value;
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public override double Function(double x, int frameSize)
        {
            // Gabor window function
            double y = x / frameSize;
            double z = Math.Pow(2 * Math.PI * y / sigma, 2);
            return Math.Exp(-z);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public override double[] GetWindow(int frameSize)
        {
            double t = (frameSize - 1) / 2.0;
            double[] x = Matrice.Compute(-t, t, 1);
            return this.Function(x, frameSize);
        }
        #endregion
    }
    /// <summary>
    /// Defines the class for window functions.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Window_function
    /// </remarks>
    /// </summary>
    public abstract class WindowBase : IWindow
    {
        #region Private data
        /// <summary>
        /// Window size.
        /// </summary>
        protected int frameSize;
        #endregion

        #region Window components
        /// <summary>
        /// Gets or sets the window size.
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
                    throw new Exception("Window size must be greater than 0");

                this.frameSize = value;
            }
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        public double Function(double x)
        {
            return this.Function(x, this.frameSize);
        }
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        public double[] GetWindow()
        {
            return this.GetWindow(this.frameSize);
        }
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
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
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <returns>Array</returns>
        public double[] Function(double[] x)
        {
            return this.Function(x, this.frameSize);
        }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Double precision floating point number</returns>
        public abstract double Function(double x, int frameSize);
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        public abstract double[] GetWindow(int frameSize);
        #endregion
    }
    #endregion

    #region Short-time Fourier analysis
    /// <summary>
    /// Defines fast short-time Fourier transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastShortTimeFourierTransform : IWindowTransform, ITransform
    {
        #region Private data
        private FastFourierTransform FFT;
        private IWindow window;
        private Direction direction;
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes fast short-time Fourier transform.
        /// </summary>
        /// <param name="function">Windows function</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
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
        /// Normalized transform or not.
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
        /// Gets or sets the processing direction.
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
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-time Fourier transform
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
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
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
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
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    /// <summary>
    /// Defines short-time Fourier transform.
    /// <remarks>
    /// More information can be found on the website:
    /// https://en.wikipedia.org/wiki/Short-time_Fourier_transform
    /// </remarks>
    /// </summary>
    [Serializable]
    public class ShortTimeFourierTransform : IWindowTransform, ITransform
    {
        #region Private data
        private FourierTransform FFT;
        private IWindow window;
        private Direction direction;
        private double[] coefs;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes short-time Fourier transform.
        /// </summary>
        /// <param name="function">Windows function</param>
        /// <param name="normalized">Normalized transform or not</param>
        /// <param name="direction">Processing direction</param>
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
        /// Normalized transform or not.
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
        /// Gets or sets the processing direction.
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
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Short-time Fourier transform
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
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
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
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
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward short-time Fourier Transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward short-time Fourier Transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion
    }
    #endregion

    #region Gabor analysis
    /// <summary>
    /// Defines a group of orthogonal bases and discrete Weyl-Heisenberg transforms.
    /// <remarks>
    /// More information can be found on the website:
    /// https://elibrary.ru/item.asp?id=29767333
    /// </remarks>
    /// </summary>
    [Serializable]
    public class WeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        /// <summary>
        /// Windows function.
        /// </summary>
        private IWindow window;
        /// <summary>
        /// Number of frequency shifts.
        /// </summary>
        private int m;
        /// <summary>
        /// Processing direction.
        /// </summary>
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes a group of orthogonal bases and Weyl-Heisenberg transformations.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public WeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
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
                    throw new Exception("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
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
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg static components
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="N">Number of samples</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] WeylHeisenberg(IWindow window, int N, int M, bool orthogonalize = true)
        {
            return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.GetPacket(window, N), M, orthogonalize);
        }
        /// <summary>
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <param name="orthogonalize">Orthogonalized matrix or not</param>
        /// <returns>Matrix</returns>
        public static Complex[,] WeylHeisenberg(double[] g0, int M, bool orthogonalize = true)
        {
            if (orthogonalize)
            {
                return WeylHeisenbergTransform.WeylHeisenberg(WeylHeisenbergTransform.Zak(g0, M), M);
            }

            return WeylHeisenbergTransform.WeylHeisenberg(g0, M);
        }
        /// <summary>
        /// Returns a vector of window function values.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="length">Number of samples</param>
        /// <returns>Array</returns>
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
        /// Returns a vector of values of a window function that satisfies the N-1 symmetry condition.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="length">Number of samples</param>
        /// <returns>Array</returns>
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
        /// Returns the complex Weyl-Heisenberg basis matrix.
        /// <remarks>
        /// Matrix dimension[N, N], where N = M * L.
        /// </remarks>
        /// </summary>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Matrix</returns>
        private static Complex[,] WeylHeisenberg(double[] g0, int M)
        {
            int N = g0.Length, L = N / M;

            if (L <= 0)
                throw new Exception("Number of frequency shifts not defined correctly");

            Complex[,] G = new Complex[N, N];
            Complex c = 2 * Maths.Pi * Maths.I;
            double a = M / 2.0;

            Parallel.For(0, N, n =>   
            {
                double phase = n - a / 2.0;            
                int k, l, u, i, j;
                Complex exp, psi;

                for (k = 0; k < M; k++)
                {
                    exp = Maths.Exp(c * k / M * phase);

                    for (l = 0; l < L; l++)
                    {
                        u = l * M + k;                
                        i = Maths.Mod(n - l * M, N);            
                        j = Maths.Mod(n + M / 2 - l * M, N);     

                        psi = new Complex(
                            (g0[i] * exp).Real,                    
                            (Maths.I * g0[j] * exp).Real);              

                        G[n, u] = psi;
                    }
                }
            });

            return G;
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            int N = A.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] B = Matrice.Dot(A, U.Hermitian());
            return B;
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            int N = B.Length;
            Complex[,] U = WeylHeisenbergTransform.WeylHeisenberg(this.window, N, this.m, true);
            Complex[] A = Matrice.Dot(B, U);
            return A;
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Zak components
        private static FourierTransform DFT = new FourierTransform(false, Direction.Vertical);
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Vertical);

        /// <summary>
        /// Implements Zak-orthogonalization of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static double[] Zak(double[] v, int M)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            double[] vort = new double[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

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

            Complex sum;
            for (i = 0; i < N; i++)
            {
                sum = 0;
                for (j = 0; j < L2; j++)
                {
                    sum += Z[j, i];
                }
                vort[i] = (sum / L2).Real;
            }

            return vort;
        }
        /// <summary>
        /// Implements Zak-orthogonalization of the vector.
        /// </summary>
        /// <param name="v">Array</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex[] Zak(Complex[] v, int M)
        {
            // Fast shaping orthogonalization algorithm
            // WH functions using a discrete Zak transform.
            // V.P. Volchkov, D.A. Petrov and V.M. Asiryan.
            // http://www.conf.mirea.ru/CD2017/pdf/p4/66.pdf

            int N = v.Length;
            Complex[] vort = new Complex[N];
            int L = N / M, L2 = L * 2, i, j;
            Complex[,] G = new Complex[L2, N];
            Complex[,] Z;

            for (i = 0; i < L2; i++)
            {
                for (j = 0; j < N; j++)
                {
                    G[i, j] = v[Maths.Mod(j + M / 2 * i, N)];
                }
            }

            if (Maths.IsPower(L2, 2))
            {
                Z = FFT.Forward(G);
            }
            else
            {
                Z = DFT.Forward(G);
            }

            double w = 2 / Math.Sqrt(M);
            double even, odd, phi;
            Complex z1, z2;

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
    /// Defines fast Weyl-Heisenberg transform.
    /// <remarks>
    /// The class represents a computationally efficient implementation of one-dimensional and two-dimensional discrete orthogonal
    /// Weyl-Heisenberg transforms.
    /// More information can be found on the website:
    /// https://elibrary.ru/title_about.asp?id=58245
    /// </remarks>
    /// </summary>
    [Serializable]
    public class FastWeylHeisenbergTransform : IWindowTransform, ITransform
    {
        #region Private data
        private static FastFourierTransform FFT = new FastFourierTransform(false, Direction.Horizontal);
        private IWindow window;
        private int m;
        private Direction direction;
        #endregion

        #region Initialize
        /// <summary>
        /// Initializes fast Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="window">Windows function</param>
        /// <param name="m">Number of frequency shifts [4, N]</param>
        /// <param name="direction">Processing direction</param>
        public FastWeylHeisenbergTransform(IWindow window, int m = 8, Direction direction = Direction.Vertical)
        {
            Window = window; M = m; Direction = direction;
        }
        /// <summary>
        /// Gets or sets number of frequency shifts [4, N].
        /// <remarks>
        /// Even number.
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
                    throw new Exception("Invalid argument value");

                this.m = value;
            }
        }
        /// <summary>
        /// Gets or sets the processing direction.
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
        /// Gets or sets the window function.
        /// </summary>
        public IWindow Window
        {
            get { return window; }
            set { window = value; }
        }
        #endregion

        #region Weyl-Heisenberg Transform
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public Complex[] Forward(Complex[] A)
        {
            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, A.Length), this.m);
            return FastWeylHeisenbergTransform.WHT(A, g0, m);
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public Complex[] Backward(Complex[] B)
        {
            double[] g0 = WeylHeisenbergTransform.Zak(WeylHeisenbergTransform.GetPacket(this.window, B.Length), this.m);
            return FastWeylHeisenbergTransform.IWHT(B, g0, m);
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
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
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Array</param>
        /// <returns>Array</returns>
        public double[] Forward(double[] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Array</param>
        /// <returns>Array</returns>
        public double[] Backward(double[] B)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="A">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Forward(double[,] A)
        {
            throw new NotSupportedException();
        }
        /// <summary>
        /// Backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="B">Matrix</param>
        /// <returns>Matrix</returns>
        public double[,] Backward(double[,] B)
        {
            throw new NotSupportedException();
        }
        #endregion

        #region Public static components
        /// <summary>
        /// Fast forward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex[] WHT(Complex[] input, double[] g0, int M)
        {
            // The function implements a fast Weil-Heisenberg direct transformation algorithm,
            // stated in the following articles:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // V.M. Asiryan, V.P. Volchkov, "EFFECTIVE IMPLEMENTATION OF THE DIRECT TRANSFORMATION OF WEIL-HEISENBERG" [2].
            // The algorithm is computationally efficient for large M.

            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);

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

            Complex x, y, z, w;

            for (l = 0; l < L; l++)
            {
                for (m = 0; m < M2; m++)
                {
                    x = A0[l, m];
                    y = A0[l, m + M2];
                    z = A0[l, M2 - m].Conjugate;
                    w = A0[l, Maths.Mod(M - m, M)].Conjugate;

                    c1re = x + y + z + w;
                    c2re = x - y - z + w;

                    x = B0[l, m];
                    y = B0[l, m + M2];
                    z = B0[l, M2 - m].Conjugate;
                    w = B0[l, Maths.Mod(M - m, M)].Conjugate;

                    c1im = x + y - z - w;
                    c2im = x - y + z - w;

                    A1[l, m] = 1.0 / (2) * (c1re + Maths.I * c2re * exp[m]);
                    B1[l, m] = 1.0 / (2 * Maths.I) * (c1im + Maths.I * c2im * exp[m]);
                }
            }

            A1 = FFT.Backward(Matrice.Conjugate(A1));
            B1 = FFT.Backward(Matrice.Conjugate(B1));

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    i = l * M + 2 * k;
                    j = l * M + 2 * k + 1;

                    x = A1[l, k];
                    y = B1[l, k];

                    output[i] = x.Real + Maths.I * y.Real;
                    output[j] = x.Imag + Maths.I * y.Imag;
                }
            }

            return output;
        }
        /// <summary>
        /// Fast backward Weyl-Heisenberg transform.
        /// </summary>
        /// <param name="input">Array</param>
        /// <param name="g0">Function</param>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
        public static Complex[] IWHT(Complex[] input, double[] g0, int M)
        {
            // The function implements a fast Weil-Heisenberg direct transformation algorithm,
            // stated in the following articles:
            // A. Vahlin, "EFFICIENT ALGORITHMS FOR MODULATION AND DEMODULATION IN OFDM-SYSTEMS" [1].
            // V.M. Asiryan, V.P. Volchkov, "EFFECTIVE IMPLEMENTATION OF THE DIRECT TRANSFORMATION OF WEIL-HEISENBERG" [2].
            // The algorithm is computationally efficient for large M.

            int N = input.Length, L = N / M, M2 = M / 2, M4 = M2 / 2;
            Complex[] output = new Complex[N];
            Complex[,] A1 = new Complex[L, M];
            Complex[,] B1 = new Complex[L, M];
            Complex[] exp = FastWeylHeisenbergTransform.GetRotation(M);
            Complex s;
            int n, k, l;

            for (k = 0; k < M; k++)
            {
                for (l = 0; l < L; l++)
                {
                    s = input[k + l * M];

                    A1[l, k] = s.Real;
                    B1[l, k] = s.Imag;
                }
            }

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

            Za = Matrice.Conjugate(FFT.Backward(Za));
            Zb = Matrice.Conjugate(FFT.Backward(Zb));

            Complex a0, a1, b0, b1;
            Complex x, y, u, v;

            for (k = 0; k < M2; k++)
            {
                for (l = 0; l < L; l++)
                {
                    a0 = Za[l, k]; a1 = Za[l, Maths.Mod(M - k, M2)].Conjugate;

                    x = 1.0 / (2) * (a0 + a1);
                    y = 1.0 / (2 * Maths.I) * (a0 - a1);
                    y *= exp[k];

                    A1[l, k] = x + y;
                    A1[l, k + M2] = x - y;

                    b0 = Zb[l, k]; b1 = Zb[l, Maths.Mod(M - k, M2)].Conjugate;

                    u = 1.0 / (2) * (b0 + b1);
                    v = 1.0 / (2 * Maths.I) * (b0 - b1);
                    v *= exp[k];

                    B1[l, k] = u + v;
                    B1[l, k + M2] = u - v;
                }
            }

            for (l = 0; l < L; l++)
            {
                for (n = 0; n < N; n++)
                {
                    output[n] += A1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M, N)] - Maths.I
                               * B1[l, Maths.Mod(n - M4, M)] * g0[Maths.Mod(n - l * M + M2, N)];
                }
            }

            return output;
        }
        /// <summary>
        /// Returns an array of phase rotations.
        /// </summary>
        /// <param name="M">Number of frequency shifts</param>
        /// <returns>Array</returns>
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
    /// Defines the interface of window functions.
    /// </summary>
    public interface IWindow
    {
        #region Interface
        /// <summary>
        /// Gets or sets the window size.
        /// </summary>
        int FrameSize { get; set; }
        /// <summary>
        /// Returns the value of a window function.
        /// </summary>
        /// <param name="x">Argument</param>
        /// <returns>Double precision floating point number</returns>
        double Function(double x);
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <returns>Array</returns>
        double[] GetWindow();
        /// <summary>
        /// Returns the window function.
        /// </summary>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        double[] GetWindow(int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <param name="frameSize">Window size</param>
        /// <returns>Array</returns>
        double[] Function(double[] x, int frameSize);
        /// <summary>
        /// Returns an array of window function values.
        /// </summary>
        /// <param name="x">Array</param>
        /// <returns>Array</returns>
        double[] Function(double[] x);
        #endregion
    }
    /// <summary>
    /// Defines the general window transform interface.
    /// </summary>
    public interface IWindowTransform
    {
        #region Interface
        /// <summary>
        /// Gets or sets the window function.
        /// </summary>
        IWindow Window { get; set; }
        #endregion
    }
    #endregion
}
