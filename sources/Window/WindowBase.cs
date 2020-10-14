using System;

namespace UMapx.Window
{
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
}
