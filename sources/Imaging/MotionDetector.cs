using System;
using SkiaDrawing;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the motion detector.
    /// <remarks>
    /// Implements IDisposable interface.
    /// </remarks>
    /// </summary>
    [Serializable]
    public class MotionDetector : IDisposable
    {
        #region Private data
        private Bitmap Frame;
        private object locker = new object();
        #endregion

        #region Class components
        /// <summary>
        /// Initializes motion detector.
        /// </summary>
        /// <param name="threshold">Threshold [0, 255]</param>
        /// <param name="useFilter">Use bitmap filter or not</param>
        public MotionDetector(byte threshold = 15, bool useFilter = false)
        {
            Threshold = threshold;
            UseFilter = useFilter;
        }
        /// <summary>
        /// Gets or sets threshold.
        /// </summary>
        public byte Threshold { get; set; }
        /// <summary>
        /// Use filter or not.
        /// </summary>
        public bool UseFilter { get; set; }
        /// <summary>
        /// Reset motion detector.
        /// </summary>
        public void Reset()
        {
            // synchronize
            lock (locker)
            {
                if (Frame != null)
                {
                    Frame.Dispose();
                    Frame = null;
                }
            }
        }
        /// <summary>
        /// Apply filter and returns motion level in range [0, 1].
        /// </summary>
        /// <param name="Data">Bitmap</param>
        /// <returns>Motion level</returns>
		public double Apply(Bitmap Data)
        {
            // synchronize
            lock (locker)
            {
                if (Frame == null)
                {
                    // create initial backgroung image
                    Frame = (Bitmap)Data.Clone();

                    // just return for the first time
                    return 0.0;
                }

                // creating clone of current frame
                var temp = (Bitmap)Data.Clone();

                // lock in memory
                var bitmapData = BitmapFormat.Lock32bpp(Data);
                var frameData = BitmapFormat.Lock32bpp(Frame);

                // calculate alarm
                var alarm = UseFilter ? ProcessFrameWithFilter(bitmapData, frameData) :
                    ProcessFrameWithoutFilter(bitmapData, frameData);

                // unlock
                Data.Unlock(bitmapData);
                Frame.Unlock(frameData);

                // update detector
                Frame.Dispose();
                Frame = temp;
                return alarm;
            }
        }
        #endregion

        #region IDisposable

        private bool _disposed;

        /// <inheritdoc/>
        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <inheritdoc/>
        protected virtual void Dispose(bool disposing)
        {
            if (!_disposed)
            {
                if (disposing)
                {
                    Reset();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~MotionDetector()
        {
            Dispose(false);
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Process frame.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe double ProcessFrameWithoutFilter(BitmapData bmData, BitmapData bmSrc)
        {
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;

            double variance = 0.0;

            for (int x = 0; x < width; x++)
            {
                int y, k;

                for (y = 0; y < height; y++, dst += 4, src += 4)
                {
                    var difference = 0.0f;

                    for (k = 0; k < 3; k++)
                    {
                        difference += Math.Abs(dst[k] - src[k]) / 3.0f;
                    }

                    var summary = Maths.Byte(difference);

                    if (summary > Threshold)
                    {
                        variance++;
                    }
                }
            };

            return variance / (width * height);
        }
        /// <summary>
        /// Process frame.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe double ProcessFrameWithFilter(BitmapData bmData, BitmapData bmSrc)
        {
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;

            double variance = 0.0;

            for (int x = 0; x < width; x++)
            {
                int y, k;

                for (y = 0; y < height; y++, dst += 4, src += 4)
                {
                    for (k = 0; k < 3; k++)
                    {
                        var difference = Math.Abs(dst[k] - src[k]);
                        dst[k] = Maths.Byte(difference);
                    }

                    var summary = RGB.Average(dst[2], dst[1], dst[0]);

                    if (summary > Threshold)
                    {
                        variance++;
                    }
                }
            };

            return variance / (width * height);
        }
        #endregion
    }
}
