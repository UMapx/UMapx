using System;
using System.Drawing;
using System.Drawing.Imaging;
using UMapx.Colorspace;
using UMapx.Core;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the motion detector.
    /// </summary>
    /// <remarks>
    /// Implements IDisposable interface.
    /// </remarks>
    [Serializable]
    public class MotionDetector : IDisposable
    {
        #region Private data
        private readonly object locker = new object();
        private Bitmap Frame;
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
        /// Applies the filter and returns motion level in range [0, 1].
        /// </summary>
        /// <param name="bitmap">Bitmap</param>
        /// <returns>Motion level</returns>
		public float Apply(Bitmap bitmap)
        {
            // synchronize
            lock (locker)
            {
                if (Frame == null)
                {
                    // create initial background image
                    Frame = (Bitmap)bitmap.Clone();

                    // just return for the first time
                    return 0.0f;
                }

                // creating clone of current frame
                var temp = (Bitmap)bitmap.Clone();

                // lock in memory
                var bitmapData = BitmapFormat.Lock32bpp(bitmap);
                var frameData = BitmapFormat.Lock32bpp(Frame);

                // calculate alarm
                var alarm = UseFilter ? ProcessFrameWithFilter(bitmapData, frameData) :
                    ProcessFrameWithoutFilter(bitmapData, frameData);

                // unlock
                bitmap.Unlock(bitmapData);
                Frame.Unlock(frameData);

                // update detector
                Frame.Dispose();
                Frame = temp;
                return alarm;
            }
        }
        /// <summary>
        /// Applies the filter and returns motion level in range [0, 1].
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <returns>Motion level</returns>
        public float Apply(BitmapData bmData)
        {
            if (bmData.PixelFormat != PixelFormat.Format32bppArgb)
                throw new NotSupportedException("Only support Format32bppArgb pixelFormat");

            // synchronize
            lock (locker)
            {
                if (Frame == null)
                {
                    // create initial background image
                    Frame = BitmapFormat.ToBitmap(bmData);

                    // just return for the first time
                    return 0.0f;
                }

                // creating clone of current frame
#pragma warning disable DF0010 // Marks undisposed local variables.
                var temp = BitmapFormat.ToBitmap(bmData);
#pragma warning restore DF0010 // Marks undisposed local variables.

                // lock in memory
                var frameData = BitmapFormat.Lock32bpp(Frame);

                // calculate alarm
                var alarm = UseFilter ? ProcessFrameWithFilter(bmData, frameData) :
                    ProcessFrameWithoutFilter(bmData, frameData);

                // unlock
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
        private unsafe float ProcessFrameWithoutFilter(BitmapData bmData, BitmapData bmSrc)
        {
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;

            float variance = 0.0f;

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
        private unsafe float ProcessFrameWithFilter(BitmapData bmData, BitmapData bmSrc)
        {
            byte* dst = (byte*)bmData.Scan0.ToPointer();
            byte* src = (byte*)bmSrc.Scan0.ToPointer();
            int width = bmData.Width, height = bmData.Height;

            float variance = 0.0f;

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
