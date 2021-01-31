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
    [Serializable]
    public class MotionDetector
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
        public MotionDetector(byte threshold = 15)
        {
            Threshold = threshold;
        }
        /// <summary>
        /// Gets or sets threshold.
        /// </summary>
        public byte Threshold { get; set; }
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
                    Frame = Data;

                    // just return for the first time
                    return 0.0;
                }

                // creating clone of current frame
                var temp = (Bitmap)Data.Clone();

                // calculate alarm
                BitmapData bitmapData = BitmapConverter.Lock32bpp(Data);
                BitmapData frameData = BitmapConverter.Lock32bpp(Frame);

                var alarm = ProcessFrame(bitmapData, frameData);

                Data.Unlock(bitmapData);
                Frame.Unlock(frameData);

                Frame.Dispose();
                Frame = temp;
                return alarm;
            }
        }
        /// <summary>
        /// Process frame.
        /// </summary>
        /// <param name="bmData">Bitmap data</param>
        /// <param name="bmSrc">Bitmap data</param>
        private unsafe double ProcessFrame(BitmapData bmData, BitmapData bmSrc)
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
