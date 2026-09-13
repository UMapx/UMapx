using System;
using System.Drawing;

namespace UMapx.Imaging
{
    /// <summary>
    /// Defines the motion event detector.
    /// </summary>
    /// <remarks>
    /// Reports the end of a motion episode after consecutive quiet frames.
    /// A single frame with motion level greater than Alarm starts an episode.
    /// Input frames are not modified. Call Reset before changing streams or frame dimensions.
    /// Implements IDisposable interface.
    /// </remarks>
    [Serializable]
    public class MotionEventDetector : IDisposable
    {
        #region Private data
        private readonly object locker = new object();
        private readonly MotionDetector detector;
        private int quietCount;
        private bool motionSeen;
        #endregion

        #region Class components
        /// <summary>
        /// Initializes motion event detector.
        /// </summary>
        /// <param name="threshold">Pixel difference threshold [0, 255].</param>
        /// <param name="alarm">Motion level threshold [0, 1]; equality counts as quiet.</param>
        /// <param name="quietFrames">Number of consecutive quiet frames required to end an episode (>0).</param>
        public MotionEventDetector(byte threshold = 15, float alarm = 0.01f, int quietFrames = 3)
        {
            Alarm = alarm;
            QuietFrames = quietFrames;
            detector = new MotionDetector(threshold, false);
        }
        /// <summary>
        /// Gets or sets the pixel difference threshold [0, 255].
        /// </summary>
        public byte Threshold
        {
            get { return detector.Threshold; }
            set { detector.Threshold = value; }
        }
        /// <summary>
        /// Gets or sets the motion level threshold [0, 1]; equality counts as quiet.
        /// </summary>
        public float Alarm { get; set; }
        /// <summary>
        /// Gets or sets the positive number of consecutive quiet frames required to end an episode.
        /// </summary>
        /// <remarks>The count measures processed frames, not elapsed time.</remarks>
        public int QuietFrames { get; set; }
        /// <summary>
        /// Reset motion event detector.
        /// </summary>
        public void Reset()
        {
            // synchronize
            lock (locker)
            {
                detector.Reset();
                motionSeen = false;
                quietCount = 0;
            }
        }
        /// <summary>
        /// Processes a frame and reports whether a motion episode has just ended.
        /// </summary>
        /// <param name="frame">Bitmap.</param>
        /// <returns>True once on the QuietFrames-th quiet frame after motion; otherwise false.</returns>
        public bool Detect(Bitmap frame)
        {
            // synchronize
            lock (locker)
            {
                // calculate motion level
                float level = detector.Apply(frame);
                if (level > Alarm)
                {
                    motionSeen = true;
                    quietCount = 0;
                    return false;
                }
                if (!motionSeen) return false;
                if (++quietCount < QuietFrames) return false;
                motionSeen = false;
                quietCount = 0;
                return true;
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
                    detector.Dispose();
                }
                _disposed = true;
            }
        }

        /// <inheritdoc/>
        ~MotionEventDetector()
        {
            Dispose(false);
        }
        #endregion
    }
}
